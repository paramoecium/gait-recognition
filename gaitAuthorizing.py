from preprocessing import *
from featureEncoding import *
from svm import svm_problem, svm_parameter
from svmutil import svm_train, svm_predict, svm_save_model, svm_read_problem
import argparse
from sklearn.utils import shuffle
from sklearn import preprocessing # preprocessing.scale(X)

TRAIN_TEST_RATIO = 0.8

if __name__=='__main__':
	argparser = argparse.ArgumentParser()
	#argparser.add_argument('fileName', type=str, help='fileName')
	args = argparser.parse_args()
	args = vars(args)
	dataSet = []
	label = []
	for i, fileName in enumerate( ['./data/arthur.csv', './data/brian.csv', './data/nofar.csv', './data/shalom.csv'] ):
		tmp = readDataset(fileName) # array of Instance
		dataSet = dataSet + tmp
		label = label + [i]*len(tmp)
		dataSet, label = shuffle(dataSet, label, random_state=0)
	cutIndex = int(TRAIN_TEST_RATIO*len(dataSet))
	## use accel_abs and tau_abs as input for encoding respectively
	data_accel = [I.accel_abs() for I in dataSet]
	data_tau = [I.tau_abs() for I in dataSet]
	RPDictionary_accel = Dictionary(50, data_accel[:cutIndex])
	RPDictionary_tau = Dictionary(50, data_tau[:cutIndex])
	aggregate_feature = [ f[0]+f[1] for f in zip( RPDictionary_accel.encoding(data_accel), RPDictionary_tau.encoding(data_tau) ) ]
	aggregate_feature_scaled = preprocessing.scale(aggregate_feature)

	writeFeature('./svm_train', aggregate_feature_scaled[:cutIndex], label[:cutIndex]) 
	writeFeature('./svm_test', aggregate_feature_scaled[cutIndex:], label[cutIndex:]) 

	## SVM training
	Y_train = label[:cutIndex]
	classNum, X_train = svm_read_problem('./svm_train')
	prob = svm_problem(Y_train, X_train)
	param = svm_parameter('-t 1 -q -d 2 -c 0.1')
	model = svm_train(prob, param)
	svm_save_model('./svm_model',model)

	## SVM predicting
	X_test = aggregate_feature_scaled[cutIndex:]
	Y_test = label[cutIndex:]
	p_labels, p_acc, p_vals = svm_predict(Y_test, X_test, model)
	print p_acc
