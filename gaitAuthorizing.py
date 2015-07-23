from preprocessing import *
from featureEncoding import *
from svm import svm_problem, svm_parameter
from svmutil import svm_train, svm_predict, svm_save_model, svm_read_problem
import argparse
from sklearn.utils import shuffle

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
        data_train = dataSet[:cutIndex]
        data_test = dataSet[cutIndex:]
        label_train = label[:cutIndex]
        label_test = label[cutIndex:]
	## use accel_abs and tau_abs as input for encoding respectively
	data_train_accel = [I.accel_abs() for I in data_train]
	data_train_tau = [I.tau_abs() for I in data_train]
	RPDictionary_accel = Dictionary(50, data_train_accel)
	RPDictionary_tau = Dictionary(50, data_train_tau)
	aggregate_feature = [ f[0]+f[1] for f in zip( RPDictionary_accel.encoding(data_train_accel), RPDictionary_tau.encoding(data_train_tau) ) ]
	writeFeature('./svm_train', aggregate_feature, label_train) 
	
	## SVM training
	Y_train = label_train
	classNum, X_train = svm_read_problem('./svm_data')
	prob = svm_problem(Y_train, X_train)
	param = svm_parameter('-q')
	model = svm_train(prob, param)
	svm_save_model('./svm_model',model)

	## SVM predicting
	data_test_accel = [I.accel_abs() for I in data_test]
	data_test_tau = [I.tau_abs() for I in data_test]
        X_test = [ f[0]+f[1] for f in zip( RPDictionary_accel.encoding(data_test_accel), RPDictionary_tau.encoding(data_test_tau) ) ]
	Y_test = label_test
	writeFeature('./svm_train', X_test, Y_test) 
	p_labels, p_acc, p_vals = svm_predict(Y_test, X_test, model)
	print p_acc
