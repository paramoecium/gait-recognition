from preprocessing import *
from featureEncoding import *
from svm import svm_problem, svm_parameter
from svmutil import svm_train, svm_predict, svm_save_model, svm_read_problem, svm_load_model
import argparse
from sklearn.utils import shuffle
from sklearn import preprocessing # preprocessing.scale(X)
from sklearn.metrics import confusion_matrix

TRAIN_SET_RATIO = 0.8
PATCH_SIZE = 20
if __name__=='__main__':
	filelist = ['./dataset2/avon.csv', './dataset2/brian_merge.csv', './dataset2/mon_merge.csv', './dataset2/nofar_merge.csv', './dataset1/shalom.csv']

	'''
	argparser = argparse.ArgumentParser()
	#argparser.add_argument('fileName', type=str, help='fileName')
	args = argparser.parse_args()
	args = vars(args)
	dataSet = []
	label = []
	#for i, fileName in enumerate( ['./dataset1/arthur.csv', './dataset1/brian.csv', './dataset1/nofar.csv', './dataset1/shalom.csv'] ):
	for i, fileName in enumerate( filelist ):
		tmp = readDataset(fileName) # array of Instance
		dataSet = dataSet + tmp
		print 'size:', len(tmp)
		label = label + [i]*len(tmp)
		dataSet, label = shuffle(dataSet, label, random_state=0)
	cutIndex = int(TRAIN_SET_RATIO*len(dataSet))
	## use accel_abs and alpha_abs as input for encoding respectively
	print 'learning dictionary...'
	data_accel = [I.accel_abs() for I in dataSet]
	data_alpha = [I.alpha_abs() for I in dataSet]
	RPDictionary_accel = Dictionary(PATCH_SIZE, data_accel[:cutIndex])
	RPDictionary_alpha = Dictionary(PATCH_SIZE, data_alpha[:cutIndex])
	print 'feature encoding...'
	aggregate_feature = [ f[0]+f[1] for f in zip( RPDictionary_accel.encoding(data_accel), RPDictionary_alpha.encoding(data_alpha) ) ]
	#aggregate_feature = preprocessing.scale(aggregate_feature) ## scale columns independently to have zero mean and unit variance
	
	print len(aggregate_feature[:cutIndex]),len(aggregate_feature[cutIndex:])
	writeFeature('./svm_train', aggregate_feature[:cutIndex], label[:cutIndex])
	writeFeature('./svm_test', aggregate_feature[cutIndex:], label[cutIndex:])
	'''
	print 'SVM training...'
	classifiers = []
	X_train, Y_train = readFeature('./svm_train',PATCH_SIZE*2)
	print Y_train[0:10]
	for i, fileName in enumerate( filelist[:-1] ):
		Y_train_i = [int(y==i) for y in Y_train]
		print Y_train_i[0:10]
		prob = svm_problem(Y_train_i, X_train)
		#param = svm_parameter('-t 1 -q -d 2')
		param = svm_parameter('-t 2 -q')
		model = svm_train(prob, param)
		svm_save_model('{}.model'.format(i), model)

	print 'SVM predicting...'
	ballot_box = []
	X_test, Y_test = readFeature('./svm_test',PATCH_SIZE*2)
	for i, fileName in enumerate( filelist[:-1] ):
		Y_test_i = [int(y==i) for y in Y_test]
		svm_load_model('{}.model'.format(i))
		p_labels, p_acc, p_vals = svm_predict(Y_test_i, X_test, model)
		ballot_box.append(p_labels)
		print p_acc
	ballot_box = np.array(ballot_box)
	Y_predict = [(sum(ballot_box[:,i])==0) for i in range(len(Y_test))]
	Y_test = [int(y==4) for y in Y_test]
	# true\predicted
	print confusion_matrix(Y_test, Y_predict)
