from preprocessing import *
from featureEncoding import *
from svm import svm_problem, svm_parameter
from svmutil import svm_train, svm_predict, svm_save_model, svm_read_problem
import argparse



if __name__=='__main__':
	argparser = argparse.ArgumentParser()
	#argparser.add_argument('fileName', type=str, help='fileName')
	args = argparser.parse_args()
	args = vars(args)
	dataSet = []
	label = []
	for i, fileName in enumerate( ['./data/arthur.csv', './data/brian.csv', './data/nofar.csv', './data/shalom.csv'] ):
		tmp = readDataset(args[fileName])
		dataSet = dataSet + tmp
		label = label + [i]*len(tmp)
	## use tau_abs as input for encoding
	1D_DataSet = [I.tau_abs() for I in dataSet]
	RPDictionary = Dictionary(20, 1D_DataSet)
	for i, a in enumerate(RPDictionary.getAtoms()):
		print i,a
	writeFeature(RPDictionary.encoding('./svm_data', 1D_DataSet, label)
