import numpy as np
import argparse

def readFile(filepath):
	print 'Reading file: ', filepath 
	data = []
	with open(filepath, 'r') as fr:
		for line in fr:
			data.append([float(num)for num in line.split(',')])
	return data

def readDataset(filepath):
	data = np.array( readFile(filepath) )
	dataSet = [] # array of Instance
	index = data[:,0]
	last_i = 0
	for i, shift in enumerate(index):
		if shift == 0:
			if(last_i != i):
				dataSet.append( Instance(data[range(last_i, i),:]) )
			last_i = i
	
	dataSet.append( Instance(data[range(last_i, len(index)),:]) )
	return dataSet # return an array of Instance

class Instance:
	def __init__(self,rawData):
		#print rawData
		self.length = len(rawData)
		self.timestamp = rawData[:,-1]
		self.accel_1 = rawData[:,1]
		self.accel_2 = rawData[:,2]
		self.accel_3 = rawData[:,3]
		self.tau_1 = rawData[:,4]
		self.tau_2 = rawData[:,5]
		self.tau_3 = rawData[:,6]
	def accel_abs(self):
		accel = [0]*self.length
		for i in range(self.length):
			accel[i] = (self.accel_1[i]**2 + self.accel_2[i]**2 + self.accel_3[i]**2)**0.5
		return accel
	
	def tau_abs(self):
		tau = [0]*self.length
		for i in range(self.length):
			tau[i] = (self.tau_1[i]**2 + self.tau_2[i]**2 + self.tau_3[i]**2)**0.5
		return tau

	def get_length(self):
		return self.length
if __name__ == '__main__':
	argparser = argparse.ArgumentParser()
	argparser.add_argument('filePath', type=str, help='filePath')
	args = argparser.parse_args()
	args = vars(args)
	dataSet = readDataset(args['filePath'])
	print 'number of lines in the file :',sum([I.get_length() for I in dataSet])
