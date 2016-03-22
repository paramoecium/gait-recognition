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

def pattern_mining():
    import os
    # compute cross-similarity and dicover repeating subsequence
    sensor_data = data = np.genfromtxt('./ICS_slipperData/Alice0105db.csv', dtype=float, delimiter=',', names=True)
    timestamp = sensor_data['Timestamp']
    signal = sensor_data['Axis2']
    SDIR='./crossmatch/CrossMatch'	# Directory of source files
    DDIR='./match_data'	            # Directory of data sequences
    RESULT='.'                      # Directory to output results
    LMIN = 1000		                # Subsequence length threshold: lmin
    EPS = 0.1                       # Distance threshold: epsilon
    '''
    axis1:epsilon:0.18 lmin:300
    axis2:epsilon:0.2  lmin:0.19
    '''
    BAND = 5000		                # Width of Sakoe-Chiba band: w
    #XLEN= len(signal)		        # Sequence length of X
    XLEN = 10000   
    SFILE = RESULT+'/result.txt'    # Result file name for similar subsequence pairs
    PFILE = RESULT+'/path.txt'      # Result file name for optimal warping paths
    #np.savetxt('{}/signal.dat'.format(DDIR),signal,fmt='%f')

    commands = []
    # Detect similar subsequence pairs
    commands.append('{}/crossmatch {}/signal.dat {}/signal.dat {} {} {} > {}'.format(SDIR,DDIR,DDIR,LMIN,EPS,BAND,SFILE))
    commands.append('wc -l {}'.format(SFILE))
    # Compute warping paths
    commands.append('perl {}/createfile.pl {} {}/exefile.sh {}/checkfile {}/signal.dat {}/signal.dat {} {} {} {}'.format(SDIR,SFILE,RESULT,RESULT,DDIR,DDIR,LMIN,EPS,BAND,PFILE))
    commands.append('sh {}/exefile.sh'.format(RESULT))
    commands.append('{}/path {}/checkfile {} {}'.format(SDIR,RESULT,EPS,XLEN))
    commands.append('perl {}/gnuplot.pl {} {} {}'.format(SDIR,RESULT,XLEN,XLEN))
    commands.append('gnuplot {}/load_dat'.format(RESULT))
    commands.append('rm {}/load_dat'.format(RESULT))
    commands.append('rm {}/checkfile'.format(RESULT))
    commands.append('rm {}/exefile.sh'.format(RESULT))
    commands.append('rm {}/temp.txt'.format(RESULT))
    commands.append("find -name 'Subseq*' -delete")
    commands.append('mv {}/path.txt {}/'.format(RESULT,DDIR))
    commands.append('mv {}/result.txt {}/'.format(RESULT,DDIR))
    for c in commands:
        print c
        os.system(c)

class Instance:
	def __init__(self,rawData):
		#print rawData
		self.length = len(rawData)
		self.timestamp = rawData[:,-1]
		self.accel_1 = rawData[:,1]
		self.accel_2 = rawData[:,2]
		self.accel_3 = rawData[:,3]
		self.alpha_1 = rawData[:,4]
		self.alpha_2 = rawData[:,5]
		self.alpha_3 = rawData[:,6]
	def accel_abs(self):
		accel = [0]*self.length
		for i in range(self.length):
			accel[i] = (self.accel_1[i]**2 + self.accel_2[i]**2 + self.accel_3[i]**2)**0.5
		return accel
	
	def alpha_abs(self):
		alpha = [0]*self.length
		for i in range(self.length):
			alpha[i] = (self.alpha_1[i]**2 + self.alpha_2[i]**2 + self.alpha_3[i]**2)**0.5
		return alpha

	def get_length(self):
		return self.length
if __name__ == '__main__':
	argparser = argparse.ArgumentParser()
	argparser.add_argument('filePath', type=str, help='filePath')
	args = argparser.parse_args()
	args = vars(args)
	dataSet = readDataset(args['filePath'])
	print 'number of lines in the file :',sum([I.get_length() for I in dataSet])
