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
    '''# USAGE:
    argparser = argparse.ArgumentParser()
    argparser.add_argument('filePath', type=str, help='filePath')
    args = argparser.parse_args()
    args = vars(args)
    dataSet = readDataset(args['filePath'])
    print 'number of lines in the file :',sum([I.get_length() for I in dataSet])
    '''
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

def cross_match():
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

def salience_vector(signal_sample, arg_func=np.argmax):
    '''# DEBUG:
    signal_sample = [5,2,1,1,1,4,1,1,1,6,4,3,7,9,6,8]
    print salience_vector(signal_sample, arg_func=np.argmax)
    '''
    signal = np.array(signal_sample)
    len_signal = len(signal)
    len_interval = 2
    salience_vector = np.ones(len_signal)
    while len_interval <= len_signal:
        update = np.zeros(len_signal)
        num_interval = np.ceil(float(len_signal)/len_interval).astype(int)
        for i in xrange(num_interval):
            local_argpeak = arg_func(signal[i*len_interval:min((i+1)*len_interval, len_signal)])
            if (i+1)*len_interval <= len_signal:
                update[i*len_interval+local_argpeak] = len_interval
            else:
                update[i*len_interval+local_argpeak] = len_signal%len_interval
        #print update.astype(list)
        salience_vector = np.maximum(salience_vector, update)
        len_interval += 1
    return salience_vector.astype(list)

def cycle_extraction(signal_sampled, sample_frequency):
    max_signal_salience = salience_vector(signal_sampled, np.argmax)
    min_signal_salience = salience_vector(signal_sampled, np.argmin)
    # assume frequency of human walking ~= 1.7Hz (1.5Hz~2.5Hz)
    walking_frequency_min = 0.7
    walking_frequency_max = 2.5
    l_min = sample_frequency/walking_frequency_max
    l_max = sample_frequency/walking_frequency_min
    h_min = 2*l_min
    signal = np.array(signal_sampled)
    len_signal = len(signal)
    peaks = []
    for t in xrange(len_signal):
        if max_signal_salience[t] < h_min: continue
        if len(peaks) == 0:
            peaks.append([t,max_signal_salience[t]])
            continue
        else:
            l = t - peaks[-1][0]
        if l > l_min and l < l_max:
            #print [t,signal[t]]
            peaks.append([t,max_signal_salience[t]])
        else:
            print '    violatino:',t
    cycle_lengths = []
    for i in xrange(1,len(peaks)):
        cycle_lengths.append(peaks[i][0]-peaks[i-1][0])
    peaks = np.array(peaks)
    cycles = np.array([peaks[:-1,0], peaks[1:,0], cycle_lengths]).astype(int)
    #print cycle_lengths
    print '(h_min, l_min, l_max) =', h_min, l_min, l_max
    print 'number of cycles detected:', cycles.shape[1]
    print 'cycle length statistics:', 'mean', np.mean(cycle_lengths), 'std', np.std(cycle_lengths)
    assert cycles.shape[1] > 0, 'No Cycle Detected!!'
    return cycles.T # (start, end, length)

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
    signal_sample = [5,2,1,1,1,4,1,1,1,6,4,3,7,9,6,8]
    print salience_vector(signal_sample, arg_func=np.argmax)
    pass
