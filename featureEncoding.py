import math
import random
 
def dynamicTimeWarp(seqA, seqB, d = lambda x,y: abs(x-y), print_flag = False):
	# create the cost matrix
	numRows, numCols = len(seqA), len(seqB)
	cost = [[0 for _ in range(numCols)] for _ in range(numRows)]

	# initialize the first row and column
	cost[0][0] = d(seqA[0], seqB[0])
	for i in xrange(1, numRows):
		cost[i][0] = cost[i-1][0] + d(seqA[i], seqB[0])

	for j in xrange(1, numCols):
		cost[0][j] = cost[0][j-1] + d(seqA[0], seqB[j])
 
	# fill in the rest of the matrix
	for i in xrange(1, numRows):
		for j in xrange(1, numCols):
			choices = cost[i-1][j], cost[i][j-1], cost[i-1][j-1]
			cost[i][j] = min(choices) + d(seqA[i], seqB[j])
	if print_flag:
		for row in cost:
			for entry in row:
				print "%03d" % entry,
			print ""
	return cost[-1][-1]

def writeFeature(fileName, instances, label=[]):
	# wtite features into libsvm format
	if len(label) == 0:
		label = ['-1']*len(instances)
	elif len(label) != len(instances):
		print 'ERROR : len(label) != len(instances)'
	with open(fileName, 'w') as fw:
		for j, features in enumerate(instances):
			feature_str = ''
			for i, f in enumerate(features):
				if f != 0:
					feature_str = ' '.join( [feature_str, ':'.join([str(i),str(f)])] )
			print >> fw, ' '.join([ str(int(label[j])), feature_str ])
	return

def readFeature(fileName):
	features = []
	labels = []
	with open(fileName, 'r') as fr:
		for line in fr:
			line = line.split()
			instance = [ float(f.split(':')[-1]) for f in line[1:] ]
			features.append(instance)
			labels.append(int(line[0]))
	return features, labels

class Dictionary:
	def __init__(self, atomNum, instances, distType = 'DTW'):
		# random patch dictionary
		# atomNum had better be much larger than the number of classes
		self.atoms = random.sample(instances, atomNum)
		self.distType = distType
	def encoding(self, instances):
		encodedInstances = []
		for i in instances:
			encodedInstances.append( [self.__dist(i, a) for a in self.atoms] )
		return encodedInstances
	def getAtoms(self):
		return self.atoms
	def __dist(self, seqA, seqB):
		if self.distType == 'DTW':
			return dynamicTimeWarp(seqA, seqB)
		else:
			print 'distType ERROR!'
			return

if __name__ == '__main__':
	seqA = [0, 0, 0, 3, 6, 13, 25, 22, 7, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
	seqB = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 5, 12, 24, 23, 8 ,3, 1, 0, 0, 0, 0, 0]
	print 'seqA', seqA
	print 'seqB', seqB
	print 'cost =', dynamicTimeWarp(seqA, seqB)
	RPDictionary = Dictionary(2, [seqA, seqB])
	writeFeature('./deleteMe', RPDictionary.encoding([seqA, seqB]))
