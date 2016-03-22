## Directories ##
SDIR=.				# Directory of source files
DDIR=./data			# Directory of data sequences
RESULT=.			# Directory to output results
## Data files ##
SEQX="samp-sines1.dat"		# File name of sequence X
SEQY="samp-sines2.dat"		# File name of sequence Y
## Result files ##
SFILE="result-sampling.txt"	# Result file name for similar subsequence pairs
## Parameters
LMIN="1500"			# Subsequence length threshold: lmin
EPS="0.05"			# Distance threshold: epsilon
BAND="5000"			# Width of Sakoe-Chiba band: w
PERIOD1="278"			# Sampling period of sequence X: Tx
PERIOD2="294"			# Sampling period of sequence Y: Ty


## Detect similar subsequence pairs ##
sband=$(expr $BAND / $PERIOD2)
$SDIR/sampling $DDIR/$SEQX $DDIR/$SEQY $LMIN $EPS $sband $PERIOD1 $PERIOD2 > $RESULT/$SFILE
