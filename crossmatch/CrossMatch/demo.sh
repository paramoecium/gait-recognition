## Directories ##
SDIR=.			# Directory of source files
DDIR=./data		# Directory of data sequences
RESULT=.		# Directory to output results
## Data files ##
SEQX="sines1.dat"	# File name of sequence X
SEQY="sines2.dat"	# File name of sequence Y
## Result files ##
SFILE="result.txt"	# Result file name for similar subsequence pairs
## Parameters ##
LMIN="1500"		# Subsequence length threshold: lmin
EPS="0.05"		# Distance threshold: epsilon
BAND="5000"		# Width of Sakoe-Chiba band: w


## Detect similar subsequence pairs ##
$SDIR/crossmatch $DDIR/$SEQX $DDIR/$SEQY $LMIN $EPS $BAND > $RESULT/$SFILE

## For naive solution ##
#$SDIR/naive $DDIR/$SEQX $DDIR/$SEQY $LMIN $EPS $BAND > $RESULT/result-naive.txt
