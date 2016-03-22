## Directories ##
SDIR=.			# Directory of source files
DDIR=./data		# Directory of data sequences
RESULT=.		# Directory to output results
## Data files ##
SEQX="sines1.dat"	# File name of sequence X
SEQY="sines2.dat"	# File name of sequence Y
## Result files ##
SFILE="result.txt"	# Result file name for similar subsequence pairs
PFILE="path.txt"	# Result file name for optimal warping paths
## Parameters ##
LMIN="1500"		# Subsequence length threshold: lmin
EPS="0.05"		# Distance threshold: epsilon
BAND="5000"		# Width of Sakoe-Chiba band: w
XLEN="10000"		# Sequence length of X
YLEN="10000"		# Sequence length of Y


## Detect similar subsequence pairs ##
$SDIR/crossmatch $DDIR/$SEQX $DDIR/$SEQY $LMIN $EPS $BAND > $RESULT/$SFILE

set -- `wc $RESULT/$SFILE`
lines=$1

if [ $lines -ge 1 ]; then

  ## Compute warping paths ##
  perl $SDIR/createfile.pl $RESULT/$SFILE $RESULT/exefile.sh $RESULT/checkfile $DDIR/$SEQX $DDIR/$SEQY $LMIN $EPS $BAND $RESULT/$PFILE
  sh $RESULT/exefile.sh
  $SDIR/path $RESULT/checkfile $EPS $YLEN

  ## Plot warping paths ##
  set -- `wc $RESULT/$PFILE`
  plines=$1
  if [ $plines -ge 1 ]; then
    perl $SDIR/gnuplot.pl $RESULT $XLEN $YLEN
    gnuplot $RESULT/load_dat
#    convert $RESULT/fig.eps $RESULT/path.jpg
    rm $RESULT/load_dat
  fi

  ## Remove temporal files ##
  rm $RESULT/Subseq*.txt
  rm $SDIR/temp.txt
  rm $RESULT/exefile.sh
  rm $RESULT/checkfile

fi
