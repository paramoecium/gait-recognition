## Directories ##
SDIR=.				# Directory of source files
DDIR=./data			# Directory of data sequences
RESULT=.			# Directory to output results
## Data files ##
SEQX="samp-sines1.dat"		# File name of sequence X
SEQY="samp-sines2.dat"		# File name of sequence Y
## Result files ##
SFILE="result-sampling.txt"	# Result file name for similar subsequence pairs
PFILE="path-sampling.txt"	# Result file name for optimal warping paths
## Parameters ##
LMIN="1500"			# Subsequence length threshold: lmin
EPS="0.05"			# Distance threshold: epsilon
BAND="5000"			# Width of Sakoe-Chiba band: w
PERIOD1="278"			# Sampling period of sequence X: Tx
PERIOD2="294"			# Sampling period of sequence Y: Ty
XLEN="10000"			# Original sequence length of X
YLEN="10000"			# Original sequence length of Y


## Detect similar subsequence pairs ##
sband=$(expr $BAND / $PERIOD2)
$SDIR/sampling $DDIR/$SEQX $DDIR/$SEQY $LMIN $EPS $sband $PERIOD1 $PERIOD2 > $RESULT/$SFILE
$SDIR/sampling-path $DDIR/$SEQX $DDIR/$SEQY $LMIN $EPS $sband $PERIOD1 $PERIOD2 > $RESULT/result-sampling-tmp.txt

set -- `wc $RESULT/result-sampling-tmp.txt`
lines=$1
length=$(expr $XLEN / $PERIOD2 + 50)

if [ $lines -ge 1 ]; then

  ## Compute warping paths ##
  perl $SDIR/createfile-sampling.pl $RESULT/result-sampling-tmp.txt $RESULT/exefile.sh $RESULT/checkfile $DDIR/$SEQX $DDIR/$SEQY $LMIN $EPS $sband $PERIOD1 $PERIOD2 $RESULT/path-tmp.txt
  sh $RESULT/exefile.sh
  $SDIR/path $RESULT/checkfile $EPS $length

  ## Plot warping paths ##
  set -- `wc $RESULT/path-tmp.txt`
  plines=$1
  if [ $plines -ge 1 ]; then
    perl $SDIR/position.pl $RESULT/path-tmp.txt $RESULT/$PFILE $PERIOD1 $PERIOD2
    perl $SDIR/gnuplot-sampling.pl $RESULT $XLEN $YLEN
    gnuplot $RESULT/load_dat
#    convert $RESULT/fig-sampling.eps $RESULT/path-sampling.jpg
  rm $RESULT/load_dat
  fi

  ## Remove temporal files ##
  rm $RESULT/Subseq*.txt
  rm $SDIR/temp.txt
  rm $RESULT/exefile.sh
  rm $RESULT/checkfile
  rm $RESULT/result-sampling-tmp.txt
  rm $RESULT/path-tmp.txt

fi
