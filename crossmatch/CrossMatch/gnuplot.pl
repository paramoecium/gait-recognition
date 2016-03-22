#!/usr/bin/perl

if (@ARGV == 3) {

  open (FPOUT, ">load_dat") || die "Cannot open write file.\n";

  printf (FPOUT "set term post eps 22\n");
  printf (FPOUT "unset key\n");
  printf (FPOUT "set xrange[0:%d]\n", $ARGV[1]);
  printf (FPOUT "set yrange[0:%d]\n", $ARGV[2]);
  printf (FPOUT "set out '%s/fig.eps'\n", $ARGV[0]);
  printf (FPOUT "plot \"%s/path.txt\"\n", $ARGV[0]);

}  else {

  printf "Usage: perl gnuplot.pl  'directory of result'  'sequence length of X'  'sequence length of Y'\n";

}
