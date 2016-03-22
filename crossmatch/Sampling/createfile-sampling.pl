#!/usr/bin/perl

if (@ARGV == 11) {

  open (FPIN, "<$ARGV[0]") || die "Cannot open read file.\n";		## Result file of CrossMatch 
  open (FPOUT, ">$ARGV[1]") || die "Cannot open write file.\n";		## Execution file for log output 
  open (FPOUT2, ">$ARGV[2]") || die "Cannot open write file2.\n";	## File for path detection 

  $i = 0;

  while ($point = <FPIN>) {

    $i++;

    ## Get starting and end positions 
    @pos = $point =~ /\d+/g;
    $x_sp = $pos[0];
    $x_ep = $pos[1];
    $y_sp = $pos[2];
    $y_ep = $pos[3];

    ## Determine log file name
    $file = join "", "Subseq", $i, ".txt";
    $break = '\n';

    ## Create execution file for log output
    printf (FPOUT "printf '%s  X[%d:%d], Y[%d:%d]%s'\n", $file, $x_sp, $x_ep, $y_sp, $y_ep, $break);
    printf (FPOUT "./sampling-log %s %s %d %f %f %d %d %d %d %d %d > %s\n",
	     $ARGV[3], $ARGV[4], $ARGV[5], $ARGV[6], $ARGV[7], $ARGV[8], $ARGV[9], $x_sp, $y_sp, $x_ep, $y_ep, $file);

    ## Create file for path detection 
    printf (FPOUT2 "%s  %s  %d  %d  %d  %d\n", $file, $ARGV[10], $x_ep, $y_ep, $x_sp, $y_sp);

  }

} else {

  printf "Usage : perl createfile-sampling.pl 'result file of CrossMatch'  'write file name for log output'  'write file name for path detection'  'sequence X file name'  'sequence Y file name'  'lmin'  'epsilon'  'band width'  'sampling period of X'  'sampling period of Y'  'result file name for path output'\n";

}
