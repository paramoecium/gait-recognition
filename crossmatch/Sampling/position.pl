#!/usr/bin/perl

if (@ARGV == 4) {

  open (FPIN, "<$ARGV[0]") || die "Cannot open read file.\n";
  open (FPOUT, ">$ARGV[1]") || die "Cannot open write file.\n";

  while ($point = <FPIN>) {

    @list = split ' ', $point;
    printf (FPOUT "%5d  %5d\n", $list[0]*$ARGV[2], $list[1]*$ARGV[3]);

  }

} else {

  printf "Usage : perl position.pl 'read file name'  'write file name'  'sampling period of X'  'sampling period of Y'\n";

}
