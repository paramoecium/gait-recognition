printf 'Subseq1.txt  X[1771:2889], Y[3658:4856]\n'
./crossmatch-log ./match_data/signal.dat ./match_data/signal.dat 300 0.175000 5000.000000 1771 3658 2889 4856 > Subseq1.txt
printf 'Subseq2.txt  X[1562:3816], Y[3155:5476]\n'
./crossmatch-log ./match_data/signal.dat ./match_data/signal.dat 300 0.175000 5000.000000 1562 3155 3816 5476 > Subseq2.txt
printf 'Subseq3.txt  X[1662:4182], Y[3004:5526]\n'
./crossmatch-log ./match_data/signal.dat ./match_data/signal.dat 300 0.175000 5000.000000 1662 3004 4182 5526 > Subseq3.txt
printf 'Subseq4.txt  X[2943:4451], Y[4076:5525]\n'
./crossmatch-log ./match_data/signal.dat ./match_data/signal.dat 300 0.175000 5000.000000 2943 4076 4451 5525 > Subseq4.txt
printf 'Subseq5.txt  X[1:5527], Y[1:5527]\n'
./crossmatch-log ./match_data/signal.dat ./match_data/signal.dat 300 0.175000 5000.000000 1 1 5527 5527 > Subseq5.txt
printf 'Subseq6.txt  X[3658:4856], Y[1771:2889]\n'
./crossmatch-log ./match_data/signal.dat ./match_data/signal.dat 300 0.175000 5000.000000 3658 1771 4856 2889 > Subseq6.txt
printf 'Subseq7.txt  X[4076:5525], Y[2943:4451]\n'
./crossmatch-log ./match_data/signal.dat ./match_data/signal.dat 300 0.175000 5000.000000 4076 2943 5525 4451 > Subseq7.txt
printf 'Subseq8.txt  X[3004:5526], Y[1662:4182]\n'
./crossmatch-log ./match_data/signal.dat ./match_data/signal.dat 300 0.175000 5000.000000 3004 1662 5526 4182 > Subseq8.txt
printf 'Subseq9.txt  X[3155:5476], Y[1562:3816]\n'
./crossmatch-log ./match_data/signal.dat ./match_data/signal.dat 300 0.175000 5000.000000 3155 1562 5476 3816 > Subseq9.txt
