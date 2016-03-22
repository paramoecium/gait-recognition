CC = gcc
CFLAGS = -O2 -Wall
#CFLAGS = -g -O -Wall
#CFLAGS = -g -Wall
CFLAGS2 = -O2 -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64
#CFLAGS2 = -g -O -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64
#CFLAGS2 = -g -Wall -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64

all : sampling sampling-path sampling-log path

sampling : crosscompo-sampling.o sampling.o
	$(CC) $(CFLAGS) -o sampling crosscompo-sampling.o sampling.o -lm
sampling-path : crosscompo-sampling.o sampling-path.o
	$(CC) $(CFLAGS) -o sampling-path crosscompo-sampling.o sampling-path.o -lm
sampling-log : crosscompo-sampling.o sampling-log.o
	$(CC) $(CFLAGS) -o sampling-log crosscompo-sampling.o sampling-log.o -lm
path : path.o
	$(CC) $(CFLAGS2) -o path path.o -lm

crosscompo-sampling.o : crosscompo-sampling.c 
	$(CC) $(CFLAGS) -c crosscompo-sampling.c
sampling.o : sampling.c sampling.h
	$(CC) $(CFLAGS) -c sampling.c
sampling-path.o : sampling-path.c sampling.h
	$(CC) $(CFLAGS) -c sampling-path.c
sampling-log.o : sampling-log.c sampling.h
	$(CC) $(CFLAGS) -c sampling-log.c
path.o : path.c
	$(CC) $(CFLAGS2) -c path.c

clean :
	rm -rf core *.o *~ sampling sampling-path sampling-log path
