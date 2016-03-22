CC = gcc
CFLAGS = -O2 -Wall
#CFLAGS = -g -O -Wall
#CFLAGS = -g -Wall
CFLAGS2 = -O2 -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64
#CFLAGS2 = -g -O -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64
#CFLAGS2 = -g -Wall -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64


##### CrossMatch #####

all : crossmatch crossmatch-log path

crossmatch : crosscompo.o crossmatch.o
	$(CC) $(CFLAGS) -o crossmatch crosscompo.o crossmatch.o -lm
crossmatch-log : crosscompo.o crossmatch-log.o
	$(CC) $(CFLAGS) -o crossmatch-log crosscompo.o crossmatch-log.o -lm
path : path.o
	$(CC) $(CFLAGS2) -o path path.o -lm

crosscompo.o : crosscompo.c 
	$(CC) $(CFLAGS) -c crosscompo.c
crossmatch.o : crossmatch.c crossmatch.h
	$(CC) $(CFLAGS) -c crossmatch.c
crossmatch-log.o : crossmatch-log.c crossmatch.h
	$(CC) $(CFLAGS) -c crossmatch-log.c
path.o : path.c
	$(CC) $(CFLAGS2) -c path.c

clean :
	rm -rf core *.o *~ crossmatch crossmatch-log path


##### Naive #####

naive : crosscompo.o naive.o
	$(CC) $(CFLAGS) -o naive crosscompo.o naive.o -lm
naive.o : naive.c crossmatch.h
	$(CC) $(CFLAGS) -c naive.c

