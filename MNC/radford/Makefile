CFLAGS=-O -ansi -pedantic -g 
CC = gcc

all:		gen denc ddec unpck

gen:		gen.o rand.o
	$(CC) gen.o rand.o -lm -o gen

gen.o:		gen.c rand.h

unpck:		unpck.o
	$(CC) unpck.o -lm -o unpck

unpck.o:	unpck.c 

denc:		denc.o
	$(CC) denc.o -lm -o denc

denc.o:		denc.c 

ddec:		ddec.o
	$(CC) ddec.o -lm -o ddec

ddec.o:		ddec.c 

rand.o:		rand.c rand.h

mod2mat-test:	mod2mat-test.o mod2mat.o
	$(CC) mod2mat-test.o mod2mat.o -lm -o mod2mat-test

mod2mat-test.o:	mod2mat-test.c mod2mat.h

invert:	invert.o mod2mat.o
	$(CC) invert.o mod2mat.o -lm -o invert

invert.o:	invert.c mod2mat.h

mod2mat-time:	mod2mat-time.o mod2mat.o rand.o
	$(CC) -p mod2mat-time.o mod2mat.o rand.o -lm -o mod2mat-time

mod2mat-time.o:	mod2mat-time.c mod2mat.h

mod2mat.o:	mod2mat.c mod2mat.h
