# MAKEFILE FOR ansi library.

CFLAGS = $(INCDIRS) \
            	 -ansi -pedantic -g\
		 -Wall -Wconversion -Wstrict-prototypes\
                 -Wformat -Wmissing-prototypes -Wshadow\
                 -Wpointer-arith -Wcast-qual -Wwrite-strings\
		 -D__USE_FIXED_PROTOTYPES__
#  -O2 -funroll-loops
#
# to profile    gprof opt2
# to not profile, remove -pg here and below. 
#
# to optimize add -O2 -funroll-loops

CC = gcc

%.tar:	%.c
	cd .. ; tar cvf /home/mackay/pub/c/$*.tar  -I ansi/$*.tar.com
	compress -f /home/mackay/pub/c/$*.tar
clean:
	rm -f *.o

first:		first.o  	cg.o	nrutil.o
	$(CC) first.o cg.o nrutil.o -lm -o first

tmp:		tmp.o  	cg.o	nrutil.o
	$(CC) tmp.o cg.o nrutil.o -lm -o tmp

test_new_cg:	test_new_cg.o test_function.o 	new_cg.o nrutil.o r.o
	$(CC) test_new_cg.o test_function.o new_cg.o nrutil.o r.o -lm -o test_new_cg

t_cg_solve:	t_cg_solve.o test_function.o 	new_cg.o nrutil.o r.o mynr.h r.h
	$(CC) t_cg_solve.o test_function.o new_cg.o nrutil.o r.o -lm -o t_cg_solve

test_mac:	test_mac.o test_function.o 	cg.o nrutil.o r.o mynr.h r.h test.h
	$(CC) test_mac.o test_function.o cg.o nrutil.o r.o -lm -o test_mac

test_macII:	test_macII.o test_function.o macopt.o nrutil.o r.o 
	$(CC) test_macII.o test_function.o macopt.o nrutil.o r.o -lm -o test_macII

rand2.o:	rand2.c	rand2.h
rand2b.o:	rand2b.c	rand2.h
smartall.o:	smartall.c smartall.h
first.o:	first.c
tmp.o:		tmp.c
cg.o:		cg.c  	mynr.h
cmatrix.o:	cmatrix.c cmatrix.h
nrutil.o:	nrutil.c	nrutil.h
new_cg.o:	new_cg.c mynr.h r.h
r.o:		r.c r.h
t_cg_solve.o:	t_cg_solve.c r.h
test_macII.o:	mynr.h 	r.h 	macopt.h 	test.h 	nrutil.h 
macopt.o:		macopt.c	macopt.h
strangeopt.o:		strangeopt.c	strangeopt.h
test_function.o:	test_function.c test.h 


