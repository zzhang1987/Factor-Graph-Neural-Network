# MAKEFILE FOR `code' 

CFLAGS = $(INCDIRS) \
	      	 -ansi -pedantic -g \
		 -Wall -Wconversion -Wstrict-prototypes\
                 -Wformat -Wmissing-prototypes -Wshadow\
                 -Wpointer-arith -Wcast-qual -Wwrite-strings
#  -O2 -funroll-loops
#
# to profile    gprof opt2
# to not profile, remove -pg here and below. 
#
# to optimize add -O2 -funroll-loops
CC = gcc
# 

# LIBS = -L/home/mackay/lib -ldbmalloc -lm 
LIBS = -lm

clean:
	rm -f *.o

all:	r2b s2t t2r t2y y2b r2z zb2x

tar:
	tar cvf ../../pub/code/MNC.tar ./*.c ./*.h ./*_var ./commander.p ./Makefile ./ansi/*.c ./ansi/*.h ./ansi/Makefile ./bnd/bnd.c ./bnd/bnd.h ./bnd/Makefile ./bnd/bnd_var* ./bnd/commander.p ./bnd/ansi ./SUMMARY _s/16 _s/4 _s/12 _r/12 _b/12 _y/12 _t/12 _z/12 _x/4 radford/mod2mat.h radford/mod2mat.c radford/Makefile codes/12.4.3.111 util/randfile
	gzip -f  ../../pub/code/MNC.tar

zb2x: zb2x.o ansi/r.o ansi/nrutil.o ansi/cmatrix.o bnd/bnd.o ansi/rand2b.o
	$(CC) zb2x.o ansi/r.o ansi/nrutil.o ansi/cmatrix.o bnd/bnd.o ansi/rand2b.o $(LIBS) -o zb2x

ansi/smartall.o:	ansi/smartall.c ansi/smartall.h
	cd ansi; make smartall.o
zb2x.o:	zb2x.c zb2x.h ansi/r.h ansi/rand2.h ansi/mynr.h ansi/cmatrix.h bnd/bnd_var_str.h zb2x_var_str.h
s2t.o:	s2t.c ansi/r.h ansi/rand2.h ansi/mynr.h ansi/cmatrix.h s2t_var_str.h radford/mod2mat.h
s2t: s2t.o ansi/r.o ansi/nrutil.o ansi/cmatrix.o ansi/rand2b.o radford/mod2mat.o
	$(CC) s2t.o ansi/r.o ansi/nrutil.o ansi/cmatrix.o ansi/rand2b.o radford/mod2mat.o $(LIBS) -o s2t
s2t_var_str.h:	s2t_var commander.p
	commander.p s2t_var
makecode.o:	makecode.c ansi/r.h ansi/mynr.h ansi/cmatrix.h mkcode_var_str.h radford/mod2mat.h
makecode: makecode.o ansi/r.o ansi/nrutil.o ansi/cmatrix.o ansi/rand2b.o radford/mod2mat.o
	$(CC) makecode.o ansi/r.o ansi/nrutil.o ansi/cmatrix.o ansi/rand2b.o radford/mod2mat.o $(LIBS) -o makecode
mkcode_var_str.h:	mkcode_var commander.p
	commander.p mkcode_var
t2r.o:	t2r.c ansi/r.h ansi/rand2.h ansi/mynr.h t2r_var_str.h
t2r: t2r.o ansi/r.o ansi/nrutil.o ansi/cmatrix.o ansi/rand2b.o
	$(CC) t2r.o ansi/r.o ansi/nrutil.o ansi/rand2b.o $(LIBS) -o t2r
t2y.o:	t2y.c ansi/r.h ansi/rand2.h ansi/mynr.h t2y_var_str.h
t2y: t2y.o ansi/r.o ansi/nrutil.o ansi/rand2b.o
	$(CC) t2y.o ansi/r.o ansi/nrutil.o ansi/rand2b.o $(LIBS) -o t2y
t2y_var_str.h:	t2y_var commander.p
	commander.p t2y_var
r2z.o:	r2z.c ansi/r.h ansi/rand2.h ansi/mynr.h r2z_var_str.h
r2z: r2z.o ansi/r.o ansi/nrutil.o ansi/cmatrix.o ansi/rand2b.o
	$(CC) r2z.o ansi/r.o ansi/nrutil.o ansi/cmatrix.o  ansi/rand2b.o $(LIBS) -o r2z
r2z_var_str.h:	r2z_var commander.p
	commander.p r2z_var
r2b.o:	r2b.c ansi/r.h ansi/rand2.h ansi/mynr.h r2b_var_str.h
r2b: r2b.o ansi/r.o ansi/nrutil.o ansi/cmatrix.o ansi/rand2b.o
	$(CC) r2b.o ansi/r.o ansi/nrutil.o ansi/cmatrix.o  ansi/rand2b.o $(LIBS) -o r2b
r2b_var_str.h:	r2b_var commander.p
	commander.p r2b_var
y2b.o:	y2b.c ansi/r.h ansi/mynr.h y2b_var_str.h
y2b: y2b.o ansi/r.o ansi/nrutil.o ansi/cmatrix.o
	$(CC) y2b.o ansi/r.o ansi/nrutil.o $(LIBS) -o y2b
y2b_var_str.h:	y2b_var commander.p
	commander.p y2b_var

radford/mod2mat.o: radford/mod2mat.c radford/mod2mat.h
	cd radford ; make mod2mat.o
bnd/bnd.o: 	bnd/bnd.c bnd/bnd.h
	cd bnd; make bnd.o
bnd/bnd_var_str.h: 	bnd/bnd_var bnd/commander.p
	cd bnd; make bnd_var_str.h
t2r_var_str.h:	t2r_var commander.p
	commander.p t2r_var
zb2x_var_str.h:	zb2x_var commander.p
	commander.p zb2x_var

ansi/cg.o:		ansi/cg.c ansi/mynr.h
	cd ansi; make cg.o
ansi/nrutil.o:	ansi/nrutil.c	ansi/nrutil.h
	cd ansi; make nrutil.o
ansi/new_cg.o:	ansi/new_cg.c ansi/mynr.h ansi/r.h
	cd ansi; make new_cg.o
ansi/r.o:		ansi/r.c ansi/r.h
	cd ansi; make r.o
ansi/cmatrix.o:	ansi/cmatrix.c ansi/cmatrix.h
	cd ansi; make cmatrix.o
ansi/rand2b.o: ansi/rand2b.c ansi/rand2.h
	cd ansi; make rand2b.o

