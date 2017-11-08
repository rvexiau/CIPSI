# Makefile
NAME=rpol
FFLAGS= -c
LFLAGS=
# FFLAGS= -c -fast
# LFLAGS= -fast
# FFLAGS= -c -fast
# LFLAGS=  -fast
# FFLAGS= -c -m32 -O3  
# LFLAGS=  -m32
COMP=ifort
CCOMP=cc
# COMP=ifort
# CCOMP=icc
CFLAGS= -c -fast
# CFLAGS= -c -fast
OBJS = \
nomf.o  nomfil.o timdat.o rpol.o  efield.o  openf.o  set.o xyztrans.o\
rot3d.o polytospher.o wigrot.o xyz0.o xyzint.o rint.o thetaint.o bessi.o 
.f.o:
	$(COMP) $(FFLAGS) $*.f

.c.o:
	$(CCOMP) $(CFLAGS) $*.c

rcut:	$(OBJS)
	$(COMP) $(LFLAGS) $(OBJS) -o $(NAME)
	cp $(NAME) ../bin

clean:
	rm $(OBJS) $(NAME)
  
