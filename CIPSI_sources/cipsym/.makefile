# Makefile
NAME=cip.X
CC=icc
CFLAGS= -c -i_dynamic -mcmodel=large
FFLAGS= -c -i_dynamic -mcmodel=large
LFLAGS= -i_dynamic -mcmodel=large
# CFLAGS= -c -O2 -m32
# FFLAGS= -c -O2 -m32
# LFLAGS= -m32
# COMP=gfortran
COMP=ifort
OBJS = \
 ai.o\
 brdnee.o\
 cas.o\
 cip.o\
 cnvmul.o\
 cotra.o\
 deter.o\
 diexcit.o\
 faifoc.o\
 hdig.o\
 hmp.o\
 hntd.o\
 ic001.o\
 igene.o\
 ijkf.o\
 initab.o\
 initsm.o\
 isto.o\
 jacscf.o\
 lec60.o\
 lecnew.o\
 mult.o\
 nomf.o\
 nomfil.o\
 openf.o\
 pert1.o\
 pertu.o\
 pie.o\
 ran.o\
 recsym.o\
 reijkl.o\
 rdener.o\
 shell.o\
 sntd.o\
 spasym.o\
 stkener.o\
 stktyp.o\
 transf.o
.f.o:
	$(COMP) $(FFLAGS) $*.f

cip:	$(OBJS)
	$(COMP) $(LFLAGS) $(OBJS) -o $(NAME)
	cp $(NAME) ../bin

clean:
	rm $(OBJS) $(NAME)
  
