# Makefile
NAME=cip
CC=icc
# CFLAGS= -c -shared-intel -mcmodel=large
# FFLAGS= -c -shared-intel -mcmodel=large
# LFLAGS= -shared-intel -mcmodel=large
# CFLAGS= -c -mcmodel=large
FFLAGS= -O2 -I$(MKLROOT)/include -c -mcmodel=large
LFLAGS= -O2 -I$(MKLROOT)/include -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -mcmodel=large
# CFLAGS= -c -i_dynamic -mcmodel=large
# FFLAGS= -c -i_dynamic -mcmodel=large
# LFLAGS= -i_dynamic -mcmodel=large
# CFLAGS= -c -O2 -m32
# FFLAGS= -c -O2 -m32
# LFLAGS= -m32
# COMP=gfortran
COMP=ifort
OBJS = \
 diagonaliser.o\
 ai.o\
 brdnee.o\
 cas.o\
 cip.o\
 cnvmul.o\
 cotra.o\
 deter.o\
 diexcit.o\
 faifoc.o\
 given.o\
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
%.o: %.f
	$(COMP) $(FFLAGS) $*.f
%.o: %.f90
	$(COMP) $(FFLAGS) $*.f90	

cip:	$(OBJS)
	$(COMP) $(LFLAGS) $(OBJS) -o $(NAME) 
	rm -f ../../bin/$(NAME)
	cp $(NAME) ../../bin/$(NAME)

clean:
	rm $(OBJS) $(NAME)
  
