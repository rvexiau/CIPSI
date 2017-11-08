# Makefile
NAME=rcut
# FFLAGS= -c -O0 -check uninit -CB
# LFLAGS= -O0 -check all -check uninit -CB
# FFLAGS= -c -fast
# LFLAGS= -fast
FFLAGS= -I$(MKLROOT)/include -c -O2 -mcmodel=large
LFLAGS= -O2 -I$(MKLROOT)/include -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -mcmodel=large
COMP=ifort
OBJS = \
 casgam.o\
 caspgam.o\
 caspgax.o\
 comps.o\
 comtab.o\
 dfbsm1.o\
 efield.o\
 gen.o\
 intab.o\
 intga.o\
 intgab.o\
 intgabx.o\
 intgax.o\
 intinf.o\
 intnumx.o\
 invers.o\
 nomf.o\
 nomfil.o\
 openf.o\
 preintn.o\
 psepin.o\
 rcut.o\
 set.o\
 varcom.o
.f.o:
	$(COMP) $(FFLAGS) $*.f

rcut:	$(OBJS)
	$(COMP) $(LFLAGS) $(OBJS) -o $(NAME)
	rm -f ../../bin/$(NAME)
	cp $(NAME) ../../bin/
	
clean:
	rm $(OBJS) $(NAME)
  

