# Makefile
NAME=rcut.X
FFLAGS= -c
LFLAGS=
# FFLAGS= -c -fast
# LFLAGS= -fast
# FFLAGS= -c -O2 -parallel -ipo -fpp -openmp
# LFLAGS= -O2 -parallel -ipo -fpp -openmp
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
	cp $(NAME) ../bin/
	
clean:
	rm $(OBJS) $(NAME)
  

