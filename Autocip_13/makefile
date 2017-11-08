# Makefile
NAME=Autocip13

FFLAGS= -I$(MKLROOT)/include -O3 -c -openmp -mkl=parallel
LFLAGS= -I$(MKLROOT)/include -O3 -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lpthread -openmp -mkl=parallel

#FFLAGS= -I$(MKLROOT)/include -g -c -openmp -mkl=parallel
#LFLAGS= -I$(MKLROOT)/include -g -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lpthread -openmp -mkl=parallel

COMP= ifort
# COMP= gfortran
OBJS = \
timestamp.o\
lhs.o\
io_unit.o\
m_mrgrnk.o\
makegrid.o\
autocip13.o\
masses.o\
tri_mono.o\
makedata.o\
makedata_pshf.o\
makedata_cip.o\
makedata_ciro.o\
calcul.o\
ffield.o\
utils.o\
sym_check.o\
spin.o

.SUFFIXES: .o .f90 .f

.f90.o:
	$(COMP) $(FFLAGS) -c $^ -o $@	
	
Autocip13:	$(OBJS)
	$(COMP) $(LFLAGS) $(OBJS) -o $(NAME)
	rm -f ../bin/$(NAME)
	cp $(NAME) ../bin
	make clean
	
clean:
	rm $(OBJS) $(NAME) *.mod

