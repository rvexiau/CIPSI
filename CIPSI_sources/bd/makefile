# Makefile
NAME=bd

# Parallèle :
# FFLAGS= -openmp -I$(MKLROOT)/include -c
# LFLAGS= -openmp -I$(MKLROOT)/include -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread

# Séquentiel
FFLAGS= -I$(MKLROOT)/include -c -O2
LFLAGS= -O2 -I$(MKLROOT)/include -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread


# FFLAGS= -c -fast
# LFLAGS= -fast
# FFLAGS= -c -O2 -parallel -ipo -fpp -openmp
# LFLAGS= -O2 -parallel -ipo -fpp -openmp -L/opt/intel/fce/10.1.015/lib -L/opt/intel/mkl/10.0.1.014/lib/em64t \
# -lmkl_lapack -lmkl -lguide -lpthread
# FFLAGS= -c
# LFLAGS= -L/opt/intel/fce/actual/lib -L/opt/intel/mkl/actual/lib/em64t \
# -lmkl_lapack -lmkl -lguide -lpthread
# LFLAGS= -L/opt/intel/fce/actual/lib -L/opt/intel/mkl/actual/lib/em64t \
# -lmkl_lapack -lmkl_em64t -lguide -lpthread
# FFLAGS= -c -O2 -parallel -ipo -fpp -openmp
# LFLAGS= -O2 -parallel -ipo -fpp -openmp -L/opt/intel/fce/actual/lib -L/opt/intel/mkl/actual/lib/em64t \
# -lmkl_lapack -lmkl_em64t -lguide -lpthread
# FFLAGS= -c -O2 -m32 
# LFLAGS= -O2 -m32 

COMP= ifort
OBJS = \
bd.o\
bigdia.o\
extr.o\
getmat.o\
impression_vecteurs.o\
jacscf.o\
low.o\
model.o\
mproj.o\
nomf.o\
nomfil.o\
openf.o\
oproj.o\
rdener.o\
rea2.o\
stkener.o\
storic.o\
trf.o\
wvec.o\
diagonaliser.o

%.o: %.f
	$(COMP) $(FFLAGS) $*.f
%.o: %.f90
	$(COMP) $(FFLAGS) $*.f90		
	
bd:	$(OBJS)
	$(COMP) $(LFLAGS) $(OBJS) -o $(NAME)
	rm -f ../../bin/$(NAME)
	cp $(NAME) ../../bin
	
clean:
	rm $(OBJS) $(NAME)

