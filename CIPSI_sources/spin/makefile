# Makefile
NAME=spin

# Parallèle :
# FFLAGS= -openmp -I$(MKLROOT)/include -c
# LFLAGS= -openmp -I$(MKLROOT)/include -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread

# Séquentiel
FFLAGS= -O2 -I$(MKLROOT)/include -c
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
spin.o\
mult.o\
nomf.o\
nomfil.o\
openf.o\
sntd.o\


%.o: %.f
	$(COMP) $(FFLAGS) $*.f
%.o: %.f90
	$(COMP) $(FFLAGS) $*.f90		
	
spin:	$(OBJS)
	$(COMP) $(LFLAGS) $(OBJS) -o $(NAME)
	rm -f ../../bin/$(NAME)
	cp $(NAME) ../../bin
	
clean:
	rm $(OBJS) $(NAME)

