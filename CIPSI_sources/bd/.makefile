# Makefile
NAME=bd
FFLAGS= -c -fast
LFLAGS= -fast
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
giveis.o\
givens.o\
impression_vecteurs.o\
jacscf.o\
low.o\
model.o\
mproj.o\
nomf.o\
nomfil.o\
openf.o\
oproj.o\
ran.o\
rdener.o\
rea2.o\
stkener.o\
storic.o\
trf.o\
wvec.o
.f.o:
	$(COMP) $(FFLAGS) $*.f

bd:	$(OBJS)
	$(COMP) $(LFLAGS) $(OBJS) -o $(NAME)
	cp $(NAME) ../bin
	
clean:
	rm $(OBJS) $(NAME)

