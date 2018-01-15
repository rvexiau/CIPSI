# Makefile
NAME=moy.exe

ifdef CLUSTER
    COMP= ifort
    OPT_LVL= -I$(MKLROOT)/include -O3 -heap-arrays -mkl=sequential -mcmodel=medium -shared-intel -axAVX,SSE4.2
    LIB= -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core  
else
    CHECK := $(shell command -V ifort 2> /dev/null)
    ifndef CHECK 				# ifort not available use gfortran
        COMP= gfortran
        OPT_LVL= -O3 -mcmodel=medium -ffree-line-length-none -std=legacy
        LIB= -lblas -llapack
    else					# use ifort
        COMP= ifort  
        OPT_LVL= -I$(MKLROOT)/include -O3 -heap-arrays -mkl=sequential -mcmodel=medium -shared-intel -xhost
        LIB= -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core 
    endif
endif

OBJSmoy = \
moy.o\
ai.o\
combin.o\
deter.o\
hdig.o\
hmp.o\
hntd.o\
ijkf.o\
lec60.o\
morue.o\
nomf.o\
nomfil.o\
openf.o\
pie.o\
reijkl.o\
spipr.o\
stkener.o\
stk62.o\
transf.o

all: moy.exe

moy.exe:	$(OBJSmoy)
	$(COMP) $(OPT_LVL) $(LIB) -o ../../bin/$@ $^ 

.SUFFIXES: .o .f90 .f

.f90.o:
	$(COMP) $(OPT_LVL) -c $^ -o $@	
.f.o:
	$(COMP) $(OPT_LVL) -c $^ -o $@		

clean:
	rm $(OBJS) ../../bin/$(NAME)

