# Makefile
NAME=rcut.exe

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

OBJSrcut = \
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

all: rcut.exe

rcut.exe:	$(OBJSrcut)
	$(COMP) $(OPT_LVL) $(LIB) -o ../../bin/$@ $^ 
	
.SUFFIXES: .o .f90 .f

.f90.o:
	$(COMP) $(OPT_LVL) -c $^ -o $@	
.f.o:
	$(COMP) $(OPT_LVL) -c $^ -o $@	
	
clean:
	rm $(OBJS) ../../bin/$(NAME)


