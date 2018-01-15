# Makefile
NAME=rpol.exe

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

OBJSrpol = \
nomf.o  nomfil.o timdat.o rpol.o  efield.o  openf.o  set.o xyztrans.o\
rot3d.o polytospher.o wigrot.o xyz0.o xyzint.o rint.o thetaint.o bessi.o 

all: rpol.exe
	
rpol.exe:	$(OBJSrpol)
	$(COMP) $(OPT_LVL) $(LIB) -o ../../bin/$@ $^ 
	
.SUFFIXES: .o .f90 .f

.f90.o:
	$(COMP) $(OPT_LVL) -c $^ -o $@	
.f.o:
	$(COMP) $(OPT_LVL) -c $^ -o $@		

clean:
	rm $(OBJS) ../../bin/$(NAME)

