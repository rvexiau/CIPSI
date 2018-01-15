# Makefile
NAME=cval.exe

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

OBJScval = \
 cval.o\
 daclos.o\
 nomf.o\
 nomfil.o\
 openf.o\
 polone.o\
 poltwo.o\
 polzer.o\
 reada.o\
 set.o\
 stkener.o\
 wrtda.o

all: cval.exe 

cval.exe:	$(OBJScval)
	$(COMP) $(OPT_LVL) $(LIB) -o ../../bin/$@ $^ 
	
.SUFFIXES: .o .f90 .f

.f90.o:
	$(COMP) $(OPT_LVL) -c $^ -o $@	
.f.o:
	$(COMP) $(OPT_LVL) -c $^ -o $@		

clean:
	rm $(OBJS) ../../bin/$(NAME)
  
