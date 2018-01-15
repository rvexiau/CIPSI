# Makefile
NAME=cipsym.exe

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

OBJScipsym = \
 diagonaliser.o\
 ai.o\
 brdnee.o\
 cas.o\
 cip.o\
 cnvmul.o\
 cotra.o\
 deter.o\
 diexcit.o\
 faifoc.o\
 hdig.o\
 hmp.o\
 hntd.o\
 ic001.o\
 igene.o\
 ijkf.o\
 initab.o\
 initsm.o\
 isto.o\
 jacscf.o\
 lec60.o\
 lecnew.o\
 mult.o\
 nomf.o\
 nomfil.o\
 openf.o\
 pert1.o\
 pertu.o\
 pie.o\
 ran.o\
 recsym.o\
 reijkl.o\
 rdener.o\
 shell.o\
 sntd.o\
 spasym.o\
 stkener.o\
 stktyp.o\
 transf.o
 
all: cipsym.exe

cipsym.exe:	$(OBJScipsym)
	$(COMP) $(OPT_LVL) $(LIB) -o ../../bin/$@ $^ 

.SUFFIXES: .o .f90 .f

.f90.o:
	$(COMP) $(OPT_LVL) -c $^ -o $@	
.f.o:
	$(COMP) $(OPT_LVL) -c $^ -o $@		

clean:
	rm $(OBJS) ../../bin/$(NAME)

