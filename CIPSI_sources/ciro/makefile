# Makefile
NAME=ciro.exe

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

OBJSciro = \
hermit.o\
aolim.o\
atpop.o\
ciro.o\
comb.o\
comps.o\
comtab.o\
denint.o\
diagonaliser.o\
dipint.o\
dipole.o\
dvint.o\
escriu.o\
given.o\
grossc.o\
indice.o\
matout.o\
mulken.o\
nomf.o\
nomfil.o\
openf.o\
ovlpop.o\
prop.o\
ran.o\
rdener.o\
reada.o\
reijkl.o\
rntd.o\
rozero.o\
scrivi.o\
selec.o\
stkener.o\
tracp.o\
vint.o\
xyzdip.o

all: ciro.exe

ciro.exe:	$(OBJSciro)
	$(COMP) $(OPT_LVL) $(LIB) -o ../../bin/$@ $^ 
	
.SUFFIXES: .o .f90 .f

.f90.o:
	$(COMP) $(OPT_LVL) -c $^ -o $@	
.f.o:
	$(COMP) $(OPT_LVL) -c $^ -o $@		
	
clean:
	rm $(OBJS) ../../bin/$(NAME)

