# Makefile, optional argument : define CLUSTER to use -axAVX,SSE4.2 flag

ifdef CLUSTER
    COMP= ifort
    OPT_LVL= -I$(MKLROOT)/include -O3 -heap-arrays -openmp -mkl=parallel -mcmodel=medium -shared-intel -axAVX,SSE4.2
    LIB= -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core  
else
    CHECK := $(shell command -V ifort 2> /dev/null)
    ifndef CHECK 				# ifort not available use gfortran
        COMP= gfortran
        OPT_LVL= -std=legacy -O3 -openmp -mcmodel=medium -ffree-line-length-none
        LIB= -lblas -llapack        
    else					# use ifort
        COMP= ifort  
        OPT_LVL= -I$(MKLROOT)/include -O3 -heap-arrays -qopenmp -mkl=parallel -mcmodel=medium -shared-intel -xhost
        LIB= -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core 
    endif
endif

Foldauto= ./
OBJSauto = ${Foldauto}timestamp.o ${Foldauto}lhs.o ${Foldauto}io_unit.o ${Foldauto}m_mrgrnk.o ${Foldauto}makegrid.o ${Foldauto}autocip13.o ${Foldauto}masses.o\
${Foldauto}tri_mono.o ${Foldauto}makedata.o ${Foldauto}makedata_pshf.o ${Foldauto}makedata_cip.o ${Foldauto}makedata_ciro.o ${Foldauto}calcul.o\
${Foldauto}ffield.o ${Foldauto}utils.o ${Foldauto}sym_check.o ${Foldauto}spin.o

all: Autocip13.exe

	
Autocip13.exe:	$(OBJSauto)
	$(COMP) $(OPT_LVL) $(LIB) -o ../bin/$@ $^ 

.SUFFIXES: .o .f90 .f

.f90.o:
	$(COMP) $(OPT_LVL) -c $^ -o $@	
.f.o:
	$(COMP) $(OPT_LVL) -c $^ -o $@		
	
	
clean:
	rm -f ../bin/Autocip13.exe Autocip13.exe *.mod

