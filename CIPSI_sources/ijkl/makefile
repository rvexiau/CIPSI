# Makefile
NAME=ijkl.exe

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

OBJSijkl = \
 iword_func.o\
 BLOCKDATA1.o\
 ex.o\
 ijkl.o\
 itijkl.o\
 nomf.o\
 nomfil.o\
 openf.o\
 otodon.o\
 read8.o\
 reada.o\
 readb.o\
 rondof.o\
 rotat1.o\
 rotat2.o\
 tijkld.o\
 tijkls.o\
 tmono.o\
 tpqkld.o\
 tpqkls.o\
 transd.o\
 transs.o\
 vector.o\
 wrtda.o\
 wrtdb.o
 
all: ijkl.exe 

ijkl.exe:	$(OBJSijkl)
	$(COMP) $(OPT_LVL) $(LIB) -o ../../bin/$@ $^ 
	
.SUFFIXES: .o .f90 .f

.f90.o:
	$(COMP) $(OPT_LVL) -c $^ -o $@	
.f.o:
	$(COMP) $(OPT_LVL) -c $^ -o $@		
	

clean:
	rm $(OBJS) ../../bin/$(NAME)

