# Makefile
NAME=hdiab
#FFLAGS= -c -O0 -check uninit -CB
#LFLAGS= -O0 -check all -check uninit -CB
# FFLAGS= -c -O0 -g -ftrapuv
# LFLAGS= -O0 -g -ftrapuv
# FFLAGS= -c -fast
# LFLAGS= -fast
FFLAGS= -I$(MKLROOT)/include -c -mcmodel=large
LFLAGS= -O2 -I$(MKLROOT)/include -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -mcmodel=large
COMP=ifort
OBJS = \
 determ.o\
 hdiab.o\
 hrot.o\
 hdiag.f\
 hvp.o\
 maxovl.o\
 nomf.o\
 nomfil.o\
 openf.o\
 readic.o\
 wrt.o\
 wrtvec.o
%.o: %.f
	$(COMP) $(FFLAGS) $*.f
%.o: %.f90
	$(COMP) $(FFLAGS) $*.f90
	
hdiab:	$(OBJS)
	$(COMP) $(LFLAGS) $(OBJS) -o $(NAME)
	$(COMP) $(LFLAGS) som.f -o som	
	rm -f ../../bin/$(NAME)
	rm -f ../../bin/som
	cp $(NAME) ../../bin/
	cp som ../../bin/
	
clean:
	rm $(OBJS) $(NAME)

