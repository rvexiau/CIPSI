# Makefile
NAME=ciro
FFLAGS= -c -O2 -shared-intel -mcmodel=large 
LFLAGS= -O2 -I$(MKLROOT)/include -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -mcmodel=large
# FFLAGS= -c -i_dynamic -mcmodel=large
# LFLAGS= -i_dynamic -mcmodel=large
# FFLAGS= -c -O2 -m32 -fbounds-check
# LFLAGS= -m32 -fbounds-check
COMP=ifort
OBJS = \
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
%.o: %.f
	$(COMP) $(FFLAGS) $*.f
%.o: %.f90
	$(COMP) $(FFLAGS) $*.f90	

ciro:	$(OBJS)
	$(COMP) $(LFLAGS) $(OBJS) -o $(NAME)
	rm -f ../../bin/$(NAME)
	cp $(NAME) ../../bin/
	
clean:
	rm $(OBJS) $(NAME)

