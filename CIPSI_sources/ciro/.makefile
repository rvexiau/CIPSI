# Makefile
NAME=ciro.X
FFLAGS= -c -i_dynamic -mcmodel=large
LFLAGS= -i_dynamic -mcmodel=large
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
.f.o:
	$(COMP) $(FFLAGS) $*.f

ciro:	$(OBJS)
	$(COMP) $(LFLAGS) $(OBJS) -o $(NAME)
	cp $(NAME) ../bin/
	
clean:
	rm $(OBJS) $(NAME)

