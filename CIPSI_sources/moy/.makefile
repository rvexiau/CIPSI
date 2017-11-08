# Makefile
NAME=moy.X
FFLAGS= -c -i_dynamic -mcmodel=large
LFLAGS= -i_dynamic -mcmodel=large
# LFLAGS= 
# FFLAGS= -c
# LFLAGS= -fast
# FFLAGS= -c -g -O2 -m32
# LFLAGS= -g -O2 -m32
COMP=ifort
OBJS = \
moy.o\
ai.o\
combin.o\
deter.o\
hdig.o\
hmp.o\
hntd.o\
ijkf.o\
lecnew.o\
morue.o\
nomf.o\
nomfil.o\
openf.o\
pie.o\
reijkl.o\
spipr.o\
stkener.o\
stk62.o\
transf.o
.f.o:
	$(COMP) $(FFLAGS) $*.f

moy:	$(OBJS)
	$(COMP) $(LFLAGS) $(OBJS) -o $(NAME)
	cp $(NAME) ../bin/

clean:
	rm $(OBJS) $(NAME)
  
