# Makefile
NAME=moysym
# FFLAGS= -c -shared-intel -mcmodel=large
# LFLAGS= -shared-intel -mcmodel=large
FFLAGS= -c -O2 -mcmodel=large
LFLAGS= -O2 -L/opt/intel/lib/intel64/ -mcmodel=large


# FFLAGS= -c -i_dynamic -mcmodel=large
# LFLAGS= -i_dynamic -mcmodel=large
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
degen.o\
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
transf.o\
wrthij.o\
wrtempt.o\
wrthijtot.o
%.o: %.f
	$(COMP) $(FFLAGS) $*.f
%.o: %.f90
	$(COMP) $(FFLAGS) $*.f90

moysym:	$(OBJS)
	$(COMP) $(LFLAGS) $(OBJS) -o $(NAME)
	rm -f ../../bin/$(NAME)
	cp $(NAME) ../../bin/

clean:
	rm $(OBJS) $(NAME)
  
