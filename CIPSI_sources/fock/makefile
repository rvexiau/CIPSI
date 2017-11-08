# Makefile
NAME=fock
FFLAGS= -c -O2 -mcmodel=large
LFLAGS= -O2 -mcmodel=large
# FFLAGS= -c -fast
# LFLAGS= -fast
# FFLAGS= -c  -O  
# LFLAGS= 
COMP=ifort
OBJS = \
 ai.o\
 fock.o\
 itijkl_new.o\
 itijkl_stan.o\
 nomf.o\
 nomfil.o\
 reijkl_new.o\
 reijkl_stan.o\
 openf.o\
 stkener.o\
 wijkl_new.o\
 wijkl_stan.o
.f.o:
	$(COMP) $(FFLAGS) $*.f

fock:	$(OBJS)
	$(COMP) $(LFLAGS) $(OBJS) -o $(NAME)
	rm -f ../../bin/$(NAME)
	cp $(NAME) ../../bin/

clean:
	rm $(OBJS) $(NAME)
  
