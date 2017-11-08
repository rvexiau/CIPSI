# Makefile
NAME=cval.X
FFLAGS= -c -fast
LFLAGS= -fast
# FFLAGS= -c -O 
# LFLAGS=  
COMP=ifort
OBJS = \
 cval.o\
 daclos.o\
 nomf.o\
 nomfil.o\
 openf.o\
 polone.o\
 poltwo.o\
 polzer.o\
 reada.o\
 set.o\
 stkener.o\
 wrtda.o
.f.o:
	$(COMP) $(FFLAGS) $*.f

cval:	$(OBJS)
	$(COMP) $(LFLAGS) $(OBJS) -o $(NAME)
	cp $(NAME) ../bin/

clean:
	rm $(OBJS) $(NAME)
  
