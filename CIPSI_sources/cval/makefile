# Makefile
NAME=cval
FFLAGS= -c -O2 -shared-intel -mcmodel=large
LFLAGS= -O2 -mcmodel=large
FFLAGS= -c -g -shared-intel -mcmodel=large
LFLAGS= -g -mcmodel=large
# FFLAGS= -c -fast
# LFLAGS= -fast
# FFLAGS= -c -O0 -check uninit -CB
# LFLAGS= -O0 -check all -check uninit -CB
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
	rm -f ../../bin/$(NAME)
	cp $(NAME) ../../bin/

clean:
	rm $(OBJS) $(NAME)
  
