# Makefile
NAME=ijkl
FFLAGS= -c -O2 -mcmodel=large
LFLAGS= -O2 -L/opt/intel/lib/intel64/ -mcmodel=large
# FFLAGS= -c -fast
# LFLAGS= -fast
# FFLAGS= -c -m32 -O2  
# LFLAGS= -m32
COMP=ifort
OBJS = \
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
.f.o:
	$(COMP) $(FFLAGS) $*.f

ijkl:	$(OBJS)
	$(COMP) $(LFLAGS) $(OBJS) -o $(NAME)
	rm -f ../../bin/$(NAME)
	cp $(NAME) ../../bin/

clean:
	rm $(OBJS) $(NAME)
  
