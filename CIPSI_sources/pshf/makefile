# Makefile
NAME=pshf
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
iword_func.o\
addmat.o\
aolim.o\
atoms.o\
atpop.o\
backoa.o\
backtr.o\
blkdta000.o\
buildz.o\
cale.o\
chelec.o\
comgrp.o\
comone.o\
comps.o\
comtab.o\
comtwo.o\
daclos.o\
dafile.o\
daopen.o\
debut.o\
denhf.o\
denint.o\
dfbsm1.o\
dfunc.o\
diagonaliser.o\
dipint.o\
dipole.o\
dmop.o\
dmtx.o\
dnode.o\
dresyl.o\
droot.o\
dsmit.o\
ecrita.o\
elimat.o\
elimv.o\
error.o\
final.o\
find.o\
flopoi.o\
focktr.o\
forms.o\
fout.o\
ftrgel.o\
ftrom.o\
gelom.o\
genral.o\
geoda.o\
getgel.o\
grossc.o\
hdiag.o\
hfprop.o\
hstar.o\
hstar2.o\
ijprim.o\
insert.o\
instat.o\
intga.o\
intgab.o\
isoin.o\
isoout.o\
jandk.o\
ligen.o\
ligen2.o\
local.o\
matbc.o\
matout.o\
mole.o\
mulken.o\
nesbet.o\
nomf.o\
nomfil.o\
opener.o\
openf.o\
orthon.o\
ortit.o\
ovlpop.o\
pert.o\
pgdpld.o\
pgdplp.o\
pgdpls.o\
pgpplp.o\
pgppls.o\
pgspls.o\
potion.o\
print.o\
psepin.o\
psepot.o\
psgrin.o\
pshf.o\
ptgrp.o\
qmat.o\
qout.o\
rdener.o\
reada.o\
readb.o\
reduc.o\
rhr.o\
root4.o\
root5.o\
roper.o\
rot.o\
rt123.o\
s0000.o\
scf.o\
scfgel.o\
scfop.o\
secnd.o\
second.o\
shells.o\
sjpd.o\
spin.o\
spind.o\
standv.o\
stkener.o\
stvint.o\
symh.o\
timit.o\
tracep.o\
tracp.o\
trans.o\
trmv.o\
trvm.o\
vec.o\
vecout.o\
vout.o\
vprod.o\
wroot4.o\
wroot5.o\
wrt.o\
wrt123.o\
wrtda.o\
wrtdb.o\
xpsgr.o\
xpsnlc.o\
xpsslc.o\
xyzdip.o\
xyzint.o
%.o: %.f
	$(COMP) $(FFLAGS) $*.f
%.o: %.f90
	$(COMP) $(FFLAGS) $*.f90
	
pshf:	$(OBJS)
	$(COMP) $(LFLAGS) $(OBJS) -o $(NAME)
	rm -f ../../bin/$(NAME)
	cp $(NAME) ../../bin/
	
clean:
	rm $(OBJS) $(NAME)
  
