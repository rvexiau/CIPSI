# Makefile
NAME=pshf.exe

ifdef CLUSTER
    COMP= ifort
    OPT_LVL= -I$(MKLROOT)/include -O3 -heap-arrays -mkl=sequential -mcmodel=medium -shared-intel -axAVX,SSE4.2
    LIB= -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core  
else
    CHECK := $(shell command -V ifort 2> /dev/null)
    ifndef CHECK 				# ifort not available use gfortran
        COMP= gfortran
        OPT_LVL= -O3 -mcmodel=medium -ffree-line-length-none -std=legacy
        LIB= -lblas -llapack
    else					# use ifort
        COMP= ifort  
        OPT_LVL= -I$(MKLROOT)/include -O3 -heap-arrays -mkl=sequential -mcmodel=medium -shared-intel -xhost
        LIB= -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core 
    endif
endif


OBJSpshf = \
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

all: pshf.exe

pshf.exe:	$(OBJSpshf)
	$(COMP) $(OPT_LVL) $(LIB) -o ../../bin/$@ $^ 

.SUFFIXES: .o .f90 .f

.f90.o:
	$(COMP) $(OPT_LVL) -c $^ -o $@	
.f.o:
	$(COMP) $(OPT_LVL) -c $^ -o $@	
	
	
clean:
	rm $(OBJS) ../../bin/$(NAME)

