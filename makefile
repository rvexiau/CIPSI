# Makefile, optional argument : define CLUSTER to use -axAVX,SSE4.2 flag

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


# Séquentiel
# FFLAGS= -I$(MKLROOT)/include -c -O2 -mcmodel=large -shared-intel
# LFLAGS= -O2 -I$(MKLROOT)/include -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -mcmodel=large

Foldauto= Autocip_13/
OBJSauto = ${Foldauto}timestamp.o ${Foldauto}lhs.o ${Foldauto}io_unit.o ${Foldauto}m_mrgrnk.o ${Foldauto}makegrid.o ${Foldauto}autocip13.o ${Foldauto}masses.o\
${Foldauto}tri_mono.o ${Foldauto}makedata.o ${Foldauto}makedata_pshf.o ${Foldauto}makedata_cip.o ${Foldauto}makedata_ciro.o ${Foldauto}calcul.o\
${Foldauto}ffield.o ${Foldauto}utils.o ${Foldauto}sym_check.o ${Foldauto}spin.o

Foldbd= CIPSI_sources/bd/
OBJSbd = ${Foldbd}bd.o ${Foldbd}bigdia.o ${Foldbd}extr.o ${Foldbd}getmat.o ${Foldbd}impression_vecteurs.o ${Foldbd}jacscf.o ${Foldbd}low.o\
${Foldbd}model.o ${Foldbd}mproj.o ${Foldbd}nomf.o ${Foldbd}nomfil.o ${Foldbd}openf.o ${Foldbd}oproj.o ${Foldbd}rdener.o ${Foldbd}rea2.o ${Foldbd}stkener.o\
${Foldbd}storic.o ${Foldbd}trf.o ${Foldbd}wvec.o ${Foldbd}diagonaliser.o

Foldbdsym= CIPSI_sources/bdsym/
OBJSbdsym = ${Foldbdsym}bd.o ${Foldbdsym}bigdia.o ${Foldbdsym}extr.o ${Foldbdsym}getmat.o ${Foldbdsym}low.o ${Foldbdsym}model.o ${Foldbdsym}mproj.o ${Foldbdsym}nomf.o\
${Foldbdsym}nomfil.o ${Foldbdsym}openf.o ${Foldbdsym}oproj.o ${Foldbdsym}rdener.o ${Foldbdsym}rea2.o ${Foldbdsym}stkener.o ${Foldbdsym}storic.o ${Foldbdsym}trf.o\
${Foldbdsym}wvec.o ${Foldbdsym}diagonaliser.o

Foldcip= CIPSI_sources/cip/
OBJScip = ${Foldcip}diagonaliser.o ${Foldcip}ai.o ${Foldcip}brdnee.o ${Foldcip}cas.o ${Foldcip}cip.o ${Foldcip}cnvmul.o ${Foldcip}cotra.o ${Foldcip}deter.o\
 ${Foldcip}diexcit.o ${Foldcip}faifoc.o ${Foldcip}hdig.o ${Foldcip}hmp.o ${Foldcip}hntd.o ${Foldcip}ic001.o ${Foldcip}igene.o ${Foldcip}ijkf.o\
 ${Foldcip}initab.o ${Foldcip}initsm.o ${Foldcip}isto.o ${Foldcip}jacscf.o ${Foldcip}lec60.o ${Foldcip}lecnew.o ${Foldcip}mult.o ${Foldcip}nomf.o ${Foldcip}nomfil.o\
 ${Foldcip}openf.o ${Foldcip}pert1.o ${Foldcip}pertu.o ${Foldcip}pie.o ${Foldcip}ran.o ${Foldcip}recsym.o ${Foldcip}reijkl.o ${Foldcip}rdener.o ${Foldcip}shell.o\
 ${Foldcip}sntd.o ${Foldcip}spasym.o ${Foldcip}stkener.o ${Foldcip}stktyp.o ${Foldcip}transf.o

Foldcipsym= CIPSI_sources/cipsym/
OBJScipsym = ${Foldcipsym}diagonaliser.o ${Foldcipsym}ai.o ${Foldcipsym}brdnee.o ${Foldcipsym}cas.o ${Foldcipsym}cip.o ${Foldcipsym}cnvmul.o ${Foldcipsym}cotra.o ${Foldcipsym}deter.o\
 ${Foldcipsym}diexcit.o ${Foldcipsym}faifoc.o ${Foldcipsym}hdig.o ${Foldcipsym}hmp.o ${Foldcipsym}hntd.o ${Foldcipsym}ic001.o ${Foldcipsym}igene.o ${Foldcipsym}ijkf.o\
 ${Foldcipsym}initab.o ${Foldcipsym}initsm.o ${Foldcipsym}isto.o ${Foldcipsym}jacscf.o ${Foldcipsym}lec60.o ${Foldcipsym}lecnew.o ${Foldcipsym}mult.o ${Foldcipsym}nomf.o ${Foldcipsym}nomfil.o\
 ${Foldcipsym}openf.o ${Foldcipsym}pert1.o ${Foldcipsym}pertu.o ${Foldcipsym}pie.o ${Foldcipsym}ran.o ${Foldcipsym}recsym.o ${Foldcipsym}reijkl.o ${Foldcipsym}rdener.o ${Foldcipsym}shell.o\
 ${Foldcipsym}sntd.o ${Foldcipsym}spasym.o ${Foldcipsym}stkener.o ${Foldcipsym}stktyp.o ${Foldcipsym}transf.o

Foldciro= CIPSI_sources/ciro/
OBJSciro = ${Foldciro}hermit.o ${Foldciro}aolim.o ${Foldciro}atpop.o ${Foldciro}ciro.o ${Foldciro}comb.o ${Foldciro}comps.o ${Foldciro}comtab.o ${Foldciro}denint.o\
${Foldciro}diagonaliser.o ${Foldciro}dipint.o ${Foldciro}dipole.o ${Foldciro}dvint.o ${Foldciro}escriu.o ${Foldciro}given.o ${Foldciro}grossc.o ${Foldciro}indice.o\
${Foldciro}matout.o ${Foldciro}mulken.o ${Foldciro}nomf.o ${Foldciro}nomfil.o ${Foldciro}openf.o ${Foldciro}ovlpop.o ${Foldciro}prop.o ${Foldciro}ran.o ${Foldciro}rdener.o\
${Foldciro}reada.o ${Foldciro}reijkl.o ${Foldciro}rntd.o ${Foldciro}rozero.o ${Foldciro}scrivi.o ${Foldciro}selec.o ${Foldciro}stkener.o ${Foldciro}tracp.o\
${Foldciro}vint.o ${Foldciro}xyzdip.o

Foldcval= CIPSI_sources/cval/
OBJScval = ${Foldcval}cval.o ${Foldcval}daclos.o ${Foldcval}nomf.o ${Foldcval}nomfil.o ${Foldcval}openf.o ${Foldcval}polone.o ${Foldcval}poltwo.o ${Foldcval}polzer.o\
 ${Foldcval}reada.o ${Foldcval}set.o ${Foldcval}stkener.o ${Foldcval}wrtda.o

Foldfock= CIPSI_sources/fock/
OBJSfock = ${Foldfock}ai.o ${Foldfock}fock.o ${Foldfock}itijkl_new.o ${Foldfock}itijkl_stan.o ${Foldfock}nomf.o ${Foldfock}nomfil.o ${Foldfock}reijkl_new.o\
 ${Foldfock}reijkl_stan.o ${Foldfock}openf.o ${Foldfock}stkener.o ${Foldfock}wijkl_new.o ${Foldfock}wijkl_stan.o

Foldijkl= CIPSI_sources/ijkl/
OBJSijkl = ${Foldijkl}iword_func.o ${Foldijkl}BLOCKDATA1.o ${Foldijkl}ex.o ${Foldijkl}ijkl.o ${Foldijkl}itijkl.o ${Foldijkl}nomf.o ${Foldijkl}nomfil.o ${Foldijkl}openf.o\
 ${Foldijkl}otodon.o ${Foldijkl}read8.o ${Foldijkl}reada.o ${Foldijkl}readb.o ${Foldijkl}rondof.o ${Foldijkl}rotat1.o ${Foldijkl}rotat2.o ${Foldijkl}tijkld.o ${Foldijkl}tijkls.o\
 ${Foldijkl}tmono.o ${Foldijkl}tpqkld.o ${Foldijkl}tpqkls.o ${Foldijkl}transd.o ${Foldijkl}transs.o ${Foldijkl}vector.o ${Foldijkl}wrtda.o ${Foldijkl}wrtdb.o

Foldmoy= CIPSI_sources/moy/
OBJSmoy = ${Foldmoy}moy.o ${Foldmoy}ai.o ${Foldmoy}combin.o ${Foldmoy}deter.o ${Foldmoy}hdig.o ${Foldmoy}hmp.o ${Foldmoy}hntd.o ${Foldmoy}ijkf.o ${Foldmoy}lec60.o\
${Foldmoy}morue.o ${Foldmoy}nomf.o ${Foldmoy}nomfil.o ${Foldmoy}openf.o ${Foldmoy}pie.o ${Foldmoy}reijkl.o ${Foldmoy}spipr.o ${Foldmoy}stkener.o ${Foldmoy}stk62.o ${Foldmoy}transf.o

Foldmoysym= CIPSI_sources/moysym/
OBJSmoysym = ${Foldmoysym}moy.o ${Foldmoysym}ai.o ${Foldmoysym}combin.o ${Foldmoysym}degen.o ${Foldmoysym}deter.o ${Foldmoysym}hdig.o ${Foldmoysym}hmp.o ${Foldmoysym}hntd.o ${Foldmoysym}ijkf.o ${Foldmoysym}lecnew.o\
${Foldmoysym}morue.o ${Foldmoysym}nomf.o ${Foldmoysym}nomfil.o ${Foldmoysym}openf.o ${Foldmoysym}pie.o ${Foldmoysym}reijkl.o ${Foldmoysym}spipr.o ${Foldmoysym}stkener.o ${Foldmoysym}stk62.o ${Foldmoysym}transf.o\
${Foldmoysym}wrthij.o ${Foldmoysym}wrtempt.o ${Foldmoysym}wrthijtot.o

Foldpshf= CIPSI_sources/pshf/
OBJSpshf = ${Foldpshf}iword_func.o ${Foldpshf}addmat.o ${Foldpshf}aolim.o ${Foldpshf}atoms.o ${Foldpshf}atpop.o ${Foldpshf}backoa.o ${Foldpshf}backtr.o ${Foldpshf}blkdta000.o ${Foldpshf}buildz.o\
${Foldpshf}cale.o ${Foldpshf}chelec.o ${Foldpshf}comgrp.o ${Foldpshf}comone.o ${Foldpshf}comps.o ${Foldpshf}comtab.o ${Foldpshf}comtwo.o ${Foldpshf}daclos.o ${Foldpshf}dafile.o ${Foldpshf}daopen.o\
${Foldpshf}debut.o ${Foldpshf}denhf.o ${Foldpshf}denint.o ${Foldpshf}dfbsm1.o ${Foldpshf}dfunc.o ${Foldpshf}diagonaliser.o ${Foldpshf}dipint.o ${Foldpshf}dipole.o ${Foldpshf}dmop.o ${Foldpshf}dmtx.o\
${Foldpshf}dnode.o ${Foldpshf}dresyl.o ${Foldpshf}droot.o ${Foldpshf}dsmit.o ${Foldpshf}ecrita.o ${Foldpshf}elimat.o ${Foldpshf}elimv.o ${Foldpshf}error.o ${Foldpshf}final.o ${Foldpshf}find.o ${Foldpshf}flopoi.o\
${Foldpshf}focktr.o ${Foldpshf}forms.o ${Foldpshf}fout.o ${Foldpshf}ftrgel.o ${Foldpshf}ftrom.o ${Foldpshf}gelom.o ${Foldpshf}genral.o ${Foldpshf}geoda.o ${Foldpshf}getgel.o ${Foldpshf}grossc.o ${Foldpshf}hdiag.o\
${Foldpshf}hfprop.o ${Foldpshf}hstar.o ${Foldpshf}hstar2.o ${Foldpshf}ijprim.o ${Foldpshf}insert.o ${Foldpshf}instat.o ${Foldpshf}intga.o ${Foldpshf}intgab.o ${Foldpshf}isoin.o ${Foldpshf}isoout.o ${Foldpshf}jandk.o\
${Foldpshf}ligen.o ${Foldpshf}ligen2.o ${Foldpshf}local.o ${Foldpshf}matbc.o ${Foldpshf}matout.o ${Foldpshf}mole.o ${Foldpshf}mulken.o ${Foldpshf}nesbet.o ${Foldpshf}nomf.o ${Foldpshf}nomfil.o ${Foldpshf}opener.o\
${Foldpshf}openf.o ${Foldpshf}orthon.o ${Foldpshf}ortit.o ${Foldpshf}ovlpop.o ${Foldpshf}pert.o ${Foldpshf}pgdpld.o ${Foldpshf}pgdplp.o ${Foldpshf}pgdpls.o ${Foldpshf}pgpplp.o ${Foldpshf}pgppls.o ${Foldpshf}pgspls.o\
${Foldpshf}potion.o ${Foldpshf}print.o ${Foldpshf}psepin.o ${Foldpshf}psepot.o ${Foldpshf}psgrin.o ${Foldpshf}pshf.o ${Foldpshf}ptgrp.o ${Foldpshf}qmat.o ${Foldpshf}qout.o ${Foldpshf}rdener.o ${Foldpshf}reada.o\
${Foldpshf}readb.o ${Foldpshf}reduc.o ${Foldpshf}rhr.o ${Foldpshf}root4.o ${Foldpshf}root5.o ${Foldpshf}roper.o ${Foldpshf}rot.o ${Foldpshf}rt123.o ${Foldpshf}s0000.o ${Foldpshf}scf.o ${Foldpshf}scfgel.o ${Foldpshf}scfop.o\
${Foldpshf}secnd.o ${Foldpshf}second.o ${Foldpshf}shells.o ${Foldpshf}sjpd.o ${Foldpshf}spin.o ${Foldpshf}spind.o ${Foldpshf}standv.o ${Foldpshf}stkener.o ${Foldpshf}stvint.o ${Foldpshf}symh.o ${Foldpshf}timit.o\
${Foldpshf}tracep.o ${Foldpshf}tracp.o ${Foldpshf}trans.o ${Foldpshf}trmv.o ${Foldpshf}trvm.o ${Foldpshf}vec.o ${Foldpshf}vecout.o ${Foldpshf}vout.o ${Foldpshf}vprod.o ${Foldpshf}wroot4.o ${Foldpshf}wroot5.o ${Foldpshf}wrt.o\
${Foldpshf}wrt123.o ${Foldpshf}wrtda.o ${Foldpshf}wrtdb.o ${Foldpshf}xpsgr.o ${Foldpshf}xpsnlc.o ${Foldpshf}xpsslc.o ${Foldpshf}xyzdip.o ${Foldpshf}xyzint.o

Foldrcut= CIPSI_sources/rcut/
OBJSrcut = ${Foldrcut}casgam.o ${Foldrcut}caspgam.o ${Foldrcut}caspgax.o ${Foldrcut}comps.o ${Foldrcut}comtab.o ${Foldrcut}dfbsm1.o ${Foldrcut}efield.o ${Foldrcut}gen.o ${Foldrcut}intab.o ${Foldrcut}intga.o ${Foldrcut}intgab.o\
 ${Foldrcut}intgabx.o ${Foldrcut}intgax.o ${Foldrcut}intinf.o ${Foldrcut}intnumx.o ${Foldrcut}invers.o ${Foldrcut}nomf.o ${Foldrcut}nomfil.o ${Foldrcut}openf.o ${Foldrcut}preintn.o ${Foldrcut}psepin.o ${Foldrcut}rcut.o\
 ${Foldrcut}set.o ${Foldrcut}varcom.o

Foldrpol= CIPSI_sources/rpol/
OBJSrpol = ${Foldrpol}nomf.o ${Foldrpol}nomfil.o ${Foldrpol}timdat.o ${Foldrpol}rpol.o ${Foldrpol}efield.o ${Foldrpol}openf.o ${Foldrpol}set.o ${Foldrpol}xyztrans.o\
${Foldrpol}rot3d.o ${Foldrpol}polytospher.o ${Foldrpol}wigrot.o ${Foldrpol}xyz0.o ${Foldrpol}xyzint.o ${Foldrpol}rint.o ${Foldrpol}thetaint.o ${Foldrpol}bessi.o 

Foldspin= CIPSI_sources/spin/
OBJSspin = ${Foldspin}spin.o ${Foldspin}mult.o ${Foldspin}nomf.o ${Foldspin}nomfil.o ${Foldspin}openf.o ${Foldspin}sntd.o

Foldspinsym= CIPSI_sources/spinsym/
OBJSspinsym = ${Foldspinsym}spin.o ${Foldspinsym}mult.o ${Foldspinsym}nomf.o ${Foldspinsym}nomfil.o ${Foldspinsym}openf.o ${Foldspinsym}sntd.o

Foldhdiab= CIPSI_sources/hdiab/
OBJShdiab = ${Foldhdiab}determ.o ${Foldhdiab}hdiab.o ${Foldhdiab}hdiag.o ${Foldhdiab}hrot.o ${Foldhdiab}hvp.o ${Foldhdiab}maxovl.o ${Foldhdiab}nomf.o ${Foldhdiab}nomfil.o ${Foldhdiab}openf.o ${Foldhdiab}readic.o\
 ${Foldhdiab}wrt.o ${Foldhdiab}wrtvec.o

all: Autocip13 bd.exe bdsym.exe cip.exe cipsym.exe ciro.exe cval.exe fock.exe ijkl.exe moy.exe moysym.exe pshf.exe rcut.exe rpol.exe spin.exe spinsym.exe hdiab.exe 

Autocip13:
	$(MAKE) -e -C Autocip_13 
	
bd.exe:
	$(MAKE) -e -C CIPSI_sources/bd 
	
bdsym.exe:
	$(MAKE) -e -C CIPSI_sources/bdsym 	
	
cip.exe:
	$(MAKE) -e -C CIPSI_sources/cip 
	
cipsym.exe:
	$(MAKE) -e -C CIPSI_sources/cipsym 	
	
ciro.exe:
	$(MAKE) -e -C CIPSI_sources/ciro 
	
cval.exe:
	$(MAKE) -e -C CIPSI_sources/cval 	
	
fock.exe:
	$(MAKE) -e -C CIPSI_sources/fock 
	
ijkl.exe:
	$(MAKE) -e -C CIPSI_sources/ijkl 	
	
moy.exe:
	$(MAKE) -e -C CIPSI_sources/moy 
	
moysym.exe:
	$(MAKE) -e -C CIPSI_sources/moysym 	
	
pshf.exe:
	$(MAKE) -e -C CIPSI_sources/pshf 
	
rcut.exe:
	$(MAKE) -e -C CIPSI_sources/rcut 	
	
rpol.exe:	
	$(MAKE) -e -C CIPSI_sources/rpol 
	
spin.exe:	
	$(MAKE) -e -C CIPSI_sources/spin 
	
spinsym.exe:
	$(MAKE) -e -C CIPSI_sources/spinsym
	
hdiab.exe:	
	$(MAKE) -e -C CIPSI_sources/hdiab 
	
	

# bd.exe:	$(OBJSbd)
# 	$(COMP) $(OPT_LVL) $(LIB) -o bin/$@ $^ 
# 	
# bdsym.exe:	$(OBJSbdsym)
# 	$(COMP) $(OPT_LVL) $(LIB) -o bin/$@ $^ 
# 
# cip.exe:	$(OBJScip)
# 	$(COMP) $(OPT_LVL) $(LIB) -o bin/$@ $^ 
# 	
# cipsym.exe:	$(OBJScipsym)
# 	$(COMP) $(OPT_LVL) $(LIB) -o bin/$@ $^ 
# 	
# ciro.exe:	$(OBJSciro)
# 	$(COMP) $(OPT_LVL) $(LIB) -o bin/$@ $^ 
# 	
# cval.exe:	$(OBJScval)
# 	$(COMP) $(OPT_LVL) $(LIB) -o bin/$@ $^ 
# 	
# fock.exe:	$(OBJSfock)
# 	$(COMP) $(OPT_LVL) $(LIB) -o bin/$@ $^ 
# 	
# ijkl.exe:	$(OBJSijkl)
# 	$(COMP) $(OPT_LVL) $(LIB) -o bin/$@ $^ 
# 	
# moy.exe:	$(OBJSmoy)
# 	$(COMP) $(OPT_LVL) $(LIB) -o bin/$@ $^ 
# 	
# moysym.exe:	$(OBJSmoysym)
# 	$(COMP) $(OPT_LVL) $(LIB) -o bin/$@ $^ 
# 	
# pshf.exe:	$(OBJSpshf)
# 	$(COMP) $(OPT_LVL) $(LIB) -o bin/$@ $^ 
# 	
# rcut.exe:	$(OBJSrcut)
# 	$(COMP) $(OPT_LVL) $(LIB) -o bin/$@ $^ 
# 	
# rpol.exe:	$(OBJSrpol)
# 	$(COMP) $(OPT_LVL) $(LIB) -o bin/$@ $^ 
# 	
# spin.exe:	$(OBJSspin)
# 	$(COMP) $(OPT_LVL) $(LIB) -o bin/$@ $^ 
# 	
# spinsym.exe:	$(OBJSspinsym)
# 	$(COMP) $(OPT_LVL) $(LIB) -o bin/$@ $^ 
# 	
# hdiab.exe:	$(OBJShdiab)
# 	$(COMP) $(OPT_LVL) $(LIB) -o bin/$@ $^ 
# 	
.SUFFIXES: .o .f90 .f

.f90.o:
	$(COMP) $(OPT_LVL) -c $^ -o $@	
.f.o:
	$(COMP) $(OPT_LVL) -c $^ -o $@		
	
clean:
	rm -f bin/*.exe *.mod ${Foldauto}*.exe ${Foldbd}*.exe ${Foldbdsym}*.exe ${Foldcip}*.exe ${Foldcipsym}*.exe ${Foldciro}*.exe ${Foldcval}*.exe ${Foldfock}*.exe ${Foldijkl}*.exe ${Foldmoy}*.exe ${Foldmoysym}*.exe ${Foldpshf}*.exe ${Foldrcut}*.exe ${Foldrpol}*.exe ${Foldspin}*.exe ${Foldspinsym}*.exe ${Foldhdiab}*.exe 

cleano:
	rm -f $(OBJSauto) $(OBJSbd) $(OBJSbdsym) $(OBJScip) $(OBJScipsym) $(OBJSciro) $(OBJScval) $(OBJSfock) $(OBJSijkl) \
	$(OBJSmoy) $(OBJSmoysym) $(OBJSpshf) $(OBJSrcut) $(OBJSrpol) $(OBJSspin) $(OBJSspinsym) $(OBJShdiab) *.mod
