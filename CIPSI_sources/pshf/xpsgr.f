      subroutine xpsgr
      implicit real*8(a-h,o-z)
      include 'pshf.prm'
      logical egp
      character*40 typegp,typf21
      common/psgrx/negp,egp,typegp(dc)
      common/output/nprint
      common/eigen2/ f(doas),v(doa,doa)
      common/iofile/ ir,iw,ip,is
      common/infoa/nat,ich,mul,numscf,nnp,ne,na,nb,zan(dc),c(3,dc)
      common/nshel/ ex(dgp),cs(dgp),cp(dgp),cd(dgp),cf(dgp),cg(dgp),
     1kstart(ds),katom(ds),ktype(ds),kng(ds),kloc(ds),kmin(ds),kmax(ds),
     2nshell,kmaz(ds)
      dimension vv(doasq),vtemp(doa,doa),ff(doa,doa),temp(doa)
      dimension exp(dgp),csp(dgp),cpp(dgp),cdp(dgp),cfp(dgp),cgp(dgp),
     1 xcp(3,dc),iatnp(dc),zap(dc),titlp(10),
     1 kstarp(ds),katop(ds),ktypp(ds),knp(ds),klop(ds)
      equivalence (vv(1),ff(1,1))
      if(.not.egp)return
      call psgrin
      itol=9
      icut=9
c
      nxscf=(numscf*(numscf+1))/2
c     boucle sur les groupes
      do   igrp=1,negp
        rewind 21
    1   read(21)typf21
        if(typf21.ne.typegp(igrp)) then
          read(21)
          read(21)
          go to 1
        else
          read(21) (idum,i=1,20),nshelp,nup,idum,idum,chdum,ngp,nep,
     1    (idum,i=1,4),natp,chdum,idum,titlp,
     2    (ktypp(k),kstarp(k),klop(k),katop(k),knp(k),k=1,nshelp),
     3    (exp(k),csp(k),cpp(k),cdp(k),cfp(k),cgp(k),k=1,ngp),
     4    ((xcp(i,k),i=1,3),xdum,idum,idum,iatnp(k),zap(k),k=1,natp)
	  write(6,*)' type de pseudo',typegp(igrp)
	  write(6,*)'nshelp,nup,ngp',nshelp,nup,ngp
          write(6,*)(ktypp(k),kstarp(k),klop(k),katop(k),
     *	  knp(k),k=1,nshelp)
          write(6,*)(exp(k),csp(k),cpp(k),cdp(k),cfp(k),
     *	  cgp(k),k=1,ngp)
	  ngaus=nup
	  ngaux=ngaus*(ngaus+1)/2
	  nbas=numscf
          call  sjpd(nshell,num,katom,c,kstart,kng,ktype,kloc,
     *              ex,cs,cp,cd,cf,cg,
     *              nshelp,nup,katop,xcp,kstarp,knp,ktypp,klop,
     *              exp,csp,cpp,cdp,cfp,cgp,
     *              itol,icut,vv)
          call trvm(vv,v,doasq,doa,doa,1,ngaus,1,nbas,.false.)
          read(21) (f(i),i=1,ngaux)
	  write(6,*)' coefficients du pseudo'
	  write(6,*) (f(i),i=1,ngaux)
	  do i=1,ngaus
	  f(i*(i+1)/2)=f(i*(i+1)/2)*2.
	  end do
          call trvm(f,vtemp,doas,doa,doa,1,ngaus,1,ngaus,.true.)
c....     calcul de l'operateur partiel dans la base atomique
          call matbc(v,doa,doa,1,ngaus,1,nbas,
     <               vtemp,doa,doa,1,ngaus,1,ngaus,
     <               v,doa,doa,1,ngaus,1,nbas,
     <               ff,doa,doa,1,nbas,1,nbas,
     <               temp)
	  call trmv(ff,f,doa,doa,doas,1,nbas,1,nbas,.true.)
	  write(6,*)' matrice de l''operateur'
	  call fout(f,nbas)
c......  on ajoute la matrice egp a l'operateur monoelectronique
          call reada(vv,nxscf,3)
          do i=1,nxscf
	    f(i)=f(i)+vv(i)
	  end do
          call wrtda(f,nxscf,3)
        end if
      end do
c
      return
      end
