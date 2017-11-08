      subroutine mulken(scftyp)
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      logical*1 last,out
      common/iofile/ir,iw,ip,is,iq,ih,iv,if2,if3,if4,ioda(19)
      common/output/nprint,itol,icut,normf,normp,nopk
      common/runlab/iflab(3,doa)
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,zan(dc),c(3,dc)
c     common/hfden/da(5050),db(5050),ia(100)
      common/hfden/da(doas),db(doas),ia(doa)
      common/atlim/limlow(dc),limsup(dc)
c     dimension s(5050),d(5050,2),goc(100,2)
      dimension s(doas),d(doas,2),goc(doa,2)
      dimension e(dc*(dc+1)/2),gac(dc,2)
      equivalence (da(1),d(1,1))
      data alpha,beta,all /8h   alpha,8h    beta,8h     all/
      data open /8huhf     /
 9999 format(//,10x,28(1h-),/,10x,'mulliken population analysis',
     1 10x,a8,' electrons',/,10x,28(1h-))
 9998 format(/,10x,'----- total gross population in ao"s -----',/)
 9997 format(10x,i5,2x,3a4,f12.5)
 9996 format(/,10x,'----- condensed to atoms -----',/)
 9995 format(/,10x,'----- total gross population on atoms ----',/)
 9994 format(/,10x,'----- ao"s spin population -----',/)
 9993 format(/,10x,'----- atomic spin population -----',/)
 9992 format(10x,i5,2x,2a4,2x,f6.1,f12.5)
      out=nprint.eq.3
c
c     ----- determine number of basis functions per atom -----
c
      call aolim
c
c     ----- read in overlap matrix -----
c
      nav=ioda(2)
      call reada(s,nx,2)
      ipass=1
      last=.false.
      elec=alpha
      if(scftyp.ne.open) elec=all
  100 continue
c
c     ----- do a mulliken population analysis ----
c           calculate overlap population
c
      write(6,*)'nx avant ovlpop',nx
      call ovlpop(d(1,ipass),s,nx)
  150 write(iw,9999) elec
      if(out) call matout(d(1,ipass),num)
c
c     ----- calculate total gross population in ao's ----
c
      call grossc(d(1,ipass),goc(1,ipass),ia,num)
      write(iw,9998)
      write(iw,9997) (i,(iflab(k,i),k=1,3),goc(i,ipass),i=1,num)
c
c     ----- compress from orbitals to atoms -----
c
      call atpop(d(1,ipass),ia,e,nat)
      write(iw,9996)
      if(out) call matout(e,nat)
c
c     ----- calculate total gross population on atoms -----
c
      call grossc(e,gac(1,ipass),ia,nat)
      write(iw,9995)
      write(iw,9992)(i,iflab(1,limlow(i)),iflab(2,limlow(i)),
     * zan(i),gac(i,ipass),i=1,nat)
      if(scftyp.ne.open) return
      if(last) return
      if(ipass.eq.2) go to 200
      ipass=2
      elec=beta
      go to 100
  200 continue
c
c     ----- calculate orbital and atomic spin densities -----
c
      do 300 i=1,num
  300 goc(i,1)=goc(i,1)-goc(i,2)
      do 400 i=1,nat
  400 gac(i,1)=gac(i,1)-gac(i,2)
      write(iw,9994)
      write(iw,9997) (i,(iflab(k,i),k=1,3),goc(i,1),i=1,num)
      write(iw,9993)
      write(iw,9992)(i,iflab(1,limlow(i)),iflab(2,limlow(i)),
     *  zan(i),gac(i,1),i=1,nat)
c
c     ----- do all electrons -----
c
      ipass=1
      last=.true.
      elec=all
      do 500 i=1,nx
  500 d(i,1)=d(i,1)+d(i,2)
      go to 150
      end
