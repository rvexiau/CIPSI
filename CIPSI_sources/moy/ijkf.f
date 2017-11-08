      subroutine ijkf
      implicit real*8 (a-h,o-x,z),logical*1 (y)
      include 'pshf.prm'
c---------lecture des integrales j, k et f
      integer*4 iorb,iwspin,itsym,its,isytr
      common/spo/ feff(2*doa),fdiag(2*doa),
     * iorb(2*doa),iwspin(2*doa),ispin(2*doa),
     *itsym(2*doa),its(nsymz,nsymz),yoc(2*doa),isytr(2*doa,nsymz)
      common  /cipsi/ norb,mnorb,noca,nocb,
     *igela,igelb,ncf,ncper,nsym,isym,isz,ntrsy,metat,
     *yprt,yion
      common/ic/c(ndetz*metz),e(ndetz),eb(metz),e2(metz),
     *emp(metz),e2en(metz),teste(metz),elevy(metz),e2mp(metz),ac,ei,
     *hdef,dpivot,fmpb(2*doa),dmpb(10),ishif,ietat
      common/jkf/ f(doa*(doa+1)/2),aj(doa*(doa+1)/2),
     *ak(doa*(doa+1)/2),empref(metz),tocom(2*doa)
      num(i)=i*(i-1)/2
      nter=norb*(norb+1)/2
      read(40) dum
      read(40) (aj(l),l=1,nter)
      read(40) (ak(l),l=1,nter)
      read(40) (f (l),l=1,nter)
      l=0
      do 3200 ll=1,norb
      l=l+ll
      fdiag(ll)=f(l)
3200  fdiag(ll+norb)=f(l)
      if(.not.yion) go to 3600
      nter=nter+1
      nn=nter+norb+norb+2
      do 3500 l=nter,nn
      f(l)=0.
c
c
      aj(l)=0.
3500  ak(l)=0.
      fdiag(mnorb+1)=0.
      feff(mnorb+1)=0.
      fdiag(mnorb+2)=0
      feff(mnorb+2)=0
3600  continue
      return
      end
