      subroutine rondof(norb,noc,enuc,yauto,group,ivecsy,nsym,ntrans,
     * naxis,nshell,nat,num,ne,yf)
      implicit real*8 (a-h,o-x,z),logical *1(y)
      include 'pshf.prm'
      character*4 iflab
      character*8 group
      common /hondo/ nt,ishell(doa),kloc(doa),ktyp(doa),
     & numc(doa),newsh(doa,nsymz),iptr(3,nsymz),idtr(6,nsymz),
     y iftr(10,nsymz),igtr(15,nsymz),yhondo,yprth
      common/nuscf/iscf(doa),iflab(3,doa)
      dimension ex(dgp),cf(dgp),cg(dgp),cs(dgp),cp(dgp),cd(dgp),
     1 c(3,dc),zan(dc),kstart(doa),iatno(dc),ks(dc),ns(dc),a(dc)
     & ,kng(doa),ptr(3,60),dtr(6,120),ftr(10,200),ktype(doa)
      dimension ivecsy(2*doa*nsymz),katom(doa)
      dimension iso(ds,12),invt(48)
      character*6 ante
      common /info/ante,yex
c
      rewind 2
      read(2) (idum,i=1,19),nt,nshell,num,norba,norbb,bname,
     1 ngauss,ne,ich,mul,na,nb,nat
      kp=3*nt
      kd=6*nt
      kf=10*nt
      kg=15*nt
      if(yprth) write(6,9999)
 9999 format(//' relecture des donnees definissant la base atomique')
c
      rewind 2
      read(2) (idum,i=1,56),
     y         (ktype(k),kstart(k),kloc(k),katom(k),kng(k),k=1,
     y nshell),(ex(k),cs(k),cp(k),cd(k),cf(k),cg(k),k=1,ngauss),
     y ((c(i,k),i=1,3),a(k),ns(k),ks(k),iatno(k),zan(k),k=1,nat),
     y ((newsh(i,k),i=1,nshell),k=1,nt),((ptr(i,k),i=1,3),k=1,kp),
     y  ((dtr(i,k),i=1,6),k=1,kd),((ftr(i,k),i=1,10),k=1,kf)
c    y ((gtr(i,k),i=1,15),k=1,kg)
     y , ((iflab(i,k),i=1,3),k=1,num),
     y ((iso(i,j),i=1,nshell),j=1,12),invt,enuc,ante
      ante='  pshf'
      do 2550 i=1,nshell
      if(yprth)  write(6,9998) i,ktype(i),katom(i)
      i1=kstart(i)
      i2=i1+kng(i)-1
      if(yprth)write(6,9997)(ex(l),cs(l),cp(l),cd(l),cf(l),cg(l),
     1 l=i1,i2)
 9998 format(/' couche',i3,' type',i3,' atome',i3,'  exposant     cs
     y     cp          cd          cf          cg')
 9997 format(24x,6f12.6)
      kty=ktype(i)
c
      go to (9001,9002,9003,9004,9005,9005,9005,9008,9009,9011),kty
 9001 nfonct=1
      go to 9010
 9002 nfonct=3
      go to 9010
 9003 nfonct=5
      if(ante.eq.'gamess') nfonct=6
       go to 9010
 9004 nfonct=7
      if(ante.eq.'gamess') nfonct=10
      go to 9010
 9005 nfonct=10
      if(ante.eq.'gamess') nfonct=15
      go to 9010
 9008 nfonct=6
      ktype(i)=3
      yex=.true.
      go to 9010
 9009 nfonct=10
      ktype(i)=4
      yex=.true.
      go to 9010
 9011 nfonct=15
      ktype(i)=5
      yex=.true.
 9010 continue
      jmin=kloc(i)
      jmax=kloc(i)+nfonct-1
      do 2545 j=jmin,jmax
      numc(j)=j-jmin
      ktyp(j)=ktype(i)
 2545 ishell(j)=i
2550  continue
      if(yprth) write(6,9996)
      if(yprth) write(6,9995) (a(i),(c(j,i),j=1,3),i=1,nat)
 9996 format(//' coordonnees des atomes')
 9995 format(/,a10,3f12.6)
      ndmax=5
      nfmax=7
      ngmax=10
      if(ante.eq.'gamess'.or.yex) then
	 ndmax=6
	 nfmax=10
	 ngmax=15
      endif
      k2p=0
      k2f=0
      k2d=0
      k2g=0
      do 2570 j=1,nt
c
      k1p=k2p+1
      k1d=k2d+1
      k1f=k2f+1
      k1g=k2g+1
      k2p=k2p+3
      k2d=k2d+6
      k2f=k2f+10
      k2g=k2g+15
      do 2565 i=1,3
2565  iptr(i,j)=ptr(i,k1p+i-1)
       do 2569 i=1,ndmax
2569  idtr(i,j)=nint(dtr(i,k1d+i-1))
      do 2567 i=1,nfmax
 2567 iftr(i,j)=nint(ftr(i,k1f+i-1))
c     do 2568 i=1,ngmax
c2568 igtr(i,j)=nint(gtr(i,k1g+i-1))
2570  continue
      write(6,*)((idtr(i,j),i=1,ndmax),j=1,nt)
62    format(1x,40i3)
      norb=num
      noc=(ne+1)/2
      if(.not.yauto) go to 200
      call otodon(group,ntrans,nshell,c,naxis,norb,ktype,katom,yf)
  200 continue
      return
      end
