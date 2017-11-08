      subroutine nesbet(d,v,e,ia,in,ehond)
      implicit real*8(a-h,o-z)
      include 'pshf.prm'
      common/infoa/nat,ich,mul,numscf,nx,ne,na,nb
      common/pseudo/qn(doa)
      common/iofile/ir,iw,ip,is,iq,ih,iv
c     common/eigen/d(doas),v(doasq),e(doa),ia(doa),in(doa)
      dimension d(*),v(*),e(*),ia(*),in(*)
      dimension ns(10),log(4),xx(lsize),ix(lsize)
      dimension g(600),akij(15),tab(10,5,2),aisp(10),amult(10)
      data tab/0.75d0,9*0.d0,2*1.d0,8*0.d0,1.75d0,2*1.d0,7*0.d0,
     # 2.d0,4*1.d0,5*0.d0,2.75d0,6*1.d0,3*0.d0, 0.75d0,9*0.d0,2*1.d0,
     # 8*0.d0,1.d0,1.75d0,1.d0,7*0.d0,1.d0,2.d0,1.d0,1.d0,0.d0,
     # 1.d0,4*0.d0,1.d0,2.75d0,3*1.d0,0.d0,0.d0,1.d0,1.d0,0.d0/                 
c
      nxsq=numscf*numscf
      nsingl=0
c     -------------------------------------------------
c     y a-t-il des corrections a apporter a l'energie ?
c     -------------------------------------------------
      do 10 i=1,na
      if(qn(i).eq.2.d0) go to 10
      if(qn(i).ne.1.d0) go to 20
      nsingl=nsingl+1
      ns(nsingl)=i
   10 continue
      if(nsingl.eq.0) return
      if(nsingl.le.5) go to 30
      write(iw,9900) nsingl
 9900 format(/1x,'energy corrections are not included for ',i3,' singly
     #occupied orbitals')
      return
   20 write(iw,9890)
 9890 format(/,1x,'energy corrections are not included for non-integer
     #occupation numbers')
      return
   30 continue
c
c     --------------------------------------------
c     corrections pour une a cinq couches ouvertes
c     --------------------------------------------
c
c     -----------------
c     calcul des k(i,j)
c     -----------------
      ij=0
      do 35 i=1,nsingl
      n1=ns(i)
      do 35 j=1,numscf
      ij=ij+1
   35 g(ij)=v(in(n1)+j)
      ij=0
      do 36 i=1,nsingl
      do 36 j=1,i
      ij=ij+1
   36 akij(ij)=0.d0
      rewind is
   40 read(is) xx,ix
      do 60 m=1,lsize
      lable=ix(m)
      i=iword1(lable)
      j=iword2(lable)
      k=iword3(lable)
      l=iword4(lable)
      val=xx(m)
      if(lable.eq.0) go to 65
      ij=0
      do 50 i1=1,nsingl
      do 50 i2=1,i1
      ij=ij+1
      x1=g(in(i1)+i)*g(in(i2)+j)
      x2=g(in(i2)+i)*g(in(i1)+j)
      x3=g(in(i1)+k)*g(in(i2)+l)
      x4=g(in(i2)+k)*g(in(i1)+l)
      ak=(x1+x2)*(x3+x4)*val
      akij(ij)=akij(ij)+ak
   50 continue
   60 continue
      go to 40
   65 continue
      ij=0
      do 68 i=1,nsingl
      do 68 j=1,i
      ij=ij+1
   68 akij(ij)=akij(ij)*2.d0
c     ------------------------------------------------------------
c     calcul de la partie commune dans les elements diagonaux de h
c     ------------------------------------------------------------
      dum=ehond
      do 75 i=1,nsingl
      ij=ia(i)+i
      dum=dum-0.25d0*akij(ij)
      if(i.eq.1) go to 75
      im=i-1
      do 70 j=1,im
      ij=ia(i)+j
   70 dum=dum+0.5d0*akij(ij)
   75 continue
c     ----------------------------
c     construction de la matrice h
c     ----------------------------
      go to (110,120,130,140,150),nsingl
c
c     ------> une couche ouverte
c
  110 g(1)=dum
      ndet=1
      go to 200
c
c     ------> deux couches ouvertes
c
  120 g(1)=dum
      g(2)=-akij(2)
      g(3)=g(1)
      ndet=2
      go to 200
c
c     ------> trois couches ouvertes
c
  130 g(1)=dum-akij(2)
      g(2)=-akij(5)
      g(3)=dum-akij(4)
      g(4)=-akij(4)
      g(5)=-akij(2)
      g(6)=dum-akij(5)
      ndet=3
      go to 200
c
c     ------> quatre couches ouvertes
c
  140 g(1)=dum-akij(2)-akij(9)
      g(3)=dum-akij(4)-akij(8)
      g(6)=dum-akij(7)-akij(5)
      g(10)=g(6)
      g(15)=g(3)
      g(21)=g(1)
c
      g(2)=-akij(5)
      g(4)=-akij(4)
      g(5)=-akij(2)
      g(7)=-akij(8)
      g(8)=-akij(9)
      g(9)=0.d0
      g(11)=-akij(7)
      g(12)=0.d0
      g(13)=-akij(9)
      g(14)=-akij(2)
      g(16)=0.d0
      g(17)=-akij(7)
      g(18)=-akij(8)
      g(19)=-akij(4)
      g(20)=-akij(5)
c
      ndet=6
      go to 200
c
c
c     ------> cinq couches ouvertes
c
  150 g(1)=dum-akij(2)-akij(4)-akij(5)-akij(14)
      g(3)=dum-akij(2)-akij(7)-akij(8)-akij(13)
      g(6)=dum-akij(4)-akij(7)-akij(9)-akij(12)
      g(10)=dum-akij(5)-akij(8)-akij(9)-akij(11)
      g(15)=dum-akij(2)-akij(11)-akij(12)-akij(9)
      g(21)=dum-akij(4)-akij(11)-akij(13)-akij(8)
      g(28)=dum-akij(5)-akij(12)-akij(13)-akij(7)
      g(36)=dum-akij(7)-akij(11)-akij(14)-akij(5)
      g(45)=dum-akij(8)-akij(12)-akij(14)-akij(4)
      g(55)=dum-akij(9)-akij(13)-akij(14)-akij(2)
c
      g(2)=-akij(9)
      g(4)=-akij(8)
      g(5)=-akij(5)
      g(7)=-akij(7)
      g(8)=-akij(4)
      g(9)=-akij(2)
      g(11)=-akij(13)
      g(12)=-akij(14)
      g(13)=0.d0
      g(14)=0.d0
      g(16)=-akij(12)
      g(17)=0.d0
      g(18)=-akij(14)
      g(19)=0.d0
      g(20)=-akij(5)
      g(22)=-akij(11)
      g(23)=0.d0
      g(24)=0.d0
      g(25)=-akij(14)
      g(26)=-akij(4)
      g(27)=-akij(2)
      g(29)=0.d0
      g(30)=-akij(12)
      g(31)=-akij(13)
      g(32)=0.d0
      g(33)=-akij(8)
      g(34)=-akij(9)
      g(35)=0.d0
      g(37)=0.d0
      g(38)=-akij(11)
      g(39)=0.d0
      g(40)=-akij(13)
      g(41)=-akij(7)
      g(42)=0.d0
      g(43)=-akij(9)
      g(44)=-akij(2)
      g(46)=0.d0
      g(47)=0.d0
      g(48)=-akij(11)
      g(49)=-akij(12)
      g(50)=0.d0
      g(51)=-akij(7)
      g(52)=-akij(8)
      g(53)=-akij(4)
      g(54)=-akij(5)
c
      ndet=10
  200 continue
c
      ndd=ndet*(ndet+1)/2
      do 205 i=1,ndd
  205 d(i)=g(i)
c
c     -------------------------------
c     diagonalisation de la matrice h
c     -------------------------------
c
c     modification du tableau in utilise dans ligen
      do 210 i=1,ndet
  210 in(i)=(i-1)*ndet
c
      call ligen(ndet)
c
c     restauration du tableau in
      do 500 i=1,numscf
  500 in(i)=(i-1)*numscf
c
c
c     -----------------------------------
c     reconnaissance des vecteurs propres
c     -----------------------------------
      ij=0
      do 240 i=1,ndet
      im=1
      if(abs(v(ij+1)).gt.1.d-8) go to 222
      im=2
      if(abs(v(ij+2)).gt.1.d-8) go to 222
      go to 224
  222 continue
      dum=0.d0
      do 230 j=1,ndet
      ij=ij+1
      if(j.eq.im) cun=v(ij)
  230 dum=dum+tab(j,nsingl,im)*v(ij)
      scarre=dum/cun
      scarre=(-1.d0+sqrt(1.d0+4.d0*scarre))
      aisp(i)=0.5d0*scarre
      amult(i)=scarre+1.d0
      go to 240
  224 write(iw,9885) i
 9885 format(/,2x,'warning...s could not be computed for vector ',i2,//)
      aisp(i)=1.d+20
      amult(i)=1.d+20
  240 continue
c
c     ------------------------
c     impression des resultats
c     ------------------------
c
      write(iw,9880) (e(i),i=1,ndet)
      write(iw,9870)
      do 250 i=1,ndet
  250 write(iw,9860) (v((j-1)*ndet+i),j=1,ndet)
      write(iw,9850) (aisp(i),i=1,ndet)
      write(iw,9840) (amult(i),i=1,ndet)
      write(iw,9830)
c
 9880 format(//,10x,18(1h-),/,10x,'open-shell systems',/,10x,
     #18(1h-),//,1x,'eigenvalues',10f12.6)
 9870 format(//,1x,'eigenvectors')
 9860 format(12x,10f12.6)
 9850 format(//,6x,'s',5x,10f12.6)
 9840 format(//,4x,'2s+1',4x,10(f9.0,3x))
 9830 format(//)
c
      return
      end
