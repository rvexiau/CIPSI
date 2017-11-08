      subroutine buildz
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      integer*2 symbl,namel(100)
      character*1 name(100,8),aster,blank
     *,symb,figure,minus,plus,equal
      logical*1 rx,ry,rz
      common /zmat/ iz(50,4),bl(50),alph(50),bet(50),nz,nzz,
     * ipar(15,5),nipar(5),npar,nstep,xx(5),dx(5),num,idum,confac,
     * name,namel,rx,ry,rz
      common/infoa/nat,ich,mul,numscf,nx,ne,na,nb,zan(dc),c(3,dc)
      dimension alpha(50),beta(50)
      dimension a(50),b(50),d(50),u1(3),u2(3),u3(3),u4(3),
     * vj(3),vp(3),v3(3)
      dimension dis(50,50)
      dimension iel(18)
      common/geod/ value(100),nsymb,symb(100,8),symbl(100)
     *,figure(16),minus,plus,equal,aster,blank
      common/opti7/ww,fin
      logical*1 ww,fin
      zero=0.0d0
      one=1.d0
      two=2.d0
      tenm5=1.0d-5
      tenm6=1.d-6
      pi=3.141592653589793d0
 1000 format(///5x,'conversion factor =',f15.8////
     * 15x,'coordinates'/20x,'x',14x,'y',14x,'z')
 1010 format(2x,9a1,3x,3f15.8)
 1020 format(//30x,'z matrix'/4x,'centre',5x,'z1',4x,'bl',11x,'z2',4x
     *'alpha',8x,'z3',4x,'beta',9x,'z4',/)
 1030 format(2x,9a1)
 1040 format(2x,9a1,i5,g14.7,'(',i3,')   ')
 1050 format(2x,9a1,i5,g14.7,'(',i3,')   ',i5,g14.7,'(',i3,')   ')
 1060 format(2x,9a1,i5,g14.7,'(',i3,')   ',i5,g14.7,'(',i3,')   ',
     *i5,g14.7,'(',i3,')   ',i5)
 1070 format(5x,'x=',f18.12,',y=',f18.12,',z=',f18.12,',')
 1080 format(4a1,3f20.10)
      torad=pi/180.d0
c
c     print z matrix
c
      if(ww)write(6,1020)
      lng=namel(1)
      if(ww)write(6,1030) (name(1,j),j=1,lng)
      if (nz.le.1) go to 6
      i1=1
      lng=namel(2)
      nbl=9-lng
      if(ww)write(6,1040) (name(2,j),j=1,lng),(blank,j=1,nbl),iz(2,1),
     *bl(2),i
      if(nz.le.2) go to 6
      i1=i1+1
      i2=nz
      lng=namel(3)
      nbl=9-lng
      if(ww)write(6,1050) (name(3,j),j=1,lng),(blank,j=1,nbl),
     *               iz(3,1),bl(3),i1,iz(3,2),alph(3),i2
      if(nz.le.3) go to 6
      i3=2*nz-3
      do 4 i=4,nz
      i1=i1+1
      i2=i2+1
      i3=i3+1
      lng=namel(i)
      nbl=9-lng
      if(ww)write(6,1060) (name(i,j),j=1,lng),(blank,j=1,nbl),
     *               iz(i,1),bl(i),i1,iz(i,2),alph(i),i2,
     *               iz(i,3),bet(i),i3,iz(i,4)
    4 continue
c     zero coordinate array c
    6 do 10 i=1,nz
      do 10 j=1,3
   10 c(j,i)=zero
c     convert angles from degrees to radians
      do 20 i=1,nz
      alpha(i)=alph(i)*torad
   20 beta(i)=bet(i)*torad
      c(3,2)=bl(2)
      if(nz-3)260,30,30
   30 c(1,3)=bl(3)*sin(alpha(3))
      if(iz(3,1)-1)50,40,50
   40 c(3,3)=bl(3)*cos(alpha(3))
      go to 60
   50 c(3,3)=c(3,2)-bl(3)*cos(alpha(3))
   60 do 80 i=4,nz
      if(abs(c(1,i-1))-tenm5)70,90,90
   70 c(1,i)=bl(i)*sin(alpha(i))
      itemp=iz(i,1)
      jtemp=iz(i,2)
   80 c(3,i)=c(3,itemp)-bl(i)*cos(alpha(i))*sign(one,c(3,itemp)
     *-c(3,jtemp))
   90 k=i
      if(k-nz)100,100,260
  100 do 250 j=k,nz
      dcaj=cos(alpha(j))
      dsaj=sin(alpha(j))
      dcbj=cos(beta(j))
      dsbj=sin(beta(j))
      if(iz(j,4))135,110,135
  110 call vec(u1,c,iz(j,2),iz(j,3))
      call vec(u2,c,iz(j,1),iz(j,2))
      call vprod(vp,u1,u2)
      r=sqrt(one-(u1(1)*u2(1)+u1(2)*u2(2)+u1(3)*u2(3))**2)
      do 120 i=1,3
  120 u3(i)=vp(i)/r
      call vprod(u4,u3,u2)
      do 130 i=1,3
      vj(i)=bl(j)*(-u2(i)*dcaj+u4(i)*dsaj*dcbj+u3(i)*dsaj*dsbj)
      itemp=iz(j,1)
  130 c(i,j)=vj(i)+c(i,itemp)
      go to 250
  135 if(iabs(iz(j,4))-1)210,140,210
  140 call vec(u1,c,iz(j,1),iz(j,3))
      call vec(u2,c,iz(j,2),iz(j,1))
      zeta=-(u1(1)*u2(1)+u1(2)*u2(2)+u1(3)*u2(3))
      a(j)=(-dcbj+zeta*dcaj)/(one-zeta*zeta)
      b(j)=(dcaj-zeta*dcbj)/(one-zeta*zeta)
      r=zero
      gamma=pi/two
      if(zeta)150,170,160
  150 r=pi
  160 gamma=atan(sqrt(one-zeta*zeta)/zeta)+r
  170 d(j)=zero
      if(abs(gamma+alpha(j)+beta(j)-two*pi)-tenm6)190,180,180
  180 d(j)=iz(j,4)*(sqrt(one+a(j)*dcbj-b(j)*dcaj))/sqrt(one-zeta*zeta)
  190 call vprod(v3,u1,u2)
      do 200 i=1,3
      u3(i)=a(j)*u1(i)+b(j)*u2(i)+d(j)*v3(i)
      vj(i)=bl(j)*u3(i)
      itemp=iz(j,1)
  200 c(i,j)=vj(i)+c(i,itemp)
      go to 250
  210 call vec(u1,c,iz(j,1),iz(j,3))
      call vec(u2,c,iz(j,2),iz(j,1))
      zeta=-(u1(1)*u2(1)+u1(2)*u2(2)+u1(3)*u2(3))
      call vprod(v3,u1,u2)
      v3mag=sqrt(v3(1)*v3(1)+v3(2)*v3(2)+v3(3)*v3(3))
      a(j)=v3mag*dcbj/(one-zeta*zeta)
      b(j)=sqrt((one-dcaj*dcaj-a(j)*dcbj*v3mag)/(one-zeta*zeta))
      if(iz(j,4)-2)220,230,220
  220 b(j)=-b(j)
  230 d(j)=b(j)*zeta+dcaj
      do 240 i=1,3
      u3(i)=b(j)*u1(i)+d(j)*u2(i)+a(j)*v3(i)
      vj(i)=bl(j)*u3(i)
      itemp=iz(j,1)
  240 c(i,j)=vj(i)+c(i,itemp)
  250 continue
  260 continue
c
c     print coordinates
c
      natoms=nz
      ione=1
      if(ww)write(6,1000) confac
      ik=0
      do 320 ii=1,natoms
      i=ii-ik
      if(ik.ne.0) then
      do 305 ij=1,3
 305  c(ij,i)=c(ij,ii)
      endif
      if (name(ii,1).eq.aster) go to 300
      lng=namel(ii)
      nbl=9-lng
      c(1,i)=confac*c(1,i)
      c(2,i)=confac*c(2,i)
      c(3,i)=confac*c(3,i)
      if (abs(c(1,i)).lt.1.d-12) c(1,i)=0.d0
      if (abs(c(2,i)).lt.1.d-12) c(2,i)=0.d0
      if (abs(c(3,i)).lt.1.d-12) c(3,i)=0.d0
      if (.not.rx) go to 270
      cy=c(2,i)
      c(2,i)=c(3,i)
      c(3,i)=-cy
  270 if (.not.ry) go to 280
      cx=c(1,i)
      c(1,i)=-c(3,i)
      c(3,i)=cx
  280 if (.not.rz) go to 290
      cx=c(1,i)
      c(1,i)=c(2,i)
      c(2,i)=-cx
  290 continue
      if(ww)write(6,1010) (name(ii,j),j=1,lng),(blank,j=1,nbl),
     *(c(j,i),j=1,3)
      if (lng.gt.4) lng=4
      nbl=4-lng
      goto 320
  300 ik=ik+1
 320  continue
      return
      end
