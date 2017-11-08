      subroutine ptgrp
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      logical stand
      character*8 drcmin,drcmaj,grpmin,grpmaj,group,direct
      character*8 title
      common/output/iop1,itol,icut,normf,normp
      common/iofile/ir,iw,ip,is
      common/symtry/invt(48),iso(ds,12),ptr(3,144),dtr(6,288)
     1,ftr(10,480),nt
      common/transf/xsmal,ysmal,zsmal,xnew,ynew,znew,xp,yp,zp
      common/frame/u1,u2,u3,v1,v2,v3,w1,w2,w3,x0,y0,z0,x1,y1,z1,
     1 x2,y2,z2,direct,group,title(10),naxis
      common/symmat/t(432)
      dimension grpmin(19),grpmaj(19),drcmin(2),drcmaj(2)
      data grpmaj/'C1','CS','CI','CN','S2N','CNH','CNV','DN','DNH',
     1 'DND','CINFV','DINFH','T','TH','TD','O','OH','I','IH'/
      data grpmin/'c1','cs','ci','cn','s2n','cnh','cnv','dn','dnh',
     1 'dnd','cinfv','dinfh','t','th','td','o','oh','i','ih'/
c
      data drcmin /'normal','parall'/
      data drcmaj /'NORMAL','PARALL'/
      data tol /1.0d-10/
      data pi2 /6.28318530717958d+00/
 9999 format (a5,i5)
 9998 format(' are you kidding... you should give up...')
 9997 format(9f10.5)
 9996 format(3f10.5,a8)
 9995 format(' linear molecule , point group is cinfv or dinfh ',/,
     1 ' group ',a8, '  will be replaced by group ',a8,' during
     2  scf calculation')
 9993 format(10a8)
 9991 format(/,' the point group of the molecule is ...',a8,/,
     1 ' the order of the principal axis is ...',i5)
 9990 format(/,' the origin of the local frame is at x =  ',f10.5,
     1 ' y = ',f10.5,' z = ',f10.5,/,
     1 ' director cosines of the new axes ',/,34x,3(5x,f10.5),/,34x,
     1 3(5x,f10.5),/,34x,3(5x,f10.5))
 9989 format(' rotations about principal axis')
 9988 format(' sigma-h followed by rotations')
 9987 format(' c2 followed by rotations ')
 9986 format(' sigma-d followed by rotations')
 9985 format(' sigma-v followed by rotations')
 9984 format(/,10x,' center of symmetry at x = ',f10.5,' y = ',f10.5,
     1 ' z = ',f10.5)
 9983 format(/,' plane of symmetry defined by its normal u = ',f10.5,
     1 ' v = ',f10.5,' w = ',f10.5)
 9982 format(/,10x,3f15.9,/,10x,3f15.9,/,10x,3f15.9)
 9981 format(' c2 followed by sigma-h followed by rotations')
 9980 format(' sigma-d followed by c2 followed by rotations')
 9979 format(' s2n rotation followed by rotations')
      stand=.true.
c     read(ir,9999) group,naxis
      do 2 i=1,19
    2 if(group.eq.grpmin(i).or.group.eq.grpmaj(i)) index=i
      write(iw,9991) grpmin(index),naxis
      if(index.ne.18.and.index.ne.19) go to 5
      write(iw,9998)
      stop
    5 rho=sqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)
      if(rho.gt.tol) stand=.false.
      if(index.le.3) go to 200
      if(index.eq.11.or.index.eq.12) go to 200
c
c     define local frame
c     read in principal axis   ( 1 card )
c     read in x-local axis   ( 1 card )
c     default option _ local frame identical to master frame
c
c     read(ir,9997) x0,y0,z0,x1,y1,z1
      if(.not.stand) go to 7
c
      x0=0.0d+00
      y0=0.0d+00
      z0=0.0d+00
      x1=0.0d+00
      y1=0.0d+00
      y2=0.0d+00
      z2=0.0d+00
      z1=1.0d+00
      x2=1.0d+00
      direct=drcmin(2)
      rho=1.0d+00
      go to 7
c   6 read(ir,9996) x2,y2,z2,direct
    7 continue
      w1=(x1-x0)/rho
      w2=(y1-y0)/rho
      w3=(z1-z0)/rho
      ww=w1*w1+w2*w2+w3*w3
      x02=x2-x0
      y02=y2-y0
      z02=z2-z0
      rho=(w1*x02+w2*y02+w3*z02)/ww
      dum=rho*w1
      x0=x0+dum
      x02=x02-dum
      dum=rho*w2
      y0=y0+dum
      y02=y02-dum
      dum=rho*w3
      z0=z0+dum
      z02=z02-dum
      uu=(x02*x02+y02*y02+z02*z02)
      u=sqrt(uu)
      u1=x02/u
      u2=y02/u
      u3=z02/u
      v3=w1*u2-w2*u1
      v2=w3*u1-w1*u3
      v1=w2*u3-w3*u2
      if(direct.eq.drcmin(2).or.direct.eq.drcmaj(2)) go to 8
      dum=u1
      u1=v1
      v1=-dum
      dum=u2
      u2=v2
      v2=-dum
      dum=u3
      u3=v3
      v3=-dum
    8 continue
      if(iop1.ne.3) go to 9
      write(iw,9990) x0,y0,z0,u1,v1,w1,u2,v2,w2,u3,v3,w3
    9 continue
      if(index.ge.13) go to 200
c
c     rotation about principal axis
c
      nn=0
      n=naxis
      alpha=0.0d+00
      alph=pi2/dble(n)
   10 nn=nn+1
      if(nn.gt.n) go to 20
      cosa=cos(alpha)
      sina=sin(alpha)
      i=9*(nn-1)
      t(i+1)=cosa
      t(i+5)=cosa
      t(i+2)=-sina
      t(i+4)=sina
      t(i+3)=0.0d+00
      t(i+6)=0.0d+00
      t(i+7)=0.0d+00
      t(i+8)=0.0d+00
      t(i+9)=1.0d+00
      alpha=alpha+alph
      go to 10
c
c     end of group 4
c
   20 nt=n
      ii=9*nt
      if(iop1.ne.3) go to 24
      write(iw,9989)
      n1=1
      n2=naxis
      call print(n1,n2)
   24 continue
      if(index.eq.4) go to 1000
      if(index.eq.5) go to 500
      if(index.eq.7) go to 115
      if(index.ne.6.and.index.ne.9) go to 55
c
c     sigma-h plane  equation (z=0) in local frame
c
      nn=0
   30 nn=nn+1
      if(nn.gt.nt) go to 50
c
c     group 6 0r 9
c
      i=ii+9*(nn-1)
      do 40 j=1,8
   40 t(i+j)=t(i+j-ii)
      t(i+9)=-t(i+9-ii)
      go to 30
   50 nt=nt+nt
      ii=9*nt
      if(iop1.ne.3) go to 54
      write(iw,9988)
      n1=n2+1
      n2=n2+naxis
      call print(n1,n2)
   54 continue
c
c     end of group 6
c
      if(index.eq.6) go to 1000
c
c     one cp2 axis is the x-axis of the local frame
c     group 8 , 9 ,10
c
   55 continue
      nn=0
   60 nn=nn+1
      if(nn.gt.nt) go to 70
      i=ii+9*(nn-1)
      t(i+1)=t(i+1-ii)
      t(i+2)=-t(i+2-ii)
      t(i+3)=-t(i+3-ii)
      t(i+4)=t(i+4-ii)
      t(i+5)=-t(i+5-ii)
      t(i+6)=-t(i+6-ii)
      t(i+7)=t(i+7-ii)
      t(i+8)=-t(i+8-ii)
      t(i+9)=-t(i+9-ii)
      go to 60
   70 nt=nt+nt
      ii=9*nt
      if(iop1.ne.3) go to 99
      write(iw,9987)
      n1=n2+1
      n2=n2+naxis
      call print(n1,n2)
      if(index.ne.9) go to 99
      write(iw,9981)
      n1=n2+1
      n2=n2+naxis
      call print(n1,n2)
   99 continue
c
c     end of group 8 and 9
c
      if(index.eq.8.or.index.eq.9) go to 1000
c
c     dnd group . equation of plane sigma-d is _
c     sin(alph/4)*x-cos(alph/4)*y=0
c     the x-axis is the cp2 axis.
c
c     group 10
c
      beta=0.5d+00*alph
      cosa=cos(beta)
      sina=sin(beta)
      nn=0
  100 nn=nn+1
      if(nn.gt.nt) go to 110
      i=ii+9*(nn-1)
      t(i+1)=cosa*t(i+1-ii) + sina*t(i+2-ii)
      t(i+2)=sina*t(i+1-ii) - cosa*t(i+2-ii)
      t(i+3)=     t(i+3-ii)
      t(i+4)=cosa*t(i+4-ii) + sina*t(i+5-ii)
      t(i+5)=sina*t(i+4-ii) - cosa*t(i+5-ii)
      t(i+6)=     t(i+6-ii)
      t(i+7)=cosa*t(i+7-ii) + sina*t(i+8-ii)
      t(i+8)=sina*t(i+7-ii) - cosa*t(i+8-ii)
      t(i+9)=     t(i+9-ii)
      go to 100
  110 nt=nt+nt
      ii=9*nt
      if(iop1.ne.3) go to 114
      write(iw,9986)
      n1=n2+1
      n2=n2+naxis
      call print(n1,n2)
      write(iw,9980)
      n1=n2+1
      n2=n2+naxis
      call print(n1,n2)
  114 continue
c
c     end of group 10
c
      go to 1000
c
c     group 7
c     sigma-v is the (x-z) plane of local frame
c
  115 continue
      nn=0
  120 nn=nn+1
      if(nn.gt.nt) go to 130
      i=ii+9*(nn-1)
      t(i+1)=t(i+1-ii)
      t(i+2)=-t(i+2-ii)
      t(i+3)=t(i+3-ii)
      t(i+4)=t(i+4-ii)
      t(i+5)=-t(i+5-ii)
      t(i+6)=t(i+6-ii)
      t(i+7)=t(i+7-ii)
      t(i+8)=-t(i+8-ii)
      t(i+9)=t(i+9-ii)
      go to 120
  130 nt=nt+nt
      ii=9*nt
c
c     end of group 7
c
      if(iop1.ne.3) go to 1000
      write(iw,9985)
      n1=n2+1
      n2=n2+naxis
      call print(n1,n2)
      go to 1000
  200 continue
      t(1)=1.0d+00
      t(5)=1.0d+00
      t(9)=1.0d+00
      t(2)=0.0d+00
      t(3)=0.0d+00
      t(4)=0.0d+00
      t(6)=0.0d+00
      t(7)=0.0d+00
      t(8)=0.0d+00
      if(index.eq.1) go to 210
      if(index.eq.2) go to 250
      if(index.eq.3) go to 300
      if(index.eq.11.or.index.eq.12) go to 400
      go to 600
  210 nt=1
      x0=0.0d+00
      y0=0.0d+00
      z0=0.0d+00
      u1=1.0d+00
      v2=1.0d+00
      w3=1.0d+00
      u2=0.0d+00
      u3=0.0d+00
      v1=0.0d+00
      v3=0.0d+00
      w1=0.0d+00
      w2=0.0d+00
      go to 1000
c
c     cs symmetry group
c     the 3 points 1,2,3 define sigma-h plane
c
  250 continue
      if(.not.stand) go to 251
c
c     default option _ plane is the (x,y) plane
      x1=0.0d+00
      y1=0.0d+00
      z1=0.0d+00
      y2=0.0d+00
      z2=0.0d+00
      x3=0.0d+00
      z3=0.0d+00
      x2=1.0d+00
      y3=1.0d+00
      go to 252
  251 x3=x2
      y3=y2
      z3=z2
      x2=x1
      y2=y1
      z2=z1
      x1=x0
      y1=y0
      z1=z0
      x0=0.d0
      y0=0.d0
      z0=0.d0
  252 continue
      nt=2
      w1=(y2-y1)*(z3-z1)-(y3-y1)*(z2-z1)
      w2=(z2-z1)*(x3-x1)-(z3-z1)*(x2-x1)
      w3=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)
      rho=sqrt(w1*w1+w2*w2+w3*w3)
      w1=w1/rho
      w2=w2/rho
      w3=w3/rho
      u1=x2-x1
      u2=y2-y1
      u3=z2-z1
      rho=sqrt(u1*u1+u2*u2+u3*u3)
      u1=u1/rho
      u2=u2/rho
      u3=u3/rho
      v1=w2*u3-w3*u2
      v2=w3*u1-w1*u3
      v3=w1*u2-w2*u1
      x0=x1
      y0=y1
      z0=z1
      t(10)=1.0d+00
      t(14)=1.0d+00
      t(18)=-1.0d+00
      t(11)=0.0d+00
      t(12)=0.0d+00
      t(13)=0.0d+00
      t(15)=0.0d+00
      t(16)=0.0d+00
      t(17)=0.0d+00
      if(iop1.ne.3) go to 1000
      write(iw,9983) w1,w2,w3
      write(iw,9982) u1,v1,w1,u2,v2,w2,u3,v3,w3
      go to 1000
c
c     ci symmetry group
c     center of inversion is (x0,y0,z0)
c
  300 continue
      if(iop1.eq.3) write(iw,9984) x0,y0,z0
      t(10)=-1.0d+00
      t(14)=-1.0d+00
      t(18)=-1.0d+00
      t(11)=0.0d+00
      t(12)=0.0d+00
      t(13)=0.0d+00
      t(15)=0.0d+00
      t(16)=0.0d+00
      t(17)=0.0d+00
      nt=2
      u1=1.0d+00
      v2=1.0d+00
      w3=1.0d+00
      u2=0.0d+00
      u3=0.0d+00
      v1=0.0d+00
      v3=0.0d+00
      w1=0.0d+00
      w2=0.0d+00
      go to 1000
  400 continue
      if(index.eq.11) then
      write(iw,9995) grpmin(11),grpmin(1)
      index=1
      go to 210
      end if
      if(index.eq.12) then
      write(iw,9995) grpmin(12),grpmin(2)
      index=2
      go to 250
      end if
      stop
  500 nn=0
      beta=0.5d+00*alph
      cosb=cos(beta)
      sinb=sin(beta)
  510 nn=nn+1
      if(nn.gt.nt) go to 520
c
c     s2n group
c     the plane of symmetry for the improper rotation
c     is the (x,y) plane of the local frame
c
      i=ii+9*(nn-1)
      t(i+1)= cosb*t(i+1-ii)+sinb*t(i+2-ii)
      t(i+2)=-sinb*t(i+1-ii)+cosb*t(i+2-ii)
      t(i+3)=     -t(i+3-ii)
      t(i+4)= cosb*t(i+4-ii)+sinb*t(i+5-ii)
      t(i+5)=-sinb*t(i+4-ii)+cosb*t(i+5-ii)
      t(i+6)=     -t(i+6-ii)
      t(i+7)= cosb*t(i+7-ii)+sinb*t(i+8-ii)
      t(i+8)=-sinb*t(i+7-ii)+cosb*t(i+8-ii)
      t(i+9)=     -t(i+9-ii)
      go to 510
  520 nt=nt+nt
      ii=9*nt
      if(iop1.ne.3) go to 1000
      write(iw,9979)
      n1=n2+1
      n2=n2+naxis
      call print(n1,n2)
      go to 1000
c
c     t group and others containing a subgroup t _
c     local x,y,z are the c2 axes
c
  600 do 610 i=10,36
  610 t(i)=0.0d+00
      t(10)=1.0d+00
      t(23)=1.0d+00
      t(36)=1.0d+00
      t(14)=-1.0d+00
      t(18)=-1.0d+00
      t(19)=-1.0d+00
      t(27)=-1.0d+00
      t(28)=-1.0d+00
      t(32)=-1.0d+00
      do 620 ii=5,12
      i=9*(ii-1)
      j=9*(ii-5)
      t(i+1)=t(j+7)
      t(i+2)=t(j+8)
      t(i+3)=t(j+9)
      t(i+4)=t(j+1)
      t(i+5)=t(j+2)
      t(i+6)=t(j+3)
      t(i+7)=t(j+4)
      t(i+8)=t(j+5)
  620 t(i+9)=t(j+6)
      nt=12
      if(index.eq.13) go to 1000
      if(index.eq.14) go to 650
      if(index.eq.15) go to 680
      go to 670
c
c     th group
c     expand group by taking direct product with ci
c
  650 i=9*nt
      do 660 j=1,i
  660 t(j+i)=-t(j)
      nt=nt+nt
      go to 1000
c
c     octahedral group is direct product of t and a c4 rotation
c     about z axis
c
  670 sign=-1.0d+00
      go to 685
c
c     td group is direct product of t and a reflection in a
c     plane ( equation of the plane   x=y  )
c
  680 sign=1.0d+00
  685 do 690 ii=13,24
      i=9*(ii-1)
      j=9*(ii-13)
      t(i+1)=t(j+4)*sign
      t(i+2)=t(j+5)*sign
      t(i+3)=t(j+6)*sign
      t(i+4)=t(j+1)
      t(i+5)=t(j+2)
      t(i+6)=t(j+3)
      t(i+7)=t(j+7)
      t(i+8)=t(j+8)
  690 t(i+9)=t(j+9)
      nt=24
      if(index.ne.17) go to 1000
c
c     oh group is direct product of o and ci
c
      i=9*nt
      do 700 j=1,i
  700 t(j+i)=-t(j)
      nt=48
 1000 continue
c
c     find the inverse transformations
c
      do 1002 itr=1,nt
      nn=9*(itr-1)
      do 1001 it=1,nt
      ii=9*(it-1)
      test= t(nn+1)*t(ii+1)+t(nn+2)*t(ii+4)+t(nn+3)*t(ii+7)
     1     +t(nn+4)*t(ii+2)+t(nn+5)*t(ii+5)+t(nn+6)*t(ii+8)
     1     +t(nn+7)*t(ii+3)+t(nn+8)*t(ii+6)+t(nn+9)*t(ii+9)
     1     -3.0d+00
      if(abs(test).gt.tol) go to 1001
      invt(itr)=it
      go to 1002
 1001 continue
 1002 continue
      return
      end
