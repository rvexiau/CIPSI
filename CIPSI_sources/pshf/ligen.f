      subroutine ligen (a,vec,eig,ia,in,n)
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
c        *****  a givens housholder matrix diagonalization   *****
c        *****  routine same as eigen but works with a       *****
c        *****  linear array.                                *****
      dimension a(*),vec(*),eig(*),ia(*),in(*)
      dimension w(doa),gamma(doa),beta(doa),betasq(doa)
      dimension p(doa),q(doa),iposv(doa),ivpos(doa),iord(doa)
      equivalence (iposv(1),gamma(1)),(ivpos(1),beta(1)),
     1 (iord(1),betasq(1)),(p(1),beta(1)),(q(1),beta(1))
      rhosq=1.d-18
      thr=9.35d-24
      if(n.eq.0) go to 560
      n1=n-1
      n2=n-2
      gamma(1)=a(1)
      if(n2) 280,270,120
  120 do 260 nr=1,n2
      ik=ia(nr+1)+nr
      b=a(ik)
      s=0.0d+00
      do 130 i=nr,n2
      ij=ia(i+2)+nr
      aij=abs(a(ij))
c     if(aij) too small  skip multiplication
c     this is linked to a specific failure of b6700 forttran
      if(aij.lt.thr) go to 130
      s=s+aij*aij
  130 continue
c        *****  prepare for possible bypass of transformation ****
      a(ik)=0.0d+00
      if(s.le.0.0d+00) go to 250
      s=s+b*b
      sgn=+1.0d+00
      if(b.ge.0.0d+00) go to 160
      sgn=-1.0d+00
  160 sqrts= sqrt(s)
      d=sgn/(sqrts+sqrts)
      temp= sqrt(.5d+00+b*d)
      w(nr)=temp
      a(ik)=temp
      d=d/temp
      b=-sgn*sqrts
c        *****  -d- is factor of proportionality. now       *****
c        *****  compute and save -w- vector. extra singly   *****
c        *****  subscripted -w- vector for speed.           *****
      do 170 i=nr,n2
      ij=ia(i+2)+nr
      temp=d*a(ij)
      w(i+1)=temp
  170 a(ij)=temp
c        *****  premultiply vector -w- by matrix -a- to     *****
c        *****  obtain -p- vector. simultaneously accumulate ****
c        *****  dot product -wp- -- scalr -k-.              *****
      wtaw=0.0d+00
      do 220 i=nr,n1
      sum=0.0d+00
      ii=ia(i+1)
      do 180 j=nr,i
      ij=ii+j+1
  180 sum=sum+a(ij)*w(j)
      i1=i+1
      if(n1.lt.i1) go to 210
      do 200 j=i1,n1
      ij=ia(j+1)+i+1
  200 sum=sum+a(ij)*w(j)
  210 p(i)=sum
  220 wtaw=wtaw+sum*w(i)
      do 230 i=nr,n1
  230 q(i)=p(i)-wtaw*w(i)
c        *****  now form -pap- matrix, required part        *****
      do 240 j=nr,n1
      qj=q(j)
      wj=w(j)
      jj=j+1
      do 240 i=j,n1
      ij=ia(i+1)+jj
  240 a(ij)=a(ij)-2.0d+00*(w(i)*qj+wj*q(i))
  250 beta(nr)=b
      betasq(nr)=b*b
      il=ik+1
  260 gamma(nr+1)=a(il)
  270 ij=ia(n)+n-1
      b=a(ij)
      beta(n-1)=b
      betasq(n-1)=b*b
      ij=ij+1
      gamma(n)=a(ij)
  280 betasq(n)=0.0d+00
c        *****  adjoin an identyty matrix to be post-       *****
c        *****  multiplied by rotations                     *****
      nn=n*n
      do 299 i=1,nn
  299 vec(i)=0.0d+00
      do 300 i=1,n
      ij=i+(i-1)*n
  300 vec(ij)=1.0d+00
      m=n
      sum=0.0d+00
      npas=1
      go to 400
  310 sum=sum+shift
      cosa=1.0d+00
      g=gamma(1)-shift
      pp=g
      ppbs=pp*pp+betasq(1)
      ppbr= sqrt(ppbs)
      do 370 j=1,m
      cosap=cosa
      if(ppbs.ne.0.0d+00) go to 320
      sina=0.0d+00
      sina2=0.0d+00
      cosa=1.0d+00
      go to 350
  320 sina=beta(j)/ppbr
      sina2=betasq(j)/ppbs
      cosa=pp/ppbr
c        *****  postmultiply identity by -p- transpose .    *****
      nt=j+npas
      if(nt.lt.n) go to 330
      nt=n
  330 ij=(j-1)*n
      ik=ij+n
      do 340 i=1,nt
      ii=i+ij
      il=i+ik
      temp=cosa*vec(ii)+sina*vec(il)
      vec(il)=-sina*vec(ii)+cosa*vec(il)
  340 vec(ii)=temp
  350 dia=gamma(j+1)-shift
      u=sina2*(g+dia)
      gamma(j)=g+u
      g=dia-u
      pp=dia*cosa-sina*cosap*beta(j)
      if(j.ne.m) go to 360
      beta(j)=sina*pp
      betasq(j)=sina2*pp*pp
      go to 380
  360 ppbs=pp*pp+betasq(j+1)
      ppbr= sqrt(ppbs)
      beta(j)=sina*ppbr
  370 betasq(j)=sina2*ppbs
  380 gamma(m+1)=g
c        *****  test for convergence of last diagonal element ****
      npas=npas+1
      if(betasq(m).gt.rhosq) go to 410
  390 eig(m+1)=gamma(m+1)+sum
  400 beta(m)=0.0d+00
      betasq(m)=0.0d+00
      m=m-1
      if(m.eq.0) go to 430
      if(betasq(m).le.rhosq) go to 390
c        *****  take root of cormer 2 by 2 nearest to       *****
c        *****  lower diagonal in value as estimate of      *****
c        *****  eigenvalue to use for shift                 *****
  410 a2=gamma(m+1)
      r2=0.5d+00*a2
      r1=0.5d+00*gamma(m)
      r12=r1+r2
      dif=r1-r2
      temp= sqrt(dif*dif+betasq(m))
      r1=r12+temp
      r2=r12-temp
      dif= abs(a2-r1)- abs(a2-r2)
      if(dif.lt.0.0d+00) go to 420
      shift=r2
      go to 310
  420 shift=r1
      go to 310
  430 eig(1)=gamma(1)+sum
      do 440 j=1,n
      iposv(j)=j
      ivpos(j)=j
  440  iord(j)=j
      m=n
      go to 470
  450 do 460 j=1,m
      if(eig(j).le.eig(j+1)) go to 460
      temp=eig(j)
      eig(j)=eig(j+1)
      eig(j+1)=temp
      itemp=iord(j)
      iord(j)=iord(j+1)
      iord(j+1)=itemp
  460 continue
  470 m=m-1
      if(m.ne.0) go to 450
      if(n1.eq.0) go to 500
      do 490 l=1,n1
      nv=iord(l)
      np=iposv(nv)
      if(np.eq.l) go to 490
      lv=ivpos(l)
      ivpos(np)=lv
      iposv(lv)=np
      il=(l-1)*n
      ik=(np-1)*n
      do 480 i=1,n
      ii=i+il
      jj=i+ik
      temp=vec(ii)
      vec(ii)=vec(jj)
  480 vec(jj)=temp
  490 continue
c        *****  back transform the vectors of the triple    *****
c        *****  diagonal matrix.                            *****
  500 do 550 nrr=1,n
      k=n1
  510 k=k-1
      if(k.le.0) go to 550
      sum=0.0d+00
      do 520 i=k,n1
      ij=ia(i+1)+k
      ik=i+1+(nrr-1)*n
  520 sum=sum+vec(ik)*a(ij)
      sum=sum+sum
      do 530 i=k,n1
      ij=ia(i+1)+k
      ik=i+1+(nrr-1)*n
  530 vec(ik)=vec(ik)-sum*a(ij)
      go to 510
  550 continue
  560 continue
      return
      end
