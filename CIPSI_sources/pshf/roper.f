      subroutine roper(r,p1,p2,p3,f1,f2,f3,pv,ia)
      implicit real*8(a-h,o-z)
      include 'pshf.prm'
      common/infoa/nat,ich,mul,n,nx
      common/open/nc,noc(3),oc(3),alpha(18),rland
      dimension r(*),pv(*),p1(*),p2(*),p3(*),f1(*),f2(*),f3(*)
      dimension ia(*),e(doa),t(doa),e2(doa),e3(doa)
      dimension t2(doa),t3(doa)
c *** r operator
c *** r0=(p+pv)+*f*(p+pv) added over all shells
c *** p=c*c+ , c: mos of the corresponding shell
      if(nc.lt.2) go to 500
      do 1 i=1,nx
1     r(i)=0.d0
      nv1=16
      nv2=5
      nv3=16+nc
      call reada(pv,nx,nv3)
      do 100 nk=1,nc
      call reada(p3,nx,nv1)
      call reada(f3,nx,nv2)
      do 10 i=1,n
      do 5 l=1,n
      e(l)=0.d0
      do 5 k=1,n
      if(i.ge.k) then
         ik=ia(i)+k
         else
         ik=ia(k)+i
      endif
      if(k.ge.l) then
         kl=ia(k)+l
         else
         kl=ia(l)+k
      endif
      pik=p3(ik)+pv(ik)
      e(l)=e(l)+pik*f3(kl)
5     continue
      do 10 j=1,i
      ij=ia(i)+j
      do 10 l=1,n
      if(j.ge.l) then
         jl=ia(j)+l
         else
         jl=ia(l)+j
      endif
      pjl=p3(jl)+pv(jl)
      r(ij)=r(ij)+e(l)*pjl
10    continue
      go to (50,60,100),nk
50    do 55 i=1,nx
      p1(i)=p3(i)
55    f1(i)=f3(i)
      go to 70
60    do 65 i=1,nx
      f2(i)=f3(i)
65    p2(i)=p3(i)
70    nv1=nv1+1
      nv2=10+(nk-1)*3
100   continue
      if(abs(rland).lt.0.d-04) go to 201
      do 110 i=1,nx
      pv(i)=f1(i)-f3(i)
      f1(i)=f1(i)-f2(i)
      f3(i)=-f3(i)+f2(i)
110   f2(i)=pv(i)
      do 200 i=1,n
      do 150 l=1,n
      e(l)=0.d0
      e2(l)=0.d0
      e3(l)=0.d0
      t(l)=0.d0
      t2(l)=0.d0
      t3(l)=0.d0
      do 150 k=1,n
      if(i.ge.k) then
         ik=ia(i)+k
         else
         ik=ia(k)+i
      endif
      if(k.ge.l)then
         kl=ia(k)+l
         else
         kl=ia(l)+k
      endif
      t(l)=t(l)+p1(ik)*f1(kl)
      e(l)=e(l)+p2(ik)*f1(kl)
      if(nc.eq.2) go to 150
      t2(l)=t2(l)+p1(ik)*f2(kl)
      t3(l)=t3(l)+p2(ik)*f3(kl)
      e2(l)=e2(l)+p3(ik)*f2(kl)
      e3(l)=e3(l)+p3(ik)*f3(kl)
150   continue
      do 200 j=1,i
      ij=ia(i)+j
      do 180 l=1,n
      if(j.ge.l) then
         jl=ia(j)+l
         else
         jl=ia(l)+j
      endif
      r(ij)=r(ij)+rland*(t(l)*p2(jl)+e(l)*p1(jl))
      if(nc.lt.3) go to 180
180   r(ij)=r(ij)+rland*(t2(l)*p3(jl)+e2(l)*p1(jl)+t3(l)*p3(jl)+
     1            e3(l)*p2(jl))
200   continue
c
c *** s*r0*s
c
201   do 210 i=1,nx
210   f1(i)=0.d0
      call reada(f2,nx,2)
      do 300 i=1,n
      do 250 l=1,n
      e(l)=0.d0
      do 250 k=1,n
      if(i.ge.k) then
         ik=ia(i)+k
         else
         ik=ia(k)+i
      endif
      if(k.ge.l) then
         kl=ia(k)+l
         else
         kl=ia(l)+k
      endif
      e(l)=e(l)+f2(ik)*r(kl)
250   continue
      do 300 j=1,i
      ij=ia(i)+j
      do 300 l=1,n
      if(j.ge.l) then
         jl=ia(j)+l
         else
         jl=ia(l)+j
      endif
      f1(ij)=f1(ij)+e(l)*f2(jl)
300   continue
      do 400 i=1,nx
      r(i)=f1(i)
400   pv(i)=p1(i)+p2(i)
      call wrtda(r,nx,5)
      if(nc.lt.3) go to 500
      do 450 i=1,nx
450   pv(i)=pv(i)+p3(i)
500   return
      end
