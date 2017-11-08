      subroutine dmop(v,p,a,b,c,p1,p2,p3,in,na)
      implicit real*8 (a-h,o-z)
      character*8 type
      common/open/nc,noc(3),oc(3),alpha(3,3),beta(3,3),rland
      common/output/nprint
      common/iofile/ir,iw
      common/infoa/nat,ich,mul,n,nx
      dimension type(3)
      dimension a(*),b(*),c(*),p(*),v(*),p1(*),p2(*),p3(*),in(*)
      data type/'closed s','open s 1','open s 2'/
9900  format(10x,17(1h-),/,10x,'a matrix ',a8,/,10x,17(1h-))
9901  format(10x,17(1h-),/,10x,'b matrix ',a8,/,10x,17(1h-))
      n1=1
      n2=noc(1)
      ipas=0
1     ipas=ipas+1
      jk=0
      do 3 k=1,n
      do 3 j=1,k
      jk=jk+1
      t=0.d0
      do 2 i=n1,n2
      ji=in(i)+j
      ki=in(i)+k
2     t=t+v(ji)*v(ki)
      p(jk)=t
3     continue
      nv1=15+ipas
      call wrtda(p,nx,nv1)
      if(ipas.gt.nc) go to 40
      n1=n2+1
      go to (5,10,15),ipas
5     do 6 i=1,nx
      c(i)=p(i)*oc(1)
6     p1(i)=p(i)
      go to 20
10    do 11 i=1,nx
      c(i)=c(i)+oc(2)*p(i)
11    p2(i)=p(i)
      go to 20
15    do 16 i=1,nx
      c(i)=c(i)+oc(3)*p(i)
16    p3(i)=p(i)
20    if(nc.gt.1.and.ipas.eq.nc) go to 30
      n2=n2+noc(ipas+1)
      go to 1
30    n2=n
      go to 1
40    continue
c *** c total density matrix
c *** a and b matrices
      do 100 k=1,nc
      do 60 i=1,nx
      a(i)=0.d0
60    b(i)=0.d0
      do 80 i=1,nx
      a(i)=a(i)+alpha(k,1)*p1(i)
      b(i)=b(i)+beta(k,1)*p1(i)
      if (nc.eq.1) go to 80
      a(i)=a(i)+alpha(k,2)*p2(i)
      b(i)=b(i)+beta(k,2)*p2(i)
      if (nc.eq.2) go to 80
      a(i)=a(i)+alpha(k,3)*p3(i)
      b(i)=b(i)+beta(k,3)*p3(i)
80    continue
      nv2=8+(k-1)*3
      call wrtda(a,nx,nv2)
      nv2=nv2+1
      call wrtda(b,nx,nv2)
      if (nprint.ne.5) go to 100
      write(iw,9900)type(k)
      call fout(a,n)
      write(iw,9901)type(k)
      call fout(b,n)
100   continue
      call wrtda(c,nx,6)
      return
      end
