      subroutine ortit (f,v,g,h,ia,in,ngel,numgel,nav)
      implicit real*8(a-h,o-z)
      include 'pshf.prm'
      logical convg
      logical gel
      dimension f(*),v(*),g(*),h(*),ia(*),in(*),ngel(*)
      common/iofile/ ir,iw,ip,is,iq,ih,iv
      common/output/nprint
      common/infoa/nat,ich,mul,numscf,nx,ne,na,nb,zan(dc),c(3,dc)
c
c
      data thresh,delmax/1.d-12,0.2/
 9996 format(5x,'orthonormalization process do''nt converge after iterat
     *ion ',i4/'          s-1/2 process is used')
 9995 format(5x,'orthonormalization converges at iteration',i4)
c
      nxsq=numscf*numscf
c     iterative orthonormalization of the basis
c     overlap in f
c     vectors in v
c
      maxit=50
      niter=0
   20 continue
      call reada(f,nx,2)
      call wrtda(v,nxsq,nav)
      niter=niter+1
      convg=.true.
      do 100 i=1,numscf
      ni=in(i)
c     calculate s*v
      do 40 j=1,numscf
   40 g(j)=0.d0
      do 70 k=1,numscf
      dum=v(ni+k)
      if(dum.eq.0.d0) go to 70
      nk=ia(k)
      do 45 j=1,k
   45 g(j)=g(j)+f(nk+j)*dum
      if(k.eq.numscf) go to 70
      kp=k+1
      do 50 j=kp,numscf
   50 g(j)=g(j)+f(ia(j)+k)*dum
   70 continue
c                +
c     calculate v *(s*v)
c
      do 90 j=i,numscf
      h(j)=0.d0
      nj=in(j)
      do 90 k=1,numscf
      dum=v(nj+k)
      if(dum.eq.0.d0) go to 80
      h(j)=h(j)+dum*g(k)
   80 continue
   90 continue
      do 95 j=i,numscf
   95 v(ni+j)=h(j)
  100 continue
      do 105 i=1,numscf
      do 105 j=i,numscf
      if(i.eq.j) g(j)=1.d0/sqrt(v(in(i)+j))
  105 f(ia(j)+i)=v(in(i)+j)
c      write(6,*) 'vsv'
c      call fout(f,numscf)
c
c
      call reada(v,nxsq,nav)
c     normalization at first iteration
      if(niter.gt.1) go to 120
      do 110 i=1,numscf
      ni=in(i)
      do 110 j=1,numscf
  110 v(ni+j)=v(ni+j)*g(i)
      ij=0
      do 115 i=1,numscf
      do 115 j=1,i
      ij=ij+1
  115 f(ij)=f(ij)*g(i)*g(j)
  120 continue
c      (n+1)             (n)    (n)   (n)
c     v     =(3/2*i-1/2*v   *s*v   )*v
c
      do 170 i=1,numscf
      gel=.false.
      if(numgel.eq.0) go to 121
      do 122 j=1,numgel
      if(ngel(j).eq.i) gel=.true.
  122 continue
  121 continue
      ni=in(i)
      do 125 k=1,numscf
  125 g(k)=0.d0
      do 140 j=1,i
      fij=f(ia(i)+j)
      nj=in(j)
      fij=-0.5d0*fij
      if(i.eq.j) fij=1.5d0+fij
      do 130 k=1,numscf
  130 g(k)=g(k)+fij*v(nj+k)
  140 continue
      if(i.eq.numscf) go to 162
      ip=i+1
      do 160 j=ip,numscf
      fij=f(ia(j)+i)
      if(fij.eq.0.d0) go to 160
      nj=in(j)
      fij=-0.5d0*fij
      do 150 k=1,numscf
  150 g(k)=g(k)+fij*v(nj+k)
  160 continue
  162 continue
      do 165 k=1,numscf
      dum=v(ni+k)
      if(gel) g(k)=dum
      if(abs(dum-g(k)).gt.thresh) convg=.false.
  165 v(ni+k)=g(k)
  170 continue
c      write(6,*) 'v,i',niter
c      call vout(v,in,numscf,numscf)
      if(convg) go to 180
      if(niter.lt.maxit) go to 20
      write(iw,9996) niter
      stop
  180 write(iw,9995) niter
      return
      end
