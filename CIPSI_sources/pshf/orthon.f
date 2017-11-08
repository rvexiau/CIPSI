      subroutine orthon (s,v,e,g,ia,in,ngel,q,numgel,nav,nvec,scfp)
      implicit real*8(a-h,o-z)
      logical pop,heff,scfp
c     orthonormalization of trial vectors
c     uses s matrix on tape (is)
c     uses temporarily tape (ih)
      common/restar/timlim,irest,nrec,intloc,ist,jst,kst,lst
c
      common/iofile/ ir,iw,ip,is,iq,ih
      common/infoa/ nat,ich,mul,numscf,nx,ne,na,nb
      dimension s(*),v(*),e(*),g(*),ia(*),in(*),ngel(*),q(*)
      common/output/ nprint
c
c
c
c
      nxa=nvec*numscf
c     read overlap matrix
c
c
      call reada(s,nx,2)
      if(nprint.eq.4) go to 150
      write(iw,1008)
      max=10
      imax=0
  110 imin=imax+1
      imax=imax+max
      if(imax.gt.nvec) imax=nvec
      write(iw,1003)
      write(iw,1004) (i,i=imin,imax)
      write(iw,1003)
      do 130 j=1,numscf
  130 write(iw,1005) j,(v(in(i)+j),i=imin,imax)
      if(imax.lt.nvec) go to 110
  150 continue
      if(scfp.and.nvec.eq.numscf) call ortit(s,v,e,g,ia,in,ngel,numgel,
     1 nav)
      if(scfp.and.nvec.eq.numscf) go to 441
c     normalization
      do 30 i=1,nvec
      ni=in(i)
      dum=0.0d+00
      do 10 k=1,numscf
      do 10 l=1,k
      nk=ni+k
      kl=ia(k)+l
      x=v(nk)*v(ni+l)*s(kl)
      if(k.ne.l) x=x*2.0d+00
   10 dum=dum+x
      if(dum.gt.1.d-4) go to 12
      write(iw,1006) i,dum
      stop
   12 continue
      dum=1.0d+00/sqrt(dum)
      do 20 k=1,numscf
   20 v(ni+k)=v(ni+k)*dum
   30 continue
c     save temporarily trial vectors on tape (ih)
      call wrtda(v,nxa,nav)
c     calculate v*s*v
      do 90 i=1,nvec
      ni=in(i)
      nn=ia(i)
      do 60 k=1,numscf
      dum=0.0d+00
      do 50 l=1,numscf
      if(k.lt.l) go to 40
      kl=ia(k)+l
      go to 50
   40 kl=ia(l)+k
   50 dum=dum+s(kl)*v(ni+l)
   60 e(k)=dum
      do 80 j=i,nvec
      m=in(j)
      dum=0.0d+00
      do 70 k=1,numscf
   70 dum=dum+v(m+k)*e(k)
      if(i.eq.j) go to 80
      if(abs(dum).lt.0.9) go to 80
      write(iw,1007) i,j,dum
      stop
   80 v(ni+j)=dum
   90 continue
      naa=(nvec*(nvec+1))/2
c     transfer of v into s
      do 95 i=1,nvec
      do 95 j=1,i
      ni=ia(i)+j
      nj=in(j)+i
   95 s(ni)=v(nj)
c
c     diagonalize s and calculate s-1/2
c      call ligen (s,v,e,ia,in,nvec)
      call diagonaliser(nvec,s,e,v,nvec)
      nasq=nvec*nvec
      dum=e(1)
      j=0
      do 200 i=1,nvec
      if(e(i).lt.1.0d-04) j=j+1
      if(e(i).lt.dum) dum=e(i)
  200 e(i)=1.0d+00/sqrt(e(i))
      if(j.ne.0) write(iw,1001) dum,j
      do 300 i=1,nvec
       do 300 j=1,i
      dum=0.0d+00
      ij=ia(i)+j
      do 250 k=1,nvec
      nk=(k-1)*nvec
  250 dum=dum+e(k)*v(nk+i)*v(nk+j)
  300 s(ij)=dum
c
c     transform trial vectors
c
      call reada(v,nxa,nav)
c
      do 440 k=1,numscf
      do 420 i=1,nvec
      dum=0.0d+00
      do 410 j=1,nvec
      if(i.lt.j) go to 400
      ij=ia(i)+j
      go to 410
  400 ij=ia(j)+i
  410 dum=dum+v(in(j)+k)*s(ij)
  420 e(i)=dum
      do 430 i=1,nvec
      if(abs(e(i)).lt.1.d-8) e(i)=0.0d+00
  430 v(in(i)+k)=e(i)
  440 continue
c
  441 continue
c     verification c
c
c
c
      call reada(s,nx,2)
      do 3050 i=1,nvec
      do 3030 k=1,numscf
      dum=0.0d+00
      do 3020 l=1,numscf
      if(k.lt.l) go to 3010
      kl=ia(k)+l
      go to 3020
 3010 kl=ia(l)+k
 3020 dum=dum+s(kl)*v(in(i)+l)
 3030 e(k)=dum
      do 3050 j=1,i
      dum=0.0d+00
      do 3040 k=1,numscf
 3040 dum=dum+v(in(j)+k)*e(k)
      if(i.eq.j) dum=dum-1.0d+00
      if(abs(dum).lt.1.0d-08) go to 3050
      write(iw,1030) i,j,dum
 1030 format(' erreur orthon',2i3,d20.8)
 3050 continue
      if(nprint.eq.4) go to 510
      write(iw,1002)
      max=10
      imax=0
  480 imin=imax+1
      imax=imax+max
      if(imax.gt.nvec) imax=nvec
      write(iw,1003)
      write(iw,1004)(i,i=imin,imax)
      write(iw,1003)
      do 500 j=1,numscf
  500 write(iw,1005) j,(v(in(i)+j),i=imin,imax)
      if(imax.lt.nvec) go to 480
  510 continue
c
c
c     form density matrix
      do 530 i=1,numscf
      do 530 j=1,i
      ij=ia(i)+j
      dum=0.0d+00
      do 520 k=1,numscf
      if(q(k).eq.0.d0) go to 520
      ni=in(k)+i
      nj=in(k)+j
      dum=dum+v(ni)*v(nj)*q(k)
  520 continue
  530 s(ij)=dum
      return
 1001 format(' ----- warning -----',/,
     1 ' the smallest eigenvalue of the overlap matrix is =',
     1 f15.8,/,' there is (are) ',i5,' eigenvalue(s) less than 1.0d-04')
 1002 format(/10x,28(1h-)/10x,'orthonormalized trial vectors'/
     1 10x,28(1h-))
 1003 format(/)
 1004 format(5x,10(4x,i3,5x))
 1005 format(i5,10f12.7)
 1006 format(/' error ... error'/,' norm of trial vector ',i3,' is too s
     1mall',d20.8,/' program stops')
 1007 format(/' error ... error',/,' overlap between trial vectors ',
     1 2i3,' is too large',d20.8,/,' program stops')
 1008 format(/10x,13(1h-)/10x,'trial vectors'/10x,13(1h-))
      end
