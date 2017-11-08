      subroutine qmat
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      common/infoa/nat,ich,mul,numscf,nx,ne,na,nb,zan(dc),c(3,dc)
      common/iofile/ir,iw,ip,is
      common/eigen2/s(doas),v(doasq),e(doa),ia(doa),in(doa)
      do 100 i=1,numscf
      in(i)=(i-1)*numscf
  100 ia(i)=i*(i-1)/2
c      call hdiag(s,v,e,ia,in,numscf)
      call diagonaliser(numscf,s,e,v,numscf)       
      dum=e(1)
      j=0
      do 200 i=1,numscf
      if(e(i).lt.1.0d-04) j=j+1
  200 e(i)=1.0d+00/sqrt(e(i))
      if(j.ne.0) write(iw,1001) dum,j
 1001 format(' ----- warning -----',/,
     1 ' the smallest eigenvalue of the overlap matrix is = ',
     1 f15.8,/,' there is(are) ',i5,' eigenvalue(s) less than 1.0d-04')
      do 300 i=1,numscf
      do 300 j=1,i
      dum=0.0d+00
      ij=ia(i)+j
      do 250 k=1,numscf
  250 dum=dum+e(k)*v(in(k)+i)*v(in(k)+j)
  300 s(ij)=dum
      call wrtda(s,nx,4)
c
      write(iw,1000)
 1000 format(' ...... end of computation of the fock',
     1 ' transformation matrix q ......')
      return
      end
