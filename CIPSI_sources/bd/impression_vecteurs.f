      subroutine impression_vecteurs(ncf,netat,cv)
      implicit real*8 (a-h,o-p,r-z),logical*4(q)
      character*1 jspin
!      character(len=20) fname
 
      integer*4 nt,np

      common/tprint/tvec
      dimension nt(8),np(8),jspin(16),jdet(16)
      dimension cv(1)
      rewind 61
!      inquire(61,name=fname)
!      write(6,*)'ATTEEEEENTION: ncf=',ncf
!      write(6,*)'file read is ',fname
      do i=1,ncf
      qpr=.false.
      do j=1,netat
        if(dabs(cv((j-1)*ncf+i)).gt.tvec) qpr=.true.
      end do
      read(61)ne1
!     write(6,*)'ne1=',ne1
      if(ne1.ne.0)then
        read(61)(nt(k),k=1,ne1),(np(k),k=1,ne1)
!          write(6,*)(nt(k),k=1,ne1),(np(k),k=1,ne1)
      endif
      if(qpr)then
      do k=1,ne1
        jspin(k)='+'
        jspin(k+ne1)='+'
        jdet(k)=nt(k)
        jdet(k+ne1)=np(k)
        if(nt(k).gt.norb)then
          jspin(k)='-'
          jdet(k)=jdet(k)-norb
          end if
        if(np(k).gt.norb)then
          jspin(k+ne1)='-'
          jdet(k+ne1)=jdet(k+ne1)-norb
        end if
        end do
      write(6,9998)i,(jdet(k),jspin(k),jdet(k+ne1),jspin(k+ne1),k=1,ne1)
 9998 format(1x,i5,3x,16(i3,a1,1x))
      end if
      end do
      return
      end
