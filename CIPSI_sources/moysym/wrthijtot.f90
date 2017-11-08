subroutine wrthijtot(nbuf,itab,jtab,vbuf,nhij,ncf,ii,jj,hij,numtab)
  implicit none
  integer,parameter :: lsize=8192
  integer,parameter :: smax=6  

! no Spin sorting. Everything is written on the singlet file/arrays
  integer,dimension(ncf,smax),intent(inout) :: numtab  
  integer,intent(inout) :: nbuf,itab(lsize,smax),jtab(lsize,smax),nhij
  real*8,intent(inout) :: vbuf(lsize,smax)
  integer,intent(in) :: ii,jj,ncf
  real*8 :: hij

  if(.not.(ii.eq.jj)) then
    if(dabs(hij).lt.1.d-10) return
  endif
  numtab(ii,1)=ii
  nhij=nhij+1          
  nbuf=nbuf+1
  itab(nbuf,1)=ii
  jtab(nbuf,1)=jj
  if(ii.eq.jj) hij=0.5d0*hij
  vbuf(nbuf,1)=hij 

  if(nbuf.eq.lsize) then
    write(71) itab(1:lsize,1),jtab(1:lsize,1),vbuf(1:lsize,1)
    nbuf=0
  endif
  return
     
end subroutine