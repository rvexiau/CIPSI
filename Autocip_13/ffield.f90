module fitting
  contains
 
  function polyfit(vx, vy, d)
    implicit none
    integer, intent(in)                   :: d
    integer, parameter                    :: dp = selected_real_kind(15, 307)
    real(dp), dimension(d+1)              :: polyfit
    real(dp), dimension(:), intent(in)    :: vx, vy
 
    real(dp), dimension(:,:), allocatable :: X
    real(dp), dimension(:,:), allocatable :: XT
    real(dp), dimension(:,:), allocatable :: XTX
 
    integer :: i, j
 
    integer     :: n, lda, lwork
    integer :: info
    integer, dimension(:), allocatable :: ipiv
    real(dp), dimension(:), allocatable :: work
 
    n = d+1
    lda = n
    lwork = n
 
    allocate(ipiv(n))
    allocate(work(lwork))
    allocate(XT(n, size(vx)))
    allocate(X(size(vx), n))
    allocate(XTX(n, n))
 
    ! prepare the matrix
    do i = 0, d
       do j = 1, size(vx)
          X(j, i+1) = vx(j)**i
       end do
    end do
 
    XT  = transpose(X)
    XTX = matmul(XT, X)
 
    ! calls to LAPACK subs DGETRF and DGETRI
    call DGETRF(n, n, XTX, lda, ipiv, info)
    if ( info /= 0 ) then
       print *, "problem"
       return
    end if
    call DGETRI(n, XTX, lda, ipiv, work, lwork, info)
    if ( info /= 0 ) then
       print *, "problem"
       return
    end if
 
    polyfit = matmul( matmul(XTX, XT), vy)
 
    deallocate(ipiv)
    deallocate(work)
    deallocate(X)
    deallocate(XT)
    deallocate(XTX)
 
  end function
 
end module

subroutine finitefield(field_dir,sym,l_accept,nfield,field,netat)
! resort autocip output for easy finite field computation
! future development : automatically do the interpolation to get the dipole/polarisability
  
  use grid_data
  use io_unit
  use fitting
  implicit none

  integer :: nfield,netat
  real(dp),dimension(nfield) :: field
  logical,dimension(nz+1,nfield) :: l_accept
  character(len=*) :: field_dir,sym
  
  integer :: ifield,i,j,k,lefield_chars
  real(dp) :: temp,dipole,pola
  real(dp),dimension(3) :: fit
  real(dp),dimension(:,:,:),allocatable :: energ
  character*80 :: efield_chars,fmt_out,state_char
  
  allocate(energ(nfield,nz,netat))
  energ=0d0
  
  do ifield=1,nfield
    write(efield_chars,'(E9.2)')field(ifield)
    efield_chars=trim(adjustl(efield_chars))
    lefield_chars=len_trim(efield_chars)

    open(unit=9,file=trim(field_dir)//efield_chars(1:lefield_chars)//'/'//'vccind_'//trim(sym),status='old')

    do i=1,nz
      if(.not.l_accept(i,ifield)) cycle
      read(9,*) (temp,j=1,ncoord),(energ(ifield,i,j),j=1,netat)
    enddo    
    close(9)
  enddo
    
  do k=1,netat
    write(state_char,'(i0)')k
    open(9,file='ffield/'//trim(state_char)//trim(sym),status='unknown')
    write(fmt_out,'(a5,i0,a8)')  '(a15,',nfield,'(f20.8))'
    write(9,fmt_out) '# field value :',(field(j),j=1,nfield)
    write(9,'(a63)') '# coordinate,PDM,polarizability,raw energy for each field value'
  
    write(fmt_out,'(a1,i0,a16,i0,a9)')  '(',ncoord,'(f12.4),2f25.12,',nfield,'(f20.12))'
    do i=1,nz
      if(.not.l_accept(i,1)) cycle
      dipole=0d0
      pola=0d0
      if(nfield.gt.2) then
        fit=polyfit(field,energ(:,i,k),2)
        dipole=-fit(2)
        pola=-fit(3)*2d0
      endif
      write(9,fmt_out) (coord(j,i),j=1,ncoord),dipole,pola,(energ(j,i,k),j=1,nfield)
    enddo
    close(9) 
  enddo
  
  deallocate(energ)

end subroutine

function polyfit(vx, vy, d)
  implicit none
  integer, intent(in)                   :: d
  integer, parameter                    :: dp = selected_real_kind(15, 307)
  real(dp), dimension(d+1)              :: polyfit
  real(dp), dimension(:), intent(in)    :: vx, vy
 
  real(dp), dimension(:,:), allocatable :: X
  real(dp), dimension(:,:), allocatable :: XT
  real(dp), dimension(:,:), allocatable :: XTX
 
  integer :: i, j
 
  integer     :: n, lda, lwork
  integer :: info
  integer, dimension(:), allocatable :: ipiv
  real(dp), dimension(:), allocatable :: work
 
  n = d+1
  lda = n
  lwork = n
 
  allocate(ipiv(n))
  allocate(work(lwork))
  allocate(XT(n, size(vx)))
  allocate(X(size(vx), n))
  allocate(XTX(n, n))
 
  ! prepare the matrix
  do i = 0, d
     do j = 1, size(vx)
        X(j, i+1) = vx(j)**i
     end do
  end do
 
  XT  = transpose(X)
  XTX = matmul(XT, X)
 
  ! calls to LAPACK subs DGETRF and DGETRI
  call DGETRF(n, n, XTX, lda, ipiv, info)
  if ( info /= 0 ) then
     print *, "problem"
     return
  end if
  call DGETRI(n, XTX, lda, ipiv, work, lwork, info)
  if ( info /= 0 ) then
     print *, "problem"
     return
  end if
 
  polyfit = matmul( matmul(XTX, XT), vy)
 
  deallocate(ipiv)
  deallocate(work)
  deallocate(X)
  deallocate(XT)
  deallocate(XTX)
 
end function