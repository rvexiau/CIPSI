subroutine diagonaliser( n, hmat, root, vect, nvec )

   integer, intent( in ) :: n
   real( kind=8 ), intent( in ), dimension(*) :: hmat
   real( kind=8 ), intent( inout ), dimension(*) :: root
   real( kind=8 ), intent( inout ), dimension(*) :: vect
    
   integer :: nb_valp, info,ij,i,j
   integer, dimension(:),allocatable :: ifail,iwork
   real( kind=8 ), dimension(:,:),allocatable :: tab_vect    
   real( kind=8 ), dimension(:),allocatable :: tab_root,work    
   real( kind=8 ), dimension(:),allocatable :: mat_lo
    
   allocate(tab_vect(n,n),tab_root(n))
   allocate(work(8*n),iwork(5*n),ifail(n))
   allocate(mat_lo(n*(n+1)/2))
   mat_lo( 1: n*(n+1)/2 ) = hmat(  1: n*(n+1)/2 )

    !    DSPEVX (JOBZ, RRANGE, UPLO, N, AP, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, IWORK, IFAIL, INFO)
   call DSPEVX ( 'V', 'I', 'L', n, mat_lo, 0._8, 0._8, 1, nvec,  1d-15, nb_valp, tab_root,&
                &tab_vect, n, work, iwork, ifail, info)
                
   ij=1
   do i=1,nvec
     do j=1,n
       vect(ij) = tab_vect(j,i)
       ij=ij+1
     enddo
   enddo  
   root(1:nvec)=tab_root(1:nvec)
   deallocate(tab_vect,tab_root,work,iwork,ifail,mat_lo)
    
    
   return      
    
end subroutine diagonaliser