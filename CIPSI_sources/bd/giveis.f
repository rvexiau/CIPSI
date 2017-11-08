!
!     DB (06/2014) : 
!          Cette subroutine n'est plus utilisée et est remplacée par 
!          la subroutine 'diagonaliser'
!


      subroutine giveis(n,nv,nd,a,b,root,vect,ierr)
Ccc   RG : modified on 05/01/09 : Use a LAPACK routine, hopefully parallelized
      implicit real*8(a-h,o,p,r-z),logical*4(q)
      dimension a(1),vect(1),root(1)
      common/info/ncf,ncper,nca,qnorm,qwa,qwb,mbuf,metat2,netat,niter
     *,qwvec,qhef,ntrsy,qrvec,seuil,qws,qsuis,qsuip,qener,njter,qgiv
     *,qictot,qsui04,qdcple

!      real*8 ala(n,n),wla(n),workla(3*n-1)
! Pour DSPEVX :
!       real*8,dimension(8*n) :: work
!       integer*4,dimension(5*n) :: iwork 
!       integer*4,dimension(n) :: ifail
!       integer*4 :: nb_valp, info
!       real*8, dimension( n, n ) :: tab_vectp
!       real*8, dimension( n ) :: tab_valp
!       real*8, dimension( n*(n+1)/2 ) :: mat_lo
!       real*8, dimension( n, n ) :: matrice
    
!      mat_lo( : ) = a( 1 : n * ( n + 1 ) / 2 )
   
c    Pour séparer une ligne de commande en 2, il faut commencer la 2e ligne après
c    le 5e caractère, merci fortran 77 :)
!       open( file = 'prout1', unit = 122, status = 'unknown',
!      +position = 'append', action = 'write' )
!       write( 122, * ) 'n = ', n
!       
!       Ça, ça sert juste à transformer la matrice triangulaire inférieure sous forme compacte
!       mat_lo en une vraie matrice pour l'imprimer :
!       do i = 1, n
!         do j = 1, i
!         
!             matrice( i, j ) = mat_lo( i + (j-1)*(2*n-j)/2 )
!             matrice( j, i ) = matrice( i, j )
!         end do
!       end do  
!       write( 122, '(<n>(e9.2,2X))' ) matrice( :, : )
!       !write( 122, '(25(e9.2,2X))' ) mat_lo( : )
!       close( 122 )
!     
!       !AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n
!       
!      
!       
! call DSPEVX( 'V', 'I', 'L', n, mat_lo, 0.d0, 0.d0, 1, n, 1.d-9, nb_valp, tab_valp, tab_vectp, n, work, iwork, ifail, info)      
! 
!       open( file = 'prout2', unit = 122, status = 'unknown',
!      +position = 'append', action = 'write' )
!
!       write( 122, * ) 'n = ', n, 'info = ', info, 'nb_valp =', nb_valp
!       write( 122, '(25(f6.3,2X))' ) tab_valp( : ) 
!       !write( 122, '(<n>(e9.2,2X))' ) tab_vectp( :, : )
!       close( 122 )
! 
!       vect( 1 : n*n ) = reshape( tab_vectp, shape( vect( 1 : n*n ) ) ) 
!       root( 1 : n ) = tab_valp
!       
!       write(*,*) 'DSPEVX : n/nev/nvec/INFO ',n,nev,nvec,info
       
      if (qgiv) then
         nev = n
         nvec = n 

c       Cela sert à ranger la matrice sous forme triangulaire supérieure, cela est nécessaire pour GIVENS,
c       car au départ elle est rangée comme triangulaire inférieure :
         ijsup = 0
         do 5  i = 1, n
            do 5  j = i, n
                ijsup = ijsup + 1
                ijinf = j * ( j - 1 ) / 2 + i
5               vect( ijsup ) = a( ijinf )

         nij = n * ( n + 1 ) / 2
                
      do 8 ijsup = 1, nij
8         a( ijsup ) = vect( ijsup )

      call givens(a,root,vect,n,nev,nvec)

      open( file = 'prout2', unit = 122, status = 'unknown',
     +position = 'append', action = 'write' )

      write( 122, * ) 'n = ', n, 'size_root =', size( root )
      write( 122, '(20(f9.6,2X))' ) root( 1:n ) 
      !write( 122, '(<n>(e9.2,2X))' ) tab_vectp( :, : )
      close( 122 )

         write(*,*) 'DSPEVX : n/nev/nvec/INFO ',n,nev,nvec,info
         
      else
         call jacscf(a,vect,root,n,-1,1.d-12)
!	JD 06/06
!	 ii=0
!	 do j=1,n
!	  do i=1,j
!	   ii=ii+1
!	   ala(i,j)=a(ii)
!	  enddo
 !        enddo
!	 do j=1,n
!	  do i=j+1,n
!	   ala(i,j)=ala(j,i)
!	  enddo
 !        enddo
!	 lwork=(3*n-1)
!	 call DSYEV('V','L',n,ala,n,wla,workla,lwork,infola)
!	 write(6,*)'call DSYEV: info',infola
!	 write(6,*)'LWORK, OPTIMAL LWORK',(3*n-1),lwork
!	 write(6,*)'eigenvalues'
!	 write(6,*)(wla(i),i=1,n)
!         write(6,*)'eigenvectors'
!	 do i=1,n
!	  write(6,*)(ala(i,j),j=1,n)
!         enddo
!	 do i=1,n
!	  root(i)=wla(i)
!	  do j=1,n
!	   vect((i-1)*n+j)=ala(i,j)
!	  enddo
!         enddo
!	JD 06/06
      endif
      return
      end
