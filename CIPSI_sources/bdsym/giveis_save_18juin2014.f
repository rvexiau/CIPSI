      subroutine giveis(n,nv,nd,a,b,root,vect,ierr)
Ccc   RG : modified on 05/01/09 : Use a LAPACK routine, hopefully parallelized
      implicit real*8(a-h,o,p,r-z),logical*4(q)
      dimension a(1),vect(1),root(1)
      common/info/ncf,ncper,nca,qnorm,qwa,qwb,mbuf,metat2,netat,niter
     *,qwvec,qhef,ntrsy,qrvec,seuil,qws,qsuis,qsuip,qener,njter,qgiv
     *,qictot,qsui04,qdcple

!      real*8 ala(n,n),wla(n),workla(3*n-1)
      real*8,dimension(8*n) :: work
      integer*4,dimension(5*n) :: iwork 
      integer*4,dimension(n) :: ifail

      if (qgiv) then
         nev=n
         nvec=n
         ijsup=0
         do 5  i=1,n
         do 5  j=i,n
      ijsup=ijsup+1
      ijinf=j*(j-1)/2+i
5     vect(ijsup)=a(ijinf)
      nij=n*(n+1)/2
      do 8 ijsup=1,nij
8     a(ijsup)=vect(ijsup)
cCc      call givens(a,root,vect,n,nev,nvec)
call dspevx('V','I','U',n,a,0.d0,0.d0,1,nev,1.d-10,m,root,vect,nvec,work,iwork,ifail,info)
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
