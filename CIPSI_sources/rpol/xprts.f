      subroutine xprts
c
c***********************************************************************
c=======================================================================
c-----------------------------------------------------------------------
c     written by: J. P. Flament
c                 LDMP
c                 Universite de Lille-1
c                 F-59655 Villeneuve d'Ascq Cedex
c***********************************************************************
c
      implicit none
c
c     ----- parameters -----
c
      integer MXSYM
      double precision ZERO,ONE
      parameter (MXSYM=48,ZERO=0.0D+00,ONE=1.0D+00)
c
c     ----- local variables -----
c
      integer it
      integer i,j,k
      double precision temp(15,15)
c
c     ----- variables in common -----
c
      integer invt,nt,ntmax,ntwd,nosym
      double precision ptr,dtr,ftr,gtr,dcoef,fcoef,gcoef,dcart,fcart,
     #                 gcart
c
      common /symtry/ invt(MXSYM),nt,ntmax,ntwd,nosym
c     common /symspd/ ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720)
      common /symspd/ ptr(3,3,48),dtr(6,6,48),ftr(10,10,48),
     #                gtr(15,15,48)
      common /dfgcof/ dcoef(6,6),fcoef(10,10),gcoef(15,15)
      common /cofdfg/ dcart(6,6),fcart(10,10),gcart(15,15)
c
c     +++++
c     start
c     +++++
c
      call dfgset
c
c     orbitales d   
c
      do it=1,nt
         do j=1,6
            do i=1,6
               temp(i,j)=ZERO
               do k=1,6
                  temp(i,j)=temp(i,j)+dcart(i,k)*dtr(k,j,it)
               enddo
            enddo
         enddo
         do j=1,6
            do i=1,6
               dtr(i,j,it)=ZERO
               do k=1,6
                  dtr(i,j,it)=dtr(i,j,it)+dcoef(k,j)*temp(i,k)
               enddo
            enddo
         enddo
      enddo
c
c     orbitales f   
c
      do it=1,nt
         do i=1,10
            do j=1,10
               temp(i,j)=ZERO
               do k=1,10
                  temp(i,j)=temp(i,j)+fcart(i,k)*ftr(k,j,it)
               enddo
            enddo
         enddo
         do j=1,10
            do i=1,10
               ftr(i,j,it)=ZERO
               do k=1,10
                  ftr(i,j,it)=ftr(i,j,it)+fcoef(k,j)*temp(i,k)
               enddo
            enddo
         enddo
      enddo
c
c     orbitales g   
c
      do it=1,nt
         do i=1,15
            do j=1,15
               temp(i,j)=ZERO
               do k=1,15
                  temp(i,j)=temp(i,j)+gcart(i,k)*gtr(k,j,it)
               enddo
            enddo
         enddo
         do j=1,15
            do i=1,15
               gtr(i,j,it)=ZERO
               do k=1,15
                  gtr(i,j,it)=gtr(i,j,it)+gcoef(k,j)*temp(i,k)
               enddo
            enddo
         enddo
      enddo
      return
      end
