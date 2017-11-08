      subroutine contr8(h,hs,num,ktype,kloc,kmin,kmax,nshell,ia,newoa,
     #                  noa,maxom)
c
c***********************************************************************
c         ----- transform to the Ylm basis HONDO8 version -----
c=======================================================================
c      H = matrice a tranformer (input)
c     HS = matrice  tranformee (output)
c    NUM = dimension du tableau IA
c  KTYPE = type des couches  (S=1 P=2 etc...)
c   KLOC = localisation de la premiere OA d'un couche dans la base
c   KMIN = numero de la premiere OA dans une couche
c   KMAX = numero de la derniere OA dans une couche
c NSHELL = nombre de couches
c     IA = tableau: ia(i) = i*(i-1)/2
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
      integer MAXL,WTDIM
      parameter (MAXL=15)
      parameter (WTDIM=MAXL*MAXL)
c
c     ----- arguments -----
c
      integer num,nshell,ktype,kloc,kmin,kmax,ia,newoa,noa,maxom
      double precision h,hs
      dimension ktype(nshell),kloc(nshell),kmin(nshell),kmax(nshell),
     #          h(noa,*),hs(maxom,*),ia(num)
c
c     ----- local variables -----
c
      integer ii,jj,loci,locj,imin,imax,jmin,jmax,lit,ljt,i,j,
     #        inew,jnew,iend,jend,nio,njo
      double precision t
      dimension t(MAXL,MAXL)
c
c     +++++
c     start
c     +++++
c
      call dfgset
c
      inew = 0
      do ii=1,nshell
         lit = ktype(ii)
         imin = kmin(ii)
         imax = kmax(ii)
         loci = kloc(ii)-imin
         iend = imax-imin+1
         nio = iend
         jnew = 0
         do jj=1,nshell
            ljt = ktype(jj)
            jmin = kmin(jj)
            jmax = kmax(jj)
            locj = kloc(jj)-jmin
            jend = jmax-jmin+1
            njo = jend
            if(lit.gt.2.or.ljt.gt.2) then
               do j=jmin,jmax
                  do i=imin,imax
                     t(i-imin+1,j-jmin+1) = h(loci+i,locj+j)
                  enddo
               enddo
               call sphrop(t,MAXL,lit,ljt,nio,njo,iend,jend)
               do i=1,iend
                  do j=1,jend
                     hs(inew+i,jnew+j) = t(i,j)
c                    hs(jnew+j,inew+i) = t(j,i)
                  enddo
               enddo
            else
               do j=jmin,jmax
                  do i=imin,imax
                     hs(inew-imin+1+i,jnew-jmin+1+j) = h(loci+i,locj+j)
c                    hs(jnew-jmin+1+j,inew-imin+1+i) = h(locj+j,loci+i)
                  enddo
               enddo
            endif
            jnew = jnew+jend
         enddo
         inew = inew+iend
      enddo
      newoa=inew
      return
      end
