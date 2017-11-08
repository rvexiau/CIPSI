      subroutine efield

      implicit double precision (a-h,o-z)

      include 'pshf.prm'
      real*8 rp(3),ra(3),rb(3)
      real*8 tab(4,10,10)
      integer lit,ljt
      real*8 gama,gamb
      real*8 rcut(dc)
      logical ysoft

      double precision i1sr4
      integer ent,sor,pl,pl1
      parameter (itypmax=3)
      logical*1 iecrit
      common/ecrit/iecrit
c     common/output/nprint,itol,icut
      common/rayon/rcut,ysoft

      common/infoa/nat,nnp,c(3,dc)
      common/infpot/iatno(dc),ipseud(dc)
      include 'NSHEL.cmm'
      common/project/pl
      common/lecrit/ent,sor
      dimension i1sr4(doas),
     1          elx(doas),
     1          ely(doas),
     1          elz(doas)
c
c
c
c
      write(*,*) 'RPOL nnp = ',nnp
      do 10000 ic=1,nat
      iato=iatno(ic)
      if(ipseud(ic).eq.0) go to 10000
      write(sor,*) 'ic=',ic,'rcut=',rcut(ic)
      do 5 i=1,nnp
      i1sr4(i)=0.0d0
      elx(i)=0.d0
      ely(i)=0.d0
      elz(i)=0.d0
    5 continue
      rp(1)=c(1,ic)
      rp(2)=c(2,ic)
      rp(3)=c(3,ic)
c
c     ----- ishell
c
      do 9000 ii=1,nshell
      i=katom(ii)
      ra(1)=c(1,i)
      ra(2)=c(2,i)
      ra(3)=c(3,i)
      gama = ex(kstart(ii))

      lit=ktype(ii)-1
      if((lit+1).gt.4) then
      write(sor,*) 'cas imprevu : lit+1 > 4 '
                       stop
      endif
      kloci=kloc(ii)
c
c     ----- jshell
c
      do 8000 jj=1,ii
      iecrit=ii.eq.11 .and. jj.eq.11 .and. pl.eq.2
c     if(iecrit) write(sor,*) 'ii=',ii,'jj=',jj
      j=katom(jj)
      rb(1)=c(1,j)
      rb(2)=c(2,j)
      rb(3)=c(3,j)
      gamb = ex(kstart(jj))

      ljt=ktype(jj)-1
      if((ljt+1).gt.4) then
      write(sor,*) 'cas imprevu : ljt+1 > 4 '
                       stop
      endif
      klocj=kloc(jj)
c
c

C now calculate 1/r^4 and r/r^3 arrays       
      call xyztrans(tab,rp,ra,rb,lit,ljt,gama,gamb,rcut(ic),ysoft)
 
C store result in tridiagonal matrix       
      do it = 1,2*lit+1
        do jt = 1,2*ljt+1
          ip = kloci+it-1
          jp = klocj+jt-1
          if (ip.ge.jp) then
C only store elements on and below diagonal 
            nij = ip*(ip-1)/2+jp
            i1sr4(nij) = i1sr4(nij) + tab(1,it,jt) 
            elx(nij)   = elx(nij)   + tab(2,it,jt) 
            ely(nij)   = ely(nij)   + tab(3,it,jt) 
            elz(nij)   = elz(nij)   + tab(4,it,jt) 
          end if
        end do
      end do

c     write (sor,*) 'atoma,atomb,la,lb, gama,gamb ',
c    &               i,j,lit,ljt,gama,gamb
       
      goto 8000
      write (sor,*) 
      write (sor,*) 
      write (sor,*) '1/r4 '
      write (sor,*) 
      do k1=1,lit*2+1
          write (sor,*) (tab(1,k1,k2), k2=1,ljt*2+1)
      end do
      write (sor,*) 
      write (sor,*) 
      write (sor,*) 'x/r3 '
      write (sor,*) 
      do k1=1,lit*2+1
          write (sor,*) (tab(2,k1,k2), k2=1,ljt*2+1)
      end do
      write (sor,*) 
      write (sor,*) 
      write (sor,*) 'y/r3 '
      write (sor,*) 
      do k1=1,lit*2+1
          write (sor,*) (tab(3,k1,k2), k2=1,ljt*2+1)
      end do
      write (sor,*) 
      write (sor,*) 
      write (sor,*) 'z/r3 '
      write (sor,*) 
      do k1=1,lit*2+1
          write (sor,*) (tab(4,k1,k2), k2=1,ljt*2+1)
      end do

 8000 continue
 9000 continue
c
      write(sor,7602)(i,elx(i),ely(i),elz(i),i=1,nnp)
c
      if(pl.eq.3)then
      write(sor,*) 'on ecrit : i, i1sr4'
      write(sor,7604)(i,i1sr4(i),i=1,nnp)
      write(sor,*) 'on ecrit : i, elx(i),ely(i),elz(i)'
 7602 format(/,4(i5,3f8.5))
 7604 format(/,4(i5,f8.5))
      endif
c
c
 9500 continue
      write(17) (elx(i),i=1,nnp)
      write(17) (ely(i),i=1,nnp)
      write(17) (elz(i),i=1,nnp)
      write(17) (i1sr4(i),i=1,nnp)
      do i=1,nnp
c      write (6,*) 'i, Ex Ey Ez 1/r**4 ',i,elx(i),ely(i),elz(i),i1sr4(i)
      end do
10000 continue
c
c
      write(sor,9400)
 9400 format(' ......  end of 1/r**4, end of electric field ......')
      return
      end
