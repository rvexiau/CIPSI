      subroutine hstar (da,fa,db,fb,ia,uhf,xx_pqrs,ix_pqrs)
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      dimension xx_pqrs(*),ix_pqrs(*)      
      logical exch
      logical uhf
c
c     form a skeleton fock matrix : hstar
c
      common/restar/timlim,irest,nrec,intloc,ist,jst,kst,lst
      common/iofile/ir,iw,ip,is
      dimension da(*),fa(*),db(*),fb(*),ia(*)
      dimension ix(lsize),xx(lsize)
      common/infoa/nat,ich,mul,numscf,nx,ne,na,nb,zan(dc),c(3,dc)
      data two /2.d0/
                
      do 50 m=1,nx
      fa(m)=0.d0
   50 fb(m)=0.d0
      i=0
      j=0
      k=0
      l=0
      ntwol=0
      imin=1
      imax=lsize        
c
c     read in one record of integrals
c
c 100 read(is) xx,ix
  100 xx=xx_pqrs(imin:imax)
      ix=ix_pqrs(imin:imax)
      imin=imin+lsize
      imax=imin+lsize-1
c
      do 300 m=1,lsize
      lable=ix(m)
      if(lable.eq.0) go to 350
      i=iword1(lable)
      j=iword2(lable)
      k=iword3(lable)
      l=iword4(lable)
c     write(6,'(2i12,4i4,e14.6)') ,m,label,i,j,k,l,xx(m)
      val=xx(m)
      val2=val+val
      val4=val2+val2
      if(.not.uhf) val2=val
      nij=ia(i)+j
      nkl=ia(k)+l
      nik=ia(i)+k
      nil=ia(i)+l
      if(j.lt.k) go to 150
      njk=ia(j)+k
      njl=ia(j)+l
      go to 250
  150 njk=ia(k)+j
      if(j.lt.l) go to 200
      njl=ia(j)+l
      go to 250
  200 njl=ia(l)+j
  250 dum1=(da(nkl)+db(nkl))*val4
      dum2=(da(nij)+db(nij))*val4
c      pour alpha partie coulombienne
      fa(nij)=fa(nij)+dum1
      fa(nkl)=fa(nkl)+dum2
c     pour alpha partie d echange
      fa(nik)=fa(nik)-val2*da(njl)
      fa(nil)=fa(nil)-val2*da(njk)
      fa(njk)=fa(njk)-val2*da(nil)
      fa(njl)=fa(njl)-val2*da(nik)
      if(.not.uhf) go to 300
      fb(nij)=fb(nij)+dum1
      fb(nkl)=fb(nkl)+dum2
      fb(nik)=fb(nik)-val2*db(njl)
      fb(nil)=fb(nil)-val2*db(njk)
      fb(njk)=fb(njk)-val2*db(nil)
      fb(njl)=fb(njl)-val2*db(nik)
  300 continue
      go to 100
  350 continue
      if(numscf.eq.1) return
      do 400 m=2,numscf
      max=m-1
      do 400 n=1,max
      nij=ia(m)+n
      fa(nij)=fa(nij)/two
  400 fb(nij)=fb(nij)/two
c      rewind is
c
      return
      end
