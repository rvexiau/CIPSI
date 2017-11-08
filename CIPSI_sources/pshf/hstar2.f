      subroutine hstar2 (f1,da1,db1,f2,da2,db2,f3,da3,db3,ia,in,nc)
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      logical*1 log
      logical exch
c
c     form a skeleton fock matrix : hstar
c
      common/restar/timlim,irest,nrec,intloc,ist,jst,kst,lst
      common/iofile/ir,iw,ip,is
      common/infoa/nat,ich,mul,numscf,nx,ne,na,nb,zan(dc),c(3,dc)
      dimension da1(*),f1(*),db1(*),f2(*),da2(*),db2(*),ia(*),in(*)
      dimension f3(*),da3(*),db3(*)
      dimension ix(lsize),xx(lsize)
      dimension log(4)
      equivalence (log(1),lable)
      equivalence (llg1,i),(llg2,j),(llg3,k),(llg4,l)
      data two /2.d0/
                     
      do 50 m=1,nx
      f1(m)=0.d0
      f2(m)=0.d0
   50 f3(m)=0.d0
      i=0
      j=0
      k=0
      l=0
      ntwol=0
c
c     read in one record of integrals
c
  100 read(is) xx,ix
c
      do 300 m=1,lsize
      lable=ix(m)
      if(lable.eq.0) go to 350
      i=iword1(lable)
      j=iword2(lable)
      k=iword3(lable)
      l=iword4(lable)
c      write(6,'(2i12,4i4,e14.6)') ,m,label,i,j,k,l,xx(m)

      val=xx(m)
      val2=val+val
      val4=val2+val2
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
  250 dum1=da1(nkl)*val4
      dum2=da1(nij)*val4
      f1(nij)=f1(nij)+dum1
      f1(nkl)=f1(nkl)+dum2
      f1(nik)=f1(nik)-val2*db1(njl)
      f1(nil)=f1(nil)-val2*db1(njk)
      f1(njk)=f1(njk)-val2*db1(nil)
      f1(njl)=f1(njl)-val2*db1(nik)
      if(nc.eq.1) go to 300
      dum1=da2(nkl)*val4
      dum2=da2(nij)*val4
      f2(nij)=f2(nij)+dum1
      f2(nkl)=f2(nkl)+dum2
      f2(nik)=f2(nik)-val2*db2(njl)
      f2(nil)=f2(nil)-val2*db2(njk)
      f2(njk)=f2(njk)-val2*db2(nil)
      f2(njl)=f2(njl)-val2*db2(nik)
      if(nc.lt.3) go to 300
      dum1=da3(nkl)*val4
      dum2=da3(nij)*val4
      f3(nij)=f3(nij)+dum1
      f3(nkl)=f3(nkl)+dum2
      f3(nik)=f3(nik)-val2*db3(njl)
      f3(nil)=f3(nil)-val2*db3(njk)
      f3(njk)=f3(njk)-val2*db3(nil)
      f3(njl)=f3(njl)-val2*db3(nik)
  300 continue
      go to 100
  350 continue
      if(numscf.eq.1) return
      do 400 m=2,numscf
      max=m-1
      do 400 n=1,max
      nij=ia(m)+n
      f1(nij)=f1(nij)/two
      f2(nij)=f2(nij)/two
  400 f3(nij)=f3(nij)/two
      rewind is
c
      return
      end
