      subroutine vector(cv)
      implicit real*8 (a-h,o-x,z),logical*1(y)
      include 'pshf.prm'
      character *4 iflab
      common norb,nsym,its(nsymz,nsymz),itsym(doa),itsyv(2*doa),
     1 nvec,isymp(nsymz),mijkl,nrec,mblk
      common isydeg(nsymz,nsymz),ntrsy,istop,iarwr(2*doa),
     1nocsy(nsymz),ymono,yprt,yprtv,yatom,ywvect,ynocsy
      common /om/ c(doa*doa*2),nbp(nsymz),nbo(nsymz),ind(nsymz+1),
     y ivecsy(2*doa*nsymz)
      common /hondo/ nt,ishell(doa),kloc(doa),ktyp(doa),
     1 numc(doa),newsh(doa,nsymz),iptr(3,nsymz),idtr(6,nsymz),
     2 iftr(10,nsymz),igtr(15,nsymz),yhondo,yprth
      common/nuscf/iscf(doa),iflab(3,doa)
      dimension ipr(doa*nsymz)
      dimension cv(doa*doa),yu(2*doa)
      do 1 i=1,2*doa
   1        yu(i)=.false.
         lsymi=0
          ii=0
         lve=0
         do 9 lsym=1,nsym
         ind(lsym)=lve+1
         k=0
        do 4 i=1,norb
          if(ivecsy(lsymi+i).eq.0) go to 4
         k=k+1
        ivecsy(lsymi+i)=k
    4        continue
      nbp(lsym)=k
         k=0
         do 6 i=1,nvec
         if(itsyv(i).ne.lsym) go to 6
      ni=(i-1)*norb
         k=k+1
        ii=ii+1
         iarwr(ii)=i
         do 5 j=1,norb
          if(ivecsy(lsymi+j).eq.0) go to 5
         lve=lve+1
         c(lve)=cv(ni+j)
    5         continue
   6       continue
           nbo(lsym)=k
9         lsymi=lsymi+norb
        ind(nsym+1)=lve+1
      nsym1=nsym+1
       nboo=0
        maxx=8
        iinf=0
         kw=0
        do 20 ii=1,nsym
c
        k=0
        do 22 i=1,norb
        if(ivecsy((ii-1)*norb+i).ne.0) then
        k=k+1
        ipr(k)=i
        end if
   22  continue
        iw1=ind(ii)-1
         iw2=ind(ii+1)-1
        if(iw1.ge.iw2) go to 20
        nboo=nboo+nbo(ii)
         write(6,14) ii
 69        kw1=kw+1
          kw=kw+maxx
        if(kw.gt.nboo) kw=nboo
        isup=kw-kw1+1
        write(6,76) (iarwr(iw),iw=kw1,kw)
        write(6,77) (iscf(iw),iw=kw1,kw)
 76       format(14x,8i9)
 77       format(18x,1h(,7(i3,1h),4x,1h(),i3,1h))
      jj=nbp(ii)
70         format(2x,3a4,1x,8f9.4,7x)
         if(.not.yhondo) jj=nbo(ii)
14        format('0om de la symetrie:',i3)
          do 15 i=1,jj
15        write(6,70)(iflab(l,ipr(i)),l=1,3),
     #      (c(iinf+(iw-1)*jj+i),iw=1,isup)
         iinf=iinf+jj*isup
         if(kw.lt.nboo) go to 69
20       continue
      if(yprtv) write(6,60)(nbo(i),nbp(i),i=1,nsym)
   60    format(' nbo,nbp',(1x,40i3))
         return
      end
