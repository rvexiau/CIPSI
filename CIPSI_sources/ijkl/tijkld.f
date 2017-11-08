      subroutine tijkld(mi,mj,mk,ml,ic,jc,kc,lc,yij,ykl,yijkl,
     &mp,mq,ypq,bijkld)
      implicit real*8 (a-h,o-x,z),logical*1 (y)
      include 'pshf.prm'
      integer*4 p,q,pq,qp,qj,pi
c version double precision
      real*4 bufs1,bufs2
      dimension bijkld(*)
      common norb,nsym,its(nsymz,nsymz),itsym(doa),itsyv(2*doa),
     1 nvec,isymp(nsymz),mijkl,nrec,mblk
      common isydeg(nsymz,nsymz),ntrsy,istop,iarwr(2*doa),
     1 nocsy(nsymz),ymono,yprt,yprtv,yatom,ywvect,ynocsy
      common /om/c(doa*doa*2), nbp(nsymz),nbo(nsymz),ind(nsymz+1),
     y ivecsy(2*doa*nsymz)
      common /hondo/ nt,ishell(doa),kloc(doa),ktyp(doa),
     1 numc(doa),newsh(doa,nsymz),iptr(3,nsymz),idtr(6,nsymz),
     2 iftr(10,nsymz),igtr(15,nsymz),yhondo,yprth
       common/tampon/bufs1(lsize),bufs2(lsize),bufd1(lsize),
     1 bufd2(lsize),monowr
      common/tabint/tab(doa*doa),tac(doa)
      numnu(i)=(i*(i-1))/2
c      
      iw1=mijkl+1
      nk=mk
      ni=mi
      nq=mq
      ns=ms
      mkl=mk*ml
      if(ykl) mkl=numnu(mk+1)
      ilk=0
      do 300 l=1,ml
      if(ykl) nk=l
      do 300 k=1,nk
      ilk=ilk+1
      nkl=0
      qp=0
      do 230 p=1,mp
      if(ypq) nq=p
      do 220 q=1,nq
      sum=bijkld(nkl+ilk)

      nkl=nkl+mkl
220   tab(qp+q)=sum
      if(.not.ypq) go to 230
      pq=0
      do 225 q=1,nq
      tab(pq+p)=tab(qp+q)
225   pq=pq+mp
230   qp=qp+mq
c
      nij=0
      qj=jc-1
      do 280 j=1,mj
      qp=0
      do 250 p=1,mp
      sum=0.d0
      do 240 q=1,mq
240   sum=sum+ c(qj+q)*tab(qp+q)
      qp=qp+mq
250   tac(p)=sum
      if(yij) ni=j
      pi=ic-1
      do 270 i=1,ni
      nij=nij+1
      sum=0.d0
      do 260 p=1,mp
260   sum=sum+ c(pi+p)*tac(p)
      if(mijkl.lt.lsize) go to 265
      mblk=mblk+1
c      write(*,*) 'ecriture'      
      write(50) bufd1
      mijkl=0
265   mijkl=mijkl+1
      bufd1(mijkl)=sum
      if(.not.yijkl) go to 270
      if(nij.eq.ilk) go to 300
270   pi=pi+mp
280   qj=qj+mq
300   continue
c
700    format(1x,15f8.4)
c
      return
      entry wtrbld
      mblk=mblk+1
       write(*,*) 'ecriture'       
      write(50) bufd1
      return
c
      end
