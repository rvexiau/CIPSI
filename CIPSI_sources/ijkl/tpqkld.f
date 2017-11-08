      subroutine tpqkld(mp,mq,mr,ms,mk,ml,ypq,yrs,ypqrs,ykl,yijkl,kc,lc
     y,bijkld)
      implicit real*8 (a-h,o-x,z),logical*1 (y)
      include 'pshf.prm'
      integer*4 p,q,r,s,rs,sr,sl,rk
c   version double precision
      real*4 bufs1,bufs2
      common norb,nsym,its(nsymz,nsymz),itsym(doa),itsyv(2*doa),
     1 nvec,isymp(nsymz),mijkl,nrec,mblk
      common /om/ c(doa*doa*2),nbp(nsymz),nbo(nsymz),ind(nsymz+1),
     y ivecsy(2*doa*nsymz)
       common/tampon/bufs1(lsize),bufs2(lsize),bufd1(lsize),
     1 bufd2(lsize),monowr
      common/tabint/tab(doa*doa),tac(doa)
      dimension bijkld(*)
      numnu(i)=(i*(i-1))/2
       nk=mk
      nl=ml
      mono=monowr
      nq=mq
        ns=ms
      ijkl=0
       nr=mr
       np=mp
      npq=0
      mpqm=np*mq
      if(ypq) mpqm=numnu(np+1)
      do 200 p=1,np
      if(ypq) go to 5
      mpq=0
      go to 6
5     nq=p
6     do 200 q=1,nq
      if(ypq) go to 10
      npq=mpq+p
      mpq=mpq+np
      go to 15
10    npq=npq+1
15    nrs=0
      mij=0
      sr=0
      do 110 r=1,nr
      if(yrs) go to 20
      mrs=0
      go to 21
20    ns=r
21    rs=0
      do 100 s=1,ns
      if(yrs) go to 60
      nrs=mrs+r
      mrs=mrs+nr
      if(ypqrs) go to 30
      mij=mpqm*(nrs-1)
      tab(sr+s)=bijkld(mij+npq)
c
      go to 100
30    if(npq.lt.nrs) go to 40
      ijkl=numnu(npq)+nrs
      tab(sr+s)=bijkld(ijkl)
      go to 100
40    tab(sr+s)=bijkld(numnu(nrs)+npq)
      go to 100
60    nrs=nrs+1
      if(ypqrs) go to 70
      pqrs=bijkld(mij+npq)
      mij=mij+mpqm
      tab(rs+r)=pqrs
      tab(sr+s)=pqrs
      go to 100
70    if(npq.lt.nrs) go to 80
      ijkl=ijkl+1
      pqrs=bijkld(ijkl)
      tab(sr+s)=pqrs
      tab(rs+r)=pqrs
      go to 100
80    pqrs=bijkld(numnu(nrs)+npq)
      tab(sr+s)=pqrs
      tab(rs+r)=pqrs
100   rs=rs+nr
110   sr=sr+ms
c
      sl=lc-1
      do 170 l=1,nl
      sr=0
      do 130 r=1,nr
      sum=0.d0
      do 120 s=1,ms
120   sum=sum+ c(sl+s)*tab(sr+s)
      sr=sr+ms
130   tac(r)=sum
      if(ykl) nk=l
      rk=kc-1
      do 160                                               k=1,nk
      sum=0.d0
      do 150 r=1,nr
150   sum=sum+ c(rk+r)*tac(r)
      rk=rk+nr
      if(mono.lt.lsize) go to 155
       write(40) bufd2
       mono=0
155   mono=mono+1
160        bufd2(mono)=sum
170   sl=sl+ms
200   continue
700    format(1x,15f8.4)
      monowr=mono
      return
      end
