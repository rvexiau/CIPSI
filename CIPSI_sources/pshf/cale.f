      subroutine cale(f,h,v,e,ipass,ending,uhf,ia,in)
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      parameter (itmm=5)
      logical ending,uhf
      character*8 type
      character*8 aname,bname
      common/iofile/ir,iw,ip,is,iq,ih,iv
      common/infoa/nat,ich,mul,numscf,nx,ne,na,nb,zan(dc),c(3,dc)
      dimension f(*),v(*),e(*),h(*),ia(*),in(*)
       common/conv/acurcy,ehf,ehf0,iter,idiis,ijump
      common/scfit/ aname,bname,maxit,nconv,npunch,npi,cextrp,amix
      common/extrp/ b(20,20),d(20),lrec,nxr,iodb(5000),fmab
      dimension bb(itmm,itmm),type(2)
      data type/' alpha  ','  beta  '/
c
      cdiver=cextrp
      if(uhf) cdiver=cdiver*0.5d0
      lrec=(nx+lsize-1)/lsize
c
cvax
cvax  open (unit=iq,access='direct',recl=2*lsize,status='scratch')
cvax
      open (unit=iq,access='direct',recl=8*lsize,status='scratch')
c    *      file='fort.9')
      iodb(1)=1
      nav1=1
      nav2=2
      nav3=5
      iodb(iter+1)=iodb(iter)+2*lrec
      if(.not.uhf) go to 6
      nav2=3
      if(ipass.eq.1) go to 6
      iodb(iter+1)=iodb(iter+1)+2*lrec
      nav1=2
      nav2=4
      nav3=8
    6 continue
c
      nxsq=numscf*numscf
c    v contient f, f contient d.
c....calcul de d*f
c
      cmix=cdiver*10.d0
      call wrtdb (f,nx,iter,nav1)
      do 60 i=1,numscf
      do 40 j=1,numscf
      dum=0.d0
      do 30 k=1,numscf
      jk=ia(j)+k
      if(j.lt.k) jk=ia(k)+j
      ki=ia(i)+k
      if(i.lt.k) ki=ia(k)+i
 30   dum=dum+h(jk)*f(ki)
 40   e(j)=dum
      do 50 j=1,numscf
 50   v(in(j)+i)=e(j)
 60   continue
c...calcul de s*d*f
      call reada(f,nx,2)
      do 160 i=1,numscf
      do 140 j=1,numscf
      dum=0.d0
      do 130 k=1,numscf
      jk=ia(j)+k
      if(j.lt.k) jk=ia(k)+j
 130  dum=dum+f(jk)*v(in(k)+i)
 140  e(j)=dum
      do 150 j=1,numscf
 150  v(in(j)+i)=e(j)
 160  continue
c
c....calcul de e=fds-sdf
      do 170 i=1,numscf
      do 170 j=1,i
      ij=ia(i)+j
      f(ij)=v(in(j)+i)-v(in(i)+j)
 170  continue
      if(ipass.eq.1) fmab=0.d0
      fmax=0.d0
      do 172 i=1,nx
       if(abs(f(i)).gt.fmax) fmax=abs(f(i))
 172  continue
      write(iw,9998) type(ipass),fmax
      if(fmax.gt.fmab) fmab=fmax
 9998 format(/' maximum of error vector for ',a8,'set',f20.8)
c
c...calcul de e*a
      call reada(h,nx,4)
      do 245 i=1,numscf
      do 240 j=1,numscf
      dum=0.d0
      do 230 k=1,numscf
      jk=ia(j)+k
      if(j.lt.k) jk=ia(k)+j
      ki=ia(i)+k
      if(i.lt.k) ki=ia(k)+i
 230  dum=dum+f(ki)*h(jk)
 240  e(j)=dum
      do 245 j=1,numscf
 245  v(in(i)+j)=e(j)
c...calcul de a+*e*a
      do 255 i=1,numscf
      do 250 j=1,numscf
      dum=0.d0
      do 260 k=1,numscf
      kj=ia(k)+j
      if(k.lt.j) kj=ia(j)+k
 260  dum=dum+h(kj)*v(in(k)+i)
 250  e(j)=dum
      do 255 j=1,numscf
 255  v(in(j)+i)=e(j)
      ij=0
      do 280 i=1,numscf
      do 280 j=1,i
      ij=ij+1
 280  f(ij)=v(in(i)+j)
      call wrtdb(f,nx,iter,nav2)
c     si calcul uhf et premier passage  fin ici
      if(ipass.eq.1.and.uhf) return
c
c     extrapolation ou interpolation suivant fmab
      if(fmab.lt.acurcy) ending=.true.
      if(iter.eq.1) return
      if(fmab.lt.cdiver) go to 179
c
      write(iw,9996) amix
 9996 format(' error vector too large... damping with amix=',f12.4)
      call readb(f,nx,iter-1,1)
      call readb(v,nx,iter,1)
c     write (6,*) 'f(i),v(i) en cale', (f(i),v(i), i= 1,10)
      alp=1.d0-amix
      do 176 i=1,nx
  176 f(i)=amix*f(i)+alp*v(i)
c     write (6,*) 'f(i),v(i) en cale desp de mix', (f(i),v(i), i= 1,10)
      call wrtda(f,nx,5)
      call wrtdb (f,nx,iter,nav1)
      if(.not.uhf) return
      call readb(f,nx,iter-1,2)
      call readb(v,nx,iter,2)
      do 177 i=1,nx
  177 f(i)=amix*f(i)+alp*v(i)
      call wrtda(f,nx,8)
      return
  179 idiis=idiis+1
c....calcul de b et de c (dans d)
      b(1,1)=0.d0
c
      itact=idiis-ijump
      if(itact.lt.itmm) go to 281
      itact=itmm-1
      ijump=ijump+1
      do 284 i=3,itmm
      do 284 j=1,itmm
      b(i-1,j)=b(i,j)
  284 continue
      itmm1=itmm-1
      do 285 i=1,itmm1
      do 285 j=3,itmm
      b(i,j-1)=b(i,j)
  285 continue
  281 continue
      ipss=1
      nav2=2
      if(uhf) nav2=3
  282 continue
      call readb(h,nx,iter,nav2)
      d(1)=-1.d0
      do 300 it=1,itact
      d(it+1)=0.d0
      call readb(f,nx,iter+it-itact,nav2)
      if(ipss.eq.1) b(it+1,itact+1)=0.d0
      b(it+1,itact+1)=b(it+1,itact+1)+tracep(f,h,numscf)
      b(itact+1,it+1)=b(it+1,itact+1)
  300 continue
      b(1,itact+1)=-1.d0
      b(itact+1,1)=-1.d0
      d(itact+1)=0.d0
      if(ipss.eq.2.or..not.uhf) go to 305
      ipss=2
      nav2=4
      go to 282
  305 it1=itact+1
c
      if(itact.eq.1) go to 345
c
      do 302 i=1,it1
      do 302 j=1,it1
  302 bb(i,j)=b(i,j)
      call dresyl(bb,itmm,itmm,d,itmm,1,it1,1,1.d-10,kkk)
c
 9999 format(' coefficients for extrapolation',8f10.6,/,(10f10.6))
      if(kkk.ne.1) goto 315
      if(itact.eq.1) goto 312
      d(1)=-1.d0
      do 306 i=2,itact
      d(i)=0.d0
      do 306 j=2,it1
  306 b(i,j)=b(i+1,j)
      do 308 i=2,itact
      do 308 j=2,itact
  308 b(i,j)=b(i,j+1)
      ijump=ijump+1
      itact=itact-1
      goto 305
 312  write(iw,9997)
 9997 format(/' matrice b reduite a un seul element')
 315  continue
      write(iw,9999) (d(i),i=2,it1)
c    calcul de la nouvelle matrice de fock
      ipss=1
      nav1=1
      nav3=5
  318 continue
      do 320 i=1,nx
  320 f(i)=0.d0
      do 330 it=1,itact
      call readb(v,nx,iter+it-itact,nav1)
      do 340 i=1,nx
  340 f(i)=f(i)+d(it+1)*v(i)
  330 continue
c
c
      call wrtda (f,nx,nav3)
  345 continue
      if(ipss.eq.2.or..not.uhf) return
      ipss=2
      nav1=2
      nav2=4
      nav3=8
      go to 318
      end
