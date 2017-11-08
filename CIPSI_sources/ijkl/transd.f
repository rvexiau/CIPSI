      subroutine transd(msi,msj,msk,msl,bijkld,xx_pqrs,ix_pqrs)
      implicit real*8 (a-h,o-x,z),logical*1 (y)
      include 'pshf.prm'
c version double precision
      real*4 bufs1,bufs2       
      dimension bijkld(*)
      logical llg1,llg2,llg3,llg4
      common norb,nsym,its(nsymz,nsymz),itsym(doa),itsyv(2*doa),
     1 nvec,isymp(nsymz),mijkl,nrec,mblk
      common/om/cvec(doa*doa*2),nbp(nsymz),nbo(nsymz),
     1 indiv(nsymz+1),ivecsy(2*doa*nsymz)
      common /hondo/ nt,ishell(doa),kloc(doa),ktyp(doa),
     1 numc(doa),newsh(doa,nsymz),iptr(3,nsymz),idtr(6,nsymz),
     2 iftr(10,nsymz),igtr(15,nsymz),yhondo,yprth
       common/tampon/bufs1(lsize),bufs2(lsize),bufd1(lsize),
     1 bufd2(lsize),monowr
      dimension xx_pqrs(*),ix_pqrs(*)
      dimension isi(nsymz),mti(nsymz),mtj(nsymz),mtk(nsymz),mtl(nsymz),
     * ysy(nsymz)
      dimension xx(lsize),ix(lsize)
      dimension yl(4)
      equivalence (lable,yl(1)),(llg1,i),(llg2,j),(llg3,k),(llg4,l)
      dimension nct(4),ncn(4)
      dimension ylsy(4),ymsy(4)
      equivalence (ylsy(1),labsy),(ymsy(1),mabsy)
      equivalence (nct(4),lct),(nct(3),kct),(nct(2),jct),(nct(1),ict),
     &            (ncn(4),lcn),(ncn(3),kcn),(ncn(2),jcn),(ncn(1),icn)
      dimension yis(doa),yjs(doa),yks(doa),yls(doa),yij(doa),
     1 ykl(doa)
      dimension yijkls(doa),nspi(doa),nspj(doa),nspk(doa),
     1 nspl(doa)
      ia(i)=(i*(i-1))/2            
                  
      nbpj=nbp(msj)
      nbpl=nbp(msl)
      nbpkl=nbp(msk)*nbp(msl)
      nbpkls=(nbpl+nbpkl)/2
      jsyi=(msj-1)*norb
      ksyi=(msk-1)*norb
      lsyi=(msl-1)*norb
      isyi=(msi-1)*norb
      ki=0
      kj=0
      kk=0
      kl=0
      
      do 850 i=1,norb
      yis(i)=ivecsy(isyi+i).eq.0
      yjs(i)=ivecsy(jsyi+i).eq.0
       yks(i)=ivecsy(ksyi+i).eq.0
       yls(i)=ivecsy(lsyi+i).eq.0
      yij (i)=yis(i).and.yjs(i)
      ykl (i)=yks(i).and.yls(i)
       if(yis(i)) go to 810
       ki=ki+1
       nspi(i)=ki
810    if(yjs(i)) go to 820
        kj=kj+1
        nspj(i)=kj
820     if(yks(i)) go to 830
        kk=kk+1
         nspk(i)=kk
830   if(yls(i)) go to 840
         kl=kl+1
         nspl(i)=kl
840   yijkls(i)=yij(i).and.ykl(i)
850   continue
      ypq=isymp(msi).eq.isymp(msj)
      yrs=isymp(msk).eq.isymp(msl)
      ypqrs=isymp(msi).eq.isymp(msk).and.isymp(msj).eq.isymp(msl)
      ypqrst=ypq.and.ypqrs
       ijcas=1
       if(ypq) ijcas=ijcas+1
       if(yrs)  ijcas=ijcas+2
              dfnt=1.d0/dfloat(nt)
c si ijkl est nul a la fin de la routine il n y a pas d int  de ce type
      ijkl=0
      qq4=1.d0
      rewind 8
      icd=0
      jcd=0
      kcd=0
      lcd=0
      imin=-lsize+1
10    imin=imin+lsize
      imax=imin+lsize-1
      xx=xx_pqrs(imin:imax)
      ix=ix_pqrs(imin:imax)
c      call read8(xx,ix,lsize)

      do 300 m=1,lsize
      lable=ix(m)
      if(lable.eq.0) go to 350
      i=iword1(lable)
      j=iword2(lable)
      k=iword3(lable)
      l=iword4(lable)
        if(yijkls(i)) go to 300
         if(yijkls(j)) go to 300
c     orbitales en ordre canonique pas les couches
           if(yijkls(k)) go to 300
            if(yijkls(l)) go to 300
      if(ypqrst) go to 190
      if(.not.(yij(i).or.yij(j))) go to 135
      if(ykl(i).or.ykl(j)) go to 300
      if(yij(k).or.yij(l)) go to 300
      icas=5
      go to 200
135   if(.not.(ykl(k).or.ykl(l))) go to 140
      if(yij(k).or.yij(l)) go to 300
      if(ykl(i).or.ykl(j)) go to 300
      icas=5
      go to 200
140     if(ypqrs) go to 160
       if(yij(k).or.yij(l)) go to 150
       if(ykl(i).or.ykl(j)) go to 150
      icas=6
      go to 200
150   icas=ijcas
      go to 200
160   icas=7
      go to 200
190    icas=8
200    continue
      if(nt.lt.2) go to 230
      ic=ishell(i)
      kc=ishell(k)
      lc=ishell(l)
      jc=ishell(j)
      ig=ic
      jg=jc
      kg=kc
      lg=lc
      if(ic.gt.jc) go to 110
      jc=ig
      ic=jg
110   if(kc.gt.lc) go to 111
      kc=lg
      lc=kg
111   if(ic-kc) 114,115,116
114   call ex(ic,kc,jc,lc)
c
      go to 116
115   if(lc.gt.jc) go to 114
116   continue
      ict=ktyp(i)
      jct=ktyp(j)
      kct=ktyp(k)
      lct=ktyp(l)
      icn=numc(i)
      jcn=numc(j)
      kcn=numc(k)
      lcn=numc(l)
      if(lc.ne.lcd.or.kc.ne.kcd.or.jc.ne.jcd.or.ic.ne.icd) go to 400
      go to 230
400   continue
      nq4=1
      isi(1)=1
      ngen=0
c
      do 500 n=2,nt
      ysy(n)=.true.
      it=newsh(ic,n)
      jt=newsh(jc,n)
      kt=newsh(kc,n)
      lt=newsh(lc,n)
      if(it.lt.jt) go to 410
      n1=it
      n2=jt
      go to 415
410   n1=jt
      n2=it
415   if(kt.lt.lt) go to 420
      n3=kt
      n4=lt
      go to 425
420   n3=lt
      n4=kt
425   if(n1-n3) 430,435,440
  435 if(n2.ge.n4) go to 440
  430 nn=n1
      n1=n3
      n3=nn
      nn=n2
      n2=n4
      n4=nn
  440 continue
      if(n1.ne.ic.or.n2.ne.jc.or.n3.ne.kc.or.n4.ne.lc) go to 460
      nq4=nq4+1
      ysy(n)=.false.
      go to 500
  460 if(ngen.eq.0) go to 800
      do 802 kx=1,ngen
      if(n1.ne.mti(kx).or.n2.ne.mtj(kx).or.n3.ne.mtk(kx).or.
     * n4.ne.mtl(kx)) go to 802
      ysy(n)=.false.
      go to 500
  802 continue
  800 ngen=ngen+1
      mti(ngen)=n1
      mtj(ngen)=n2
      mtk(ngen)=n3
      mtl(ngen)=n4
500   continue
      lcd=lc
      kcd=kc
      jcd=jc
      icd=ic
      qq4=dfloat(nq4)*dfnt
230   val=xx(m)*qq4
      yijg=i.eq.j
      yijgr=yijg
      if(yijg) val=val+val
      yklg=k.eq.l
      yklgr=yklg
      if(yklg) val=val+val
      yijklg=i.eq.k.and.j.eq.l
      if(yijklg)val=val+val
      ind=0
       vali=val
235   if(ind.ge.nt) go to 300
      ind=ind+1
      ysig=.false.
      if(ind.eq.1) go to 260
      if(.not.ysy(ind)) go to 235
      do 490 mn=1,4
      if(nct(mn).eq.1) go to 490
      mcn=ncn(mn)+1
      if(nct(mn)-3) 470,475,480
470   if(iptr(mcn,ind).lt.0) ysig=.not.ysig
      go to 490
475   if(idtr(mcn,ind).lt.0) ysig=.not.ysig
      go to 490
480   if(iftr(mcn,ind).lt.0) ysig=.not.ysig
490   continue
        val=vali
        if(ysig)val=-vali
      i=kloc(newsh(ig,ind))+icn
      j=kloc(newsh(jg,ind))+jcn
      k=kloc(newsh(kg,ind))+kcn
      l=kloc(newsh(lg,ind))+lcn
      yijg=yijgr
      yklg=yklgr
6566     format(1x,4i4,f20.8,i10,9i4)
260   yflem=.false.
        go to (900,910,920,930,940,980,990,1000),icas
900   ij2=0
      if(yis(i).or.yjs(j)) go to 901
      ij1=nbpj*(nspi(i)-1)+nspj(j)
      if(yis(j).or.yjs(i).or.yijg) go to 902
      ij2=nbpj*(nspi(j)-1)+nspj(i)
      go to 902
901   if(yis(j).or.yjs(i)) go to 1050
      ij1=nbpj*(nspi(j)-1)+nspj(i)
902   if(yks(k).or.yls(l)) go to 903
       kl1=nbpl*(nspk(k)-1)+nspl(l)
       ijkl=nbpkl*(ij1-1)+kl1
      bijkld(ijkl)=val
      if(ij2.eq.0) go to 903
       ijkl=nbpkl*(ij2-1)+kl1
      bijkld(ijkl)=val
903   if(yks(l).or.yls(k).or.yklg) go to 1050
904     kl1=nbpl*(nspk(l)-1)+nspl(k)
       ijkl=nbpkl*(ij1-1)+kl1
      bijkld(ijkl)=val
      if(ij2.eq.0) go to 1050
       ijkl=nbpkl*(ij2-1)+kl1
      bijkld(ijkl)=val
      go to 1050
c i ou j sur 1  k ou l sur 2
910     if(i.ge.j) go to 911
      ij1=ia(nspj(j))+nspi(i)
      go to 912
911   ij1=ia(nspi(i))+nspj(j)
912   if(yks(k).or.yls(l)) go to 914
       kl1=nbpl*(nspk(k)-1)+nspl(l)
      ijkl=nbpkl*(ij1-1)+kl1
      bijkld(ijkl)=val
      if(yklg) go to 1050
914   if(yls(k).or.yks(l)) go to 1050
      kl1=nbpl*(nspk(l)-1)+nspl(k)
      ijkl=nbpkl*(ij1-1)+kl1
      bijkld(ijkl)=val
      go to 1050
c ks=ls  is ne js
920   if(k.ge.l) go to 921
      kl1=ia(nspl(l))+nspk(k)
      go to 922
921   kl1=ia(nspk(k))+nspl(l)
922   if(yis(i).or.yjs(j)) go to 924
      ij1=nbpj*(nspi(i)-1)+nspj(j)
      ijkl=nbpkls*(ij1-1)+kl1
      bijkld(ijkl)=val
      if(yijg) go to 1050
924   if(yis(j).or.yjs(i)) go to 1050
      ij1=nbpj*(nspi(j)-1)+nspj(i)
      ijkl=nbpkls*(ij1-1)+kl1
      bijkld(ijkl)=val
      go to 1050
c is=js ks=ls ********
930   if(i.ge.j) go to 931
      ij1=ia(nspj(j))+nspi(i)
      go to 932
931   ij1=ia(nspi(i))+nspj(j)
932   if(k.ge.l) go to 934
      kl1=ia(nspl(l))+nspk(k)
      go to 935
934    kl1=ia(nspk(k))+nspl(l)
935   ijkl=nbpkls*(ij1-1)+kl1
      bijkld(ijkl)=val
      go to 1050
c  i,j sur ks ls  k,l sur is js
940      call ex(i,k,j,l)
         yflem=.false.
      ycon=yijg
      yijg=yklg
      yklg=ycon
        go to (900,910,920,930),ijcas
c i,j sur is js ks ls  k,l sur is js ks ls
980       yflem=.true.
      if(yijklg) yflem=.false.
        go to (900,910,920,930),ijcas
1050     if(yflem) go to 940
          go to 235
c  is=ks js=ls is#js
990    yflem=.true.
      if(yis(i).or.yjs(j)) go to 999
991   ij1=nbpj*(nspi(i)-1)+nspj(j)
      if(yks(k).or.yls(l)) go to 995
      kl1=nbpl*(nspk(k)-1)+nspl(l)
      if(k-i) 992,993,994
992   ijkl=ij1*(ij1-1)/2+kl1
      bijkld(ijkl)=val
      if(yklg) go to 999
      go to 995
993   if(j.ge.l) go to 992
994   ijkl=kl1*(kl1-1)/2+ij1
      bijkld(ijkl)=val
      if(yklg) go to 999
995   if(yks(l).or.yls(k)) go to 999
      if(.not.yflem.and.yijklg) go to 235
       kl1=nbpl*(nspk(l)-1)+nspl(k)
      if(l-i) 996,997,998
996     ijkl=ij1*(ij1-1)/2+kl1
      bijkld(ijkl)=val
      go to 999
997      if(j.ge.k) go to 996
998       ijkl=kl1*(kl1-1)/2+ij1
      bijkld(ijkl)=val
999   if(.not.yflem.or.yijg) go to 235
         yflem=.false.
      if(yis(j).or.yjs(i)) go to 235
       call ex(i,j,k,l)
       go to 991
c is=js=ks=ls
1000      if(ind.eq.1) go to 1010
          if(i.ge.j) go to 1001
           is=i
      i=j
      j=is
1001    if(k.ge.l) go to 1002
         ks=k
         k=l
       l=ks
1002         if(i-k) 1003,1004,1010
1003     call ex(i,k,j,l)
          go to 1010
1004      if(j.lt.l) go to 1003
1010       ij1=ia(nspi(i))+nspj(j)
           kl1=ia(nspk(k))+nspl(l)
         ijkl=ij1*(ij1-1)/2+kl1
      bijkld(ijkl)=val
      go to 235
300   continue
      go to 10
350   if(ijkl.eq.0) msi=0

      return
      end
      
