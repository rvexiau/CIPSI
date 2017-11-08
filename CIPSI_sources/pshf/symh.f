      subroutine symh(f,h,ia)
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      logical*1 iandj
c
c     symmetrize the skeleton fock matrix
c
      dimension f(*),h(*),ia(*)
      common/nshel/ex(dgp),cs(dgp),cp(dgp),cd(dgp),cf(dgp),cg(dgp),
     1kstart(ds),katom(ds),ktype(ds),kng(ds),kloc(ds),kmin(ds),kmax(ds),
     2nshell,kmaz(ds)
      common/symtry/invt(48),iso(ds,12),ptr(3,144),dtr(6,288)
     1,ftr(10,480),nt
      common/output/nprint,itol,icut,normf,normp
      common/infoa/nat,ich,mul,numscf,nx,ne,na,nb,zan(dc),c(3,dc)
      common/hsym/t(35,35),mini,maxi,lit,minj,maxj,ljt,ntr
      common/isopac/indin(48),indout(12)
      dimension mi(48)
c      if(nt.eq.1) return
      ntwd=(nt+3)/4
      do 10 i=1,nx
   10 h(i)=0.0d+00
c
c     ----- find a block (i,j)
c
      do 2000 ii=1,nshell
      do 20 itr=1,ntwd
   20 indout(itr)=iso(ii,itr)
      call isoout(nt)
      do 30 itr=1,nt
      ish=indin(itr)
      if(ish.gt.ii) go to 2000
   30 mi(itr)=ish
      lit=ktype(ii)
      mini=kmin(ii)
      maxi=kmaz(ii)
      loci=kloc(ii)-mini
      do 1000 jj=1,ii
      do 40 itr=1,ntwd
   40 indout(itr)=iso(jj,itr)
      call isoout(nt)
      do 60 itr=1,nt
      jsh=indin(itr)
      if(jsh.gt.ii) go to 1000
      ish=mi(itr)
      if(ish.ge.jsh) go to 50
      n=ish
      ish=jsh
      jsh=n
   50 if(ish.eq.ii.and.jsh.gt.jj) go to 1000
   60 continue
      ljt=ktype(jj)
      minj=kmin(jj)
      maxj=kmaz(jj)
      locj=kloc(jj)-minj
      iandj=ii.eq.jj
      jmax=maxj
c
c     ----- find the equivalent blocks
c     ----- transfer equivalent block into t-matrix
c     ----- compute (r) t (r)
c     ----- put the result back into the (i,j) block of the h-matrix
c
      do 110 itr=1,nt
      ntr=itr
      kk=mi(itr)
      ll=indin(itr)
      lock=kloc(kk)-kmin(kk)
      locl=kloc(ll)-kmin(ll)
      do 90 k=mini,maxi
      lck=lock+k
      if(iandj) jmax=k
      do 90 l=minj,jmax
      if(ll.gt.kk) go to 70
      kl=ia(lck)+locl+l
      go to 80
   70 kl=ia(locl+l)+lck
   80 t(k,l)=f(kl)
      if(iandj) t(l,k)=f(kl)
   90 continue
      call rhr
      do 100 i=mini,maxi
      lci=ia(loci+i)+locj
      if(iandj) jmax=i
      do 100 j=minj,jmax
      ij=lci+j
  100 h(ij)=h(ij)+t(i,j)
  110 continue
c
c     ----- for each block (k,l) equivalent to (i,j)
c     ----- find the transformation that maps (k,l) into (i,j)
c     ----- compute (r) t (r)
c     ----- put the result back into the (k,l) block of the h-matrix
c
      do 200 itr=1,nt
      if(itr.eq.1) go to 200
      kk=mi(itr)
      ll=indin(itr)
      if(kk.ge.ll) go to 120
      k=ll
      l=kk
      go to 130
  120 k=kk
      l=ll
  130 if(k.eq.ii.and.l.eq.jj) go to 200
      ntr=itr+1
      if(ntr.gt.nt) go to 160
      do 150 it=ntr,nt
      i=mi(it)
      j=indin(it)
      if(i.ge.j) go to 140
      ij=i
      i=j
      j=ij
  140 if(i.eq.k.and.j.eq.l) go to 200
  150 continue
  160 continue
      ntr=invt(itr)
      do 170 i=mini,maxi
      lci=ia(loci+i)+locj
      if(iandj) jmax=i
      do 170 j=minj,jmax
      t(i,j)=h(lci+j)
      if(iandj) t(j,i)=h(lci+j)
  170 continue
      call rhr
      lock=kloc(kk)-kmin(kk)
      locl=kloc(ll)-kmin(ll)
      do 190 k=mini,maxi
      lck=lock+k
      if(iandj) jmax=k
      do 190 l=minj,jmax
      if(ll.gt.kk) go to 180
      kl=ia(lck)+locl+l
      go to 190
  180 kl=ia(locl+l)+lck
  190 h(kl)=t(k,l)
  200 continue
 1000 continue
 2000 continue
      ntdum=nt
c
      tn=dble(ntdum)
      do 3000 i=1,nx
 3000 f(i)=h(i)/tn
c
c     transfer the whole fock matrix into the v -array
c
c
      return
      end
