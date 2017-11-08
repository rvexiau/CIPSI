      subroutine shells(index)
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      logical*1 iandj,kandl,same
      common/output/nprint,itol,icut,normf,normp
      common/infoa/nat,ich,mul,numscf,nnp,ne,na,nb,zan(dc),c(3,dc)
      common/nshel/ex(dgp),cs(dgp),cp(dgp),cd(dgp),cf(dgp),cg(dgp),
     1kstart(ds),katom(ds),ktype(ds),kng(ds),kloc(ds),kmin(ds),kmax(ds),
     2nshell
      common/root/xx,u(9),w(9),nroots
      common/shlnos/qq4,ii,jj,kk,ll,lit,ljt,lkt,llt,
     1 mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     1 loci,locj,lock,locl,ij,kl,ijkl,nij
     1,ni,nj,nk,nl
      common/index/ijx(225),ijy(225),ijz(225),ik(225)
     1            ,klx(225),kly(225),klz(225)
      common/shlinf/ga(20),csa(20),cpa(20),cda(20),cfa(20),cga(20),
     1              gb(20),csb(20),cpb(20),cdb(20),cfb(20),cgb(20),
     1              gc(20),csc(20),cpc(20),cdc(20),cfc(20),cgc(20),
     1              gd(20),csd(20),cpd(20),cdd(20),cfd(20),cgd(20),
     1              ax,ay,az,bx,by,bz,rab,cx,cy,cz,dx,dy,dz,rcd,
     1              nga,ngb,ngc,ngd
      common/misc/maxll,iandj,kandl,same
      dimension ix(35),iy(35),iz(35),jx(35),jy(35),jz(35),
     1 kx(35),ky(35),kz(35),lx(35),ly(35),lz(35)
      data lx /   0,  1,  0,  0,  2,  0,  0,  1,  1,  0,
     1            3,  0,  0,  2,  2,  1,  0,  1,  0,  1,
     2            4,  0,  0,  3,  3,  1,  0,  1,  0,  2,
     3            2,  0,  2,  1,  1/
      data kx /   0,  5,  0,  0, 10,  0,  0,  5,  5,  0,
     1           15,  0,  0, 10, 10,  5,  0,  5,  0,  5,
     2           20,  0,  0, 15, 15,  5,  0,  5,  0, 10,
     3           10,  0, 10,  5,  5/
      data jx /   0, 25,  0,  0, 50,  0,  0, 25, 25,  0,
     1           75,  0,  0, 50, 50, 25,  0, 25,  0, 25,
     2          100,  0,  0, 75, 75, 25,  0, 25,  0, 50,
     3           50,  0, 50, 25, 25/
      data ix /   1,126,  1,  1,251,  1,  1,126,126,  1,
     1          376,  1,  1,251,251,126,  1,126,  1,126,
     2          501,  1,  1,376,376,126,  1,126,  1,251,
     3          251,  1,251,126,126/
      data ly /   0,  0,  1,  0,  0,  2,  0,  1,  0,  1,
     1            0,  3,  0,  1,  0,  2,  2,  0,  1,  1,
     2            0,  4,  0,  1,  0,  3,  3,  0,  1,  2,
     3            0,  2,  1,  2,  1/
      data ky /   0,  0,  5,  0,  0, 10,  0,  5,  0,  5,
     1            0, 15,  0,  5,  0, 10, 10,  0,  5,  5,
     2            0, 20,  0,  5,  0, 15, 15,  0,  5, 10,
     3            0, 10,  5, 10,  5/
      data jy /   0,  0, 25,  0,  0, 50,  0, 25,  0, 25,
     1            0, 75,  0, 25,  0, 50, 50,  0, 25, 25,
     2            0,100,  0, 25,  0, 75, 75,  0, 25, 50,
     3            0, 50, 25, 50, 25/
      data iy /   1,  1,126,  1,  1,251,  1,126,  1,126,
     1            1,376,  1,126,  1,251,251,  1,126,126,
     2            1,501,  1,126,  1,376,376,  1,126,251,
     3            1,251,126,251,126/
      data lz /   0,  0,  0,  1,  0,  0,  2,  0,  1,  1,
     1            0,  0,  3,  0,  1,  0,  1,  2,  2,  1,
     2            0,  0,  4,  0,  1,  0,  1,  3,  3,  0,
     3            2,  2,  1,  1,  2/
      data kz /   0,  0,  0,  5,  0,  0, 10,  0,  5,  5,
     1            0,  0, 15,  0,  5,  0,  5, 10, 10,  5,
     2            0,  0, 20,  0,  5,  0,  5, 15, 15,  0,
     3           10, 10,  5,  5, 10/
      data jz /   0,  0,  0, 25,  0,  0, 50,  0, 25, 25,
     1            0,  0, 75,  0, 25,  0, 25, 50, 50, 25,
     2            0,  0,100,  0, 25,  0, 25, 75, 75,  0,
     3           50, 50, 25, 25, 50/
      data iz /   1,  1,  1,126,  1,  1,251,  1,126,126,
     1            1,  1,376,  1,126,  1,126,251,251,126,
     2            1,  1,501,  1,126,  1,126,376,376,  1,
     3          251,251,126,126,251/
      if(index.eq.2) go to 200
      iandj=ii.eq.jj
c
c     ----- permute ii and jj shells, for their type
c
      if(ktype(ii).ge.ktype(jj)) go to 5
      in=jj
      jn=ii
      go to 10
    5 in=ii
      jn=jj
   10 continue
c
c     ----- ishell
c
      i=katom(in)
      ax=c(1,i)
      ay=c(2,i)
      az=c(3,i)
      i1=kstart(in)
      i2=i1+kng(in)-1
      lit=ktype(in)
      ni=lit
      if(ni.ge.6) ni=ni-5
      mini=kmin(in)
      maxi=kmax(in)
      loci=kloc(in)-mini
      nga=0
      do 50 i=i1,i2
      nga=nga+1
      ga(nga)=ex(i)
      csa(nga)=cs(i)
      cpa(nga)=cp(i)
      cda(nga)=cd(i)
      cfa(nga)=cf(i)
      cga(nga)=cg(i)
   50 continue
c
c     ----- jshell
c
      j=katom(jn)
      bx=c(1,j)
      by=c(2,j)
      bz=c(3,j)
      j1=kstart(jn)
      j2=j1+kng(jn)-1
      ljt=ktype(jn)
      nj=ljt
      if(nj.ge.6) nj=nj-5
      minj=kmin(jn)
      maxj=kmax(jn)
      locj=kloc(jn)-minj
      ngb=0
      do 100 j=j1,j2
      ngb=ngb+1
      gb(ngb)=ex(j)
      csb(ngb)=cs(j)
      cpb(ngb)=cp(j)
      cdb(ngb)=cd(j)
      cfb(ngb)=cf(j)
      cgb(ngb)=cg(j)
  100 continue
      rab=((ax-bx)**2+(ay-by)**2+(az-bz)**2)
c
c     ----- prepare indices for pairs of (i,j) functions
c
      ij=0
      max=maxj
      do 150 i=mini,maxi
      nx=ix(i)
      ny=iy(i)
      nz=iz(i)
      if(iandj) max=i
      do 150 j=minj,max
      ij=ij+1
      ijx(ij)=nx+jx(j)
      ijy(ij)=ny+jy(j)
      ijz(ij)=nz+jz(j)
  150 continue
      return
  200 continue
      kandl=kk.eq.ll
      same=ii.eq.kk.and.jj.eq.ll
c
c     ----- permute kk and ll shells, for their type
c
      if(ktype(kk).ge.ktype(ll)) go to 205
      kn=ll
      ln=kk
      go to 210
  205 kn=kk
      ln=ll
  210 continue
c
c     ----- kshell
c
      k=katom(kn)
      cx=c(1,k)
      cy=c(2,k)
      cz=c(3,k)
      k1=kstart(kn)
      k2=k1+kng(kn)-1
      lkt=ktype(kn)
      nk=lkt
      if(nk.ge.6) nk=nk-5
      mink=kmin(kn)
      maxk=kmax(kn)
      lock=kloc(kn)-mink
      ngc=0
      do 250 k=k1,k2
      ngc=ngc+1
      gc(ngc)=ex(k)
      csc(ngc)=cs(k)
      cpc(ngc)=cp(k)
      cdc(ngc)=cd(k)
      cfc(ngc)=cf(k)
      cgc(ngc)=cg(k)
  250 continue
c
c     ----- lsheln
c
      l=katom(ln)
      dx=c(1,l)
      dy=c(2,l)
      dz=c(3,l)
      l1=kstart(ln)
      l2=l1+kng(ln)-1
      llt=ktype(ln)
      nl=llt
      if(nl.ge.6) nl=nl-5
      minl=kmin(ln)
      maxl=kmax(ln)
      locl=kloc(ln)-minl
      ngd=0
      do 300 l=l1,l2
      ngd=ngd+1
      gd(ngd)=ex(l)
      csd(ngd)=cs(l)
      cpd(ngd)=cp(l)
      cdd(ngd)=cd(l)
      cfd(ngd)=cf(l)
      cgd(ngd)=cg(l)
  300 continue
      nroots=(ni+nj+nk+nl-2)/2
      rcd=((cx-dx)**2+(cy-dy)**2+(cz-dz)**2)
c
c     ----- prepare indices for pairs of (k,l) functions
c
      kl=0
      max=maxl
      do 350 k=mink,maxk
      nx=kx(k)
      ny=ky(k)
      nz=kz(k)
      if(kandl) max=k
      do 350 l=minl,max
      kl=kl+1
      klx(kl)=nx+lx(l)
      kly(kl)=ny+ly(l)
      klz(kl)=nz+lz(l)
  350 continue
      max=kl
      do 400 i=1,ij
      if(same) max=i
  400 ik(i)=max
      ijkl=ij*kl
      if(same) ijkl=ij*(ij+1)/2
      return
      end
