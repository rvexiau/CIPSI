      subroutine final(index)
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      common/output/nprint,itol,icut,normf,normp
      common/restar/timlim,irest,nrec,intloc,ist,jst,kst,lst
      common/buf/tei(lsize),ila(lsize)
      common/shlt/tol,cutoff,icount,out
      common/iofile/ir,iw,ip,is
      common/nshel/ex(dgp),cs(dgp),cp(dgp),cd(dgp),cf(dgp),cg(dgp),
     1kstart(ds),katom(ds),ktype(ds),kng(ds),kloc(ds),kmin(ds),kmax(ds),
     2nshell
      common/shlnos/qq4,ii,jj,kk,ll,lit,ljt,lkt,llt,
     1 mini,minj,mink,minl,maxi,maxj,maxk,maxl,
     1 loci,locj,lock,locl,ij,kl,ijkl,nij
      common/misc/maxll,iandj,kandl,same
      if(index.eq.1) go to 400
      irest=2
      ist=ii
      jst=jj
      kst=kk
      lst=ll+1
      if(lst.le.maxll) go to 200
      lst=1
      kst=kk+1
      if(kst.le.ii) go to 200
      kst=1
      jst=jj+1
      if(jst.le.ii) go to 200
      jst=1
      ist=ii+1
      if(ist.gt.nshell) go to 400
  200 write(iw,1001) timlim,nprint,itol,icut,normf,normp,irest,
     1 ist,jst,kst,lst,nrec,icount
      write(ip,1001) timlim,nprint,itol,icut,normf,normp,irest,
     1 ist,jst,kst,lst,nrec,icount
 1001 format(f10.0,10i3,i10,i5)
      write(is) tei,ila
      write(iw,1000) nrec
 1000 format(' there is (are) ',i7,' record(s) of 2e-integrals',
     1 ' written on the integral file (is) and ',i5,' integrals in'
     2' the last buffer')
      write(iw,'(a27)')'nombre de blocs sur la file'
      write(iw,*) is,nrec
      write(iw,1002)
 1002 format(' warning   .......   warning   .......',/,
     1 20x,'this job must be restarted ..... ')
      return
  400 continue
      ila(icount)=0
      write(is) tei,ila
      irest=3
      intloc=icount
      write(iw,1000) nrec,icount-1
      write(iw,'(a27)')'nombre de blocs sur la file'
      write(iw,*) is,nrec
      write(iw,1003)
 1003 format(' ...... end of two-electron integrals ......')
      return
      end
