      subroutine morue
      implicit real*8(a-h,o-x,z),logical*1(y)
      include 'pshf.prm'
      real*8 fac(2,6),val3d1(3),val3q1(3),val4dquintet(6),val4d(6)
      
      integer*8 nhntd
      integer*4 trou,part,nexst,nd
      integer*4 ibid
      character*20 title
      character*40 typener 
      integer*4 iorb,iwspin,itsym,its,isytr,kkk,lll    
      integer,dimension(:,:),allocatable :: deg,numtab,ordertab    
      logical,dimension(:),allocatable :: ygood
      common/dinde/trou(8*ndimh),part(8*ndimh),nexst(ndimh),
     * nd(ndimh)
      common/e0cip/ebmp(metz),eben(metz),tau
      common/spo/ feff(2*doa),fdiag(2*doa),
     * iorb(2*doa),iwspin(2*doa),ispin(2*doa),
     *itsym(2*doa),its(nsymz,nsymz),yoc(2*doa),isytr(2*doa,nsymz)      
      common/ic/c(ndetz*metz),e(ndetz),eb(metz),e2(metz),
     *emp(metz),e2en(metz),teste(metz),elevy(metz),e2mp(metz),ac,ei,
     *hdef,dpivot,fmpb(2*doa),dmpb(10),ishif,ietat
      common  /cipsi/ norb,mnorb,noca,nocb,
     *igela,igelb,ncf,ncper,nsym,isym,isz,ntrsy,metat,
     *isymat(smax),nelac,nprint,yprt,yion,ystkmp
      common/nom33/title
      common/mod90/noa,ydp      
      data zero/0.d0/
      dimension cunt(ndimh,metz)
      dimension numero(metz),hpert(metz),hefmp(nhefz),
     *hefen(nhefz),hefenb(nhefz),
     *hef(nhefz)
      integer itab(lsize,smax),jtab(lsize,smax)
      integer itabbuf(ndimh,3),jtabbuf(ndimh,3)      
      dimension ntst(8,2),npst(8,2),nest(2)      
      dimension vbuf(lsize,smax),vbufbuf(ndimh,3)
      dimension nhij(smax),nbuf(smax),nbufbuf(3)
      dimension it(smax),jt(smax),jtbuf(3)
      dimension stefmp(ndimh),stefen(ndimh)
      dimension hijsym(10,10)
c
c
      call ijkf
      call reijkl(norb,nsym,ydp)
   
      read(4) norb,noca,nocb,metat,ndims,isz,yion,ii,jj
      
c      read(4) norb,noca,nocb,metat,ndims,isz,yion,ii,jj,
c     * (ibid,ibid,i=1,ndims),
c     & (ibid,ibid,i=1,ii),(c(i),i=1,jj),(e(i),i=1,metat),
c     &  (ebmp(i),i=1,metat),(eben(i),i=1,metat),(numero(i),i=1,metat)
c      mmm=metat*(metat+1)/2
c      read(4)
c      read(4) (hef(i),i=1,mmm)

      
      do 3000 i=1,mnorb
 3000 yoc(i)=.false.
      yoc(mnorb+1)=.false.
       
      nvirtelac=nelac
      if(yion) nvirtelac=nvirtelac+1
c     deg(*,10+1) : allocation for up to 5e-
c     for 5e- there are 10 possible combinaisons of ms=-1/2    
      nsmax=smax
      allocate(deg(ncf,10+1),numtab(ncf,nsmax),ordertab(ncf,10),
     * ygood(ncf))
      ordertab=0
      ygood=.true.
      yspin=.false.
      do i=1,nsmax
        if(isymat(i).eq.1) then
          yspin=.true.
          exit
        endif  
      enddo
      if(nelac.eq.1.or.nelac.gt.4) yspin=.false.
      if(.not.yspin) then

        do i=1,ncf
          deg(i,1)=1
          deg(i,2)=i
        enddo           
      else
      
        call degen(nelac,trou,part,nd,nexst,ygood,
     *               deg,ordertab,ncf,norb)
      
        write(*,*) 'Spin symmetrization done'
      endif

c      do i=1,ncf
c      if(deg(i,1).eq.6)write(*,*) 'quint ',i,nexst(i)
c      enddo
      
c     ************test determinant*********      
c      do i=1,50
c      write(*,*) '*********************'
c      write(*,*) ygood(i)
c      write(*,'(5i7)') deg(i,1),(deg(i,j),j=2,deg(i,1)+1)
c      write(*,*) (ordertab(i,j),j=1,3)
c      write(*,'(8i7)') (part(nd(deg(i,2))+jj),jj=1,
c     *nexst(deg(i,2))),
c     *(trou(nd(deg(i,2))+jj),jj=1,nexst(deg(i,2)))        
c      if(deg(i,1).gt.1) then
c      write(*,'(8i7)') (part(nd(deg(i,3))+jj),jj=1,
c     *nexst(deg(i,3))),
c     *(trou(nd(deg(i,3))+jj),jj=1,nexst(deg(i,3)))        
c      write(*,'(8i7)') (part(nd(deg(i,4))+jj),jj=1,
c     *nexst(deg(i,4))),
c     *(trou(nd(deg(i,4))+jj),jj=1,nexst(deg(i,4)))  
c      endif
c      enddo     


c     start computing the H matrix 
      nhntd=0
      cunt=0d0
      numtab=0
      nhij=0
      nbuf=0
      nbufbuf=0
      itab=0
      itabbuf=0
      jtab=0
      jtabbuf=0
      it=0
      jt=0
      jtbuf=0
      ntot=0
      do nss=1,smax
        rewind(70+nss)
      enddo

      do i=1,ncf
        if(.not.ygood(i)) cycle
        ne1=nexst(deg(i,2))
        nest(1)=ne1
        do j=1,i
          if(.not.ygood(j)) cycle
          ne2=nexst(deg(j,2))
          ndif=ne1-ne2
          if(ndif.gt.2) then
              call wrtempt(nelac,i,it,jt,jtbuf,deg(i,1),deg(j,1),
     *       ncf,numtab,isymat) 
            cycle
          endif
          nest(2)=ne2
c         non-zero couplings, compute the different hij
          hijsym=0d0

          do ii=1,deg(i,1)
            nind1=deg(i,ii+1)
            ndd1=nd(nind1)
            do k=1,ne1
              ntk=trou(ndd1+k)
              ntst(k,1)=ntk
              yoc(ntk)=.true.
              npk=part(ndd1+k)
              npst(k,1)=npk
              yoc(npk)=.true.
            end do
            do jj=1,deg(j,1)
              nind2=deg(j,jj+1)
              ndd2=nd(nind2)
c             diagonal element              
              if(nind1.eq.nind2) then
                stfn=hdig(nest,ntst,npst)
                hijsym(ii,jj)=stfn
c                if(nind1.gt.ndims) then
c                  stefen(nind1)=stfn
c                  stefmp(nind1)=hmp(nest,ntst,npst)                
c                endif
                cycle
              endif
              
              ndif2=ndif
              do k=1,ne2
                ntk=trou(ndd2+k)
                ntst(k,2)=ntk
                if(.not.yoc(ntk))ndif2=ndif2+1
                npk=part(ndd2+k)
                npst(k,2)=npk
                if(.not.yoc(npk))ndif2=ndif2+1
              end do
              if(ndif2.gt.2) cycle 
              hijsym(ii,jj)=hntd(nest,ntst,npst) 
              nhntd=nhntd+1

c      if element couple space S and M : add to 2nd order pert. calc 

c              if((nind2.gt.ndims.and.nind1.le.ndims).or.
c     *         (nind1.gt.ndims.and.nind2.le.ndims)) then
c                mmind=min(nind1,nind2)
c                mpind=max(nind1,nind2)
c                mi=mmind-ndims
c                do m=1,metat
c                  mi=mi+ndims
c                  cunt(mpind,m)=cunt(mpind,m)+c(mi)*hijsym(ii,jj)
c                enddo
c              endif
            enddo
c           reinitialize yoc 
            do k=1,ne1
              ntk=trou(ndd1+k)
              yoc(ntk)=.false.
              npk=part(ndd1+k)
              yoc(npk)=.false.
            end do
          enddo

          if(yspin) then
          
c          if(i.eq.112.and.j.eq.112) then
c          write(*,'(6i7)')(ordertab(i,ii),ii=1,6)
c          write(*,'(6i7)') (deg(i,jj+1),jj=1,deg(i,1))
c          write(*,'(8i7)') (part(nd(deg(j,2))+jj),jj=1,
c     *     nexst(deg(j,2))),
c     *     (trou(nd(deg(j,2))+jj),jj=1,nexst(deg(j,2)))
c          write(*,'(8i7)') (part(nd(deg(j,3))+jj),jj=1,
c     *     nexst(deg(j,3))),
c     *     (trou(nd(deg(j,3))+jj),jj=1,nexst(deg(j,3)))
c          write(*,'(8i7)') (part(nd(deg(j,4))+jj),jj=1,
c     *     nexst(deg(j,4))),
c     *     (trou(nd(deg(j,4))+jj),jj=1,nexst(deg(j,4)))      
c          write(*,'(8i7)') (part(nd(deg(j,5))+jj),jj=1,
c     *     nexst(deg(j,5))),
c     *     (trou(nd(deg(j,5))+jj),jj=1,nexst(deg(j,5)))   
c          write(*,'(8i7)') (part(nd(deg(j,6))+jj),jj=1,
c     *     nexst(deg(j,6))),
c     *     (trou(nd(deg(j,6))+jj),jj=1,nexst(deg(j,6)))   
c          write(*,'(8i7)') (part(nd(deg(j,7))+jj),jj=1,
c     *     nexst(deg(j,7))),
c     *     (trou(nd(deg(j,7))+jj),jj=1,nexst(deg(j,7)))        
c          write(*,'(6f20.12)') (hijsym(1,ii),ii=1,6)
c          write(*,'(6f20.12)') (hijsym(2,ii),ii=1,6)
c          write(*,'(6f20.12)') (hijsym(3,ii),ii=1,6)
c          write(*,'(6f20.12)') (hijsym(4,ii),ii=1,6)
c          write(*,'(6f20.12)') (hijsym(5,ii),ii=1,6)
c          write(*,'(6f20.12)') (hijsym(6,ii),ii=1,6)
          
          
c          endif
          
c           check couplings between quintet and singlet           
          if(deg(i,1).eq.1.and.deg(j,1).eq.6) then
          hij=0d0
          val4d=(/1d0,1d0,0d0,0d0,0d0,0d0/)
          val4d=val4d/sqrt(2d0)
          fac(1,1)=val4d(ordertab(i,1))
    !      fac(1,2)=val4d(ordertab(i,2))   
          val4dquintet=(/0d0,-1d0,-1d0,1d0,1d0,0d0/)
          val4dquintet=val4dquintet/sqrt(6d0)     
          fac(2,1)=val4dquintet(ordertab(j,1))
          fac(2,2)=val4dquintet(ordertab(j,2)) 
          fac(2,3)=val4dquintet(ordertab(j,3))
          fac(2,4)=val4dquintet(ordertab(j,4))  
          fac(2,5)=val4dquintet(ordertab(j,5))
          fac(2,6)=val4dquintet(ordertab(j,6))            
           hij=0d0        
          do ii=1,deg(i,1)
            do jj=1,deg(j,1)
              hij = hij + fac(1,ii)*fac(2,jj)* hijsym(ii,jj)
            enddo
          enddo  
          if(dabs(hij).gt.10d-10)then 
          write(*,*) 'ERROR  ' ,i,j
          write(*,*) hij
          write(*,*) ne1,ne2
          write(*,'(6i7)') (deg(i,jj+1),jj=1,deg(i,1))
          write(*,'(6i7)') (deg(j,jj+1),jj=1,deg(j,1))
          write(*,'(6i7)') (ordertab(i,jj),jj=1,deg(i,1)) 
          write(*,'(6i7)') (ordertab(j,jj),jj=1,deg(j,1)) 
          write(*,'(8i7)') (part(nd(deg(j,2))+jj),jj=1,
     *     nexst(deg(j,2))),
     *     (trou(nd(deg(j,2))+jj),jj=1,nexst(deg(j,2)))
          write(*,'(8i7)') (part(nd(deg(j,3))+jj),jj=1,
     *     nexst(deg(j,3))),
     *     (trou(nd(deg(j,3))+jj),jj=1,nexst(deg(j,3)))
          write(*,'(8i7)') (part(nd(deg(j,4))+jj),jj=1,
     *     nexst(deg(j,4))),
     *     (trou(nd(deg(j,4))+jj),jj=1,nexst(deg(j,4)))  
          write(*,*)hijsym(1,1),hijsym(1,2),hijsym(1,3)
          write(*,*)hijsym(1,4),hijsym(1,5),hijsym(1,6)         
          endif
          endif
            call wrthij(nbuf,itab,jtab,vbuf,nbufbuf,itabbuf,jtabbuf,
     *       vbufbuf,nhij,nelac,ncf,i,j,it,jt,jtbuf,deg(i,1),deg(j,1),
     *       hijsym,ordertab,numtab,isymat)
          else
            call wrthijtot(nbuf(1),itab,jtab,vbuf,nhij(1),ncf,
     *     i,j,hijsym(1,1),numtab)
          endif
        enddo
      enddo
      
      write(*,*) '************'
      write(*,*) ' call to hntd, nb of element written on disk : '
      write(*,'(i12,4x,6i12)') nhntd,(nhij(i),i=1,6)
      
c     write last buffer  
      do nss=1,smax
        itab(lsize,nss)=-1
        jtab(lsize,nss)=nbuf(nss)
        if(nhij(nss).eq.0) cycle
        write(70+nss) itab(:,nss),jtab(:,nss),vbuf(:,nss)
      enddo       
      if(.not.(yspin)) it(1)=ncf
        write(6,846) 
  846   format(//' *** fin de la construction de h ***')

c      write information on the mat file
       rewind(70)
       write(*,*) 'Total number of determinant : ',ncf
       write(*,*) 'Number of determinant in each 2S+1 matrix :'
       write(*,'(6i9)') (it(i),i=1,smax)      
       write(70) ncf,nelac,nsmax,yspin
       write(70) (it(i),i=1,nsmax),((deg(i,j),i=1,ncf),j=1,11),
     *  ((ordertab(i,j),i=1,ncf),j=1,10),
     *  ((numtab(i,j),i=1,ncf),j=1,nsmax)
      
      deallocate(deg,numtab,ordertab,ygood)
      
      if(.true.) go to 8300
       
c      write(6,3456)
c      write(6,3457) (stefen(i),i=1,ncf)
 3456 format(//' elements diagonaux')
 3457 format(12f10.6)
      if(ystkmp) write(63) (stefmp(i),i=1,ncf)
      if(ystkmp) write(6,2378)
 2378 format(/, ' elements diagonaux m.p. sur file 63')
 
 
c      Calcul des perturbations 
 
      write(6,9220)
9220  format(//)
      write(6,9344)
 9344 format(/' caracteristique des etats a l''ordre zero')
      write(6,9230)
9230  format(1x,40('*'),/)
      write(6,9233)
9233  format('  etat',4x,'partition m.p.b',10x,'partition e.n.vp',10x,
     *'partition e.n.b',/)
      do 9345 i=1,metat
 9345 write(6,9346) i,ebmp(i),e(i),eben(i)
 9346 format(i4,f18.8,2f25.8)
      ij=0
      do 25 i=1,metat
      do 25 j=1,i
      ij=ij+1
       hefenb(ij)=0.d0
       hefmp(ij)=0.d0
25     hefen(ij)=0.d0
       
      if(ndims.eq.ncf) go to 8300

      do 2 k=ndims+1,ncf
      stenk=stefen(k)
      stmpk=stefmp(k)
      enk2=stenk+stenk
      empk2=stmpk+stmpk
      ii=0
      do 103 m=1,metat
      do 103 n=1,m
      ii=ii+1
      denb=eben(m)+eben(n)-enk2
      dmp=ebmp(m)+ebmp(n)-empk2
      den=e(m)+e(n)-enk2
      hmn=cunt(k,m)*cunt(k,n)
      hmn=hmn+hmn
      hefmp(ii)=hefmp(ii)+hmn/dmp
      hefen(ii)=hefen(ii)+hmn/den
      hefenb(ii)=hefenb(ii)+hmn/denb
103   continue
2     continue

      write(6,9000)
9000  format(//,1x,'contribution 2n ordre des determinants moyens',/)
      write(6,9234)
9234  format(1x,45('*'))
      write(6,9001)
9001  format(1x,'perturbation mp',/)
      do 120 m=1,metat
      i1=m*(m-1)/2+1
      i2=m*(m+1)/2
120   write(6,9004) (hefmp(ii),ii=i1,i2)
      write(6,9220)
      write(6,9003)
9003  format(1x,'perturbation envp',/)
      do 140 m=1,metat
      i1=m*(m-1)/2+1
      i2=m*(m+1)/2
140   write(6,9004) (hefen(ii),ii=i1,i2)
      write(6,9220)
      write(6,9002)
9002  format(1x,'perturbation enb',/)
      do 130 m=1,metat
      i1=m*(m-1)/2+1
      i2=m*(m+1)/2
130   write(6,9004) (hefenb(ii),ii=i1,i2)
9004  format(11f11.6)
      write(6,9220)
      do m=1,metat
	hpert(m)=hefmp(m*(m+1)/2)
      end do 
      typener='moyens mbp '//title
      call stkener(typener,hpert,metat)
      do m=1,metat
	hpert(m)=hefen(m*(m+1)/2)
      end do 
      typener='moyens envp'//title
      call stkener(typener,hpert,metat)
      do m=1,metat
	hpert(m)=hefenb(m*(m+1)/2)
      end do 
      typener='moyens enb '//title
      call stkener(typener,hpert,metat)
      write(6,9005)
      write(6,9005)
9005  format(/,2x,128('*'))
 8300 rewind 4
      read(4)
      read(4)
      read(4)
      write(4) (hefmp(i),i=1,mmm),
     *         (hefen(i),i=1,mmm),
     *         (hefenb(i),i=1,mmm)
      return
      end
