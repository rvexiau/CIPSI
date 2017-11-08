      subroutine itijkl_stan(bijkls,bijkld,ydp)
      implicit real*8 (a-h,o-x,z),logical*1 (y)                         
      include 'pshf.prm'                                  
      real*4 bufs(lsize)
      dimension bufd(lsize)
      real*4 bijkls(*)
      dimension bijkld(*)                                         
c                                                                       
      common/truc/indic(doa*(doa+1)/2),jndic(doa*(doa+1)/2),
     & kt(10000),nbo(99),ndeb(500),nad(10000),
     & lndic(doa*(doa+1)/2)                                   
      integer*2 indic,jndic,nad,kt                                      
      common/sy/itsym(doa),its(nsymz,nsymz),mdegen(5),isydeg(nsymz,nsy
     *mz),itsyp(doa)                                                  
      common nos(nsymz),ns(doa),nus(nsymz),nbp(99),npair,ipair(2000),  
     y noms(doa,nsymz)                                                
      integer*2 noms                                                    
      num(i)=i*(i-1)/2
      ijkl=0                                                            
      mijkl=0                                                           
      mbuf=lsize                                                        
      do 100 n=1,npair                                                  
      ij=ipair(n+n-1)                                                   
      nij=nbp(ij)                                                       
      kl=ipair(n+n)                                                     
      nkl=nbp(kl)                                                       
c      write(2,*)'n=',n,' kl=',kl,'  ij=',ij
      if(ij.lt.kl) go to 2                                              
      ijkls=num(ij)+kl                                                  
      go to 3                                                           
    2 ijkls=num(kl)+ij                      
    3 ysam=kt(ijkls).eq.1.or.kt(ijkls).eq.10                            
c      write(2,*)'ijkls=',ijkls                                    
      do 15 m=1,20                                                      
      if(num(m).lt.ij) go to 15                                         
      is=m-1                                                            
      ni=nos(is)                                                        
      js=ij-num(is)                                                     
      nj=nos(js)                                                        
      go to 16                                                          
15    continue                                                          
16    do 17 m=1,20                                                      
      if(num(m).lt.kl) go to 17                                         
      ks=m-1                                                            
      nk=nos(ks)                                                        
      ls=kl-num(ks)                                                     
      nl=nos(ls)                                                        
      go to 20                                                          
17    continue                                                          
20    continue                                                          
c                                                                       
      ml=nl                                                             
      mj=nj                                                             
      nij=0                                                             
c      write(2,*) n,ij,kl,ijkls,kt(ijkls),ni,nj,nk,nl,mbuf              
      do 100 i=1,ni                                                     
c     ii=noms(i,is)                                                     
      if(is.eq.js) mj=i                                                 
      do 100 j=1,mj                                                     
c      jj=noms(j,js)                                                    
      nij=nij+1                                                         
      nkl=0                                                             
      do 90 k=1,nk                                                      
c      kk=noms(k,ks)                                                    
      if(ks.eq.ls) ml=k                                                 
       do 90 l=1,ml                                                     
c      ll=noms(l,ls)                                                    
      nkl=nkl+1                                                         
      if(nkl.gt.nij.and.ysam) go to 100                                 
      mbuf=mbuf+1                                                       
      if(mbuf.le.lsize) go to 40                                        
      if(.not.ydp) read(50) bufs                                        
      if (ydp) read(50) bufd                                            
      mbuf=1                                                            
40    continue                                                          
c     if(yocs(i).or.yocs(j).or.yocs(k).or.yocs(l)) go to 70             
c     if(yocs(i+n).or.yocs(j.+n).or.yocs(k+n).or.yocs(l+n) ) go to 70   
c     if(i.le.noc1.and.j.le.noc1) go to 90                              
c     if(i.gt.noc2.and.j.gt.noc2) go to 90                              
c     if(k.le.noc1.and.l.le.noc1) go to 90                              
c     if(k.gt.noc2.and.l.gt.noc2) go to 90                              
70    ijkl=ijkl+1                                                       
      if(.not.ydp ) bijkls(ijkl)=bufs(mbuf)                             
      if( ydp ) bijkld(ijkl)=bufd(mbuf)                                 
90    continue                                                          
100   continue                                                          
600   format(1x,15f8.4)                                                 
       return                                                           
        end                                                             
