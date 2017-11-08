      subroutine mulken                                                 
      implicit real*8 (a-h,o-z)                                         
      include 'pshf.prm'
      character*4 iflab                                                 
      logical*1 last,out                                                
      logical*1 yiet,yjet,ycou
      common/cd/r(ncouz*(doa**2)),ndd,ndt,ncf1,ncf2,metat,icf2,ncou,  
     &ibeg(ncouz),ietat(ncouz),jetat(ncouz),ycou(metz,metz),      
     &yiet(metz),yjet(metz)                                             
      common/iofile/ir,iw,ip,is,iq,ih,iv,if2,if3,if4,ioda(19)           
      common/output/nprint,itol,icut,normf,normp                        
      common/runlab/iflab(3,doa)                                      
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,zan(dc),c(3,dc)          
      common/atlim/limlow(dc),limsup(dc)                                
      dimension s(doa*(doa+1)/2),d(doa*doa),            
     1 goc(doa,2)                                                       
      dimension e(dc*(dc+1)/2),gac(dc,2)                                        
      data alpha,beta,all /8h   alpha,8h    beta,8h     all/            
 9999 format(//,10x,28(1h-),/,10x,'ci mulliken population analysis',    
     1 /,10x,28(1h-))                                                   
 9998 format(/,10x,'----- total gross population in ao"s -----',/)      
 9997 format(10x,i5,2x,3a4,f12.5)                                       
 9996 format(/,10x,'----- condensed to atoms -----',/)                  
 9995 format(/,10x,'----- total gross population on atoms ----',/)      
 9994 format(/,10x,'----- ao"s spin population -----',/)                
 9993 format(/,10x,'----- atomic spin population -----',/)              
 9992 format(10x,i5,2x,2a4,2x,f6.1,f12.5)                               
      out=nprint.eq.3                                                   
c                                                                       
c     ----- determine number of basis functions per atom -----          
c                                                                       
      call aolim                                                        
c                                                                       
c     ----- read in overlap matrix -----                                
c                                                                       
      nav=ioda(2)                                                       
      read(2)(s(i),i=1,nx)                                              
      last=.true.                                                       
c     last=.false.                                                      
c 100 continue                                                          
      do 2000 icou=1,ncou                                               
      m1=ietat(icou)                                                    
      m2=jetat(icou)                                                    
      ind=ibeg(icou)
      ij=0
      do i=1,num
        do j=1,i
          ij=ij+1
            d(ij)=r(ind+(i-1)*num+j)
        enddo     
      enddo     
      ind=ind+1
c                                                                       
c     ----- do a mulliken population analysis ----                      
c           calculate overlap population                                
c                                                                       
      call ovlpop(d,s,nx)                                          
      write(iw,9999)                                                    

c                                                                       
c     ----- calculate total gross population in ao's ----               
c                                                                       
      call grossc(d,goc(1,1),num)                                  
      write(iw,9998)                                                    
      write(iw,9997) (i,(iflab(k,i),k=1,3),goc(i,1),i=1,num)            
c                                                                       
c     ----- compress from orbitals to atoms -----                       
c                                                                       
      call atpop(d,e,nat)                                          
      write(iw,9996)                                                    
      if(out) call matout(e,nat)                                        
c                                                                       
c     ----- calculate total gross population on atoms -----             
c                                                                       
      call grossc(e,gac(1,1),nat)                                       
      write(iw,9995)                                                    
      write(iw,9992)(i,iflab(1,limlow(i)),iflab(2,limlow(i)),           
     * zan(i),gac(i,1),i=1,nat)                                         
 2000 continue                                                          
      return
      end                                                               
