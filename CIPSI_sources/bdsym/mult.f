      subroutine mult(fmul,ncf,c)                                       
      implicit real*8(a-h,o-x,z),logical*1 (y)                          
      include 'bd.prm'
      integer*4 ne,nd,trou,part,iorb,iwspin,itsym,its,isytr,nesp,nespo  
      integer*4 imul                                                    
      character*4 ibla2                                                 
      character*8 deg          
      logical*1 qion
      dimension c(*)
      common/spob/isz,noca,nvr,qion
      common/spo/ne(ndetz+1),nd(ndetz+1),
     1trou(8*ndetz),part(8*ndetz),
     2iorb(2*doa),ispin(2*doa),yocd(2*doa,ndetz)              
      dimension fmul(metz)                      
      nborb=2*noca                                                      
      szmult=dfloat(isz)                                                
      if (qion) then                                                    
         nborb=nborb-1                                                  
         szmult=dfloat(isz)-0.5                                         
      endif                                      
c      write(6,*) ' isz ',isz,' szmult',szmult                          
      ng=1                                                              
      do ijk=1,nvr                                                    
         fmul(ijk)=0.                                                   
      enddo                                                             
      do 350 i=1,ncf                                                    
      nek=ne(i)                                                         
      nosh=2*nek                                                        
      if(qion) nosh=nosh-1                                              
      do ii=1,nek                                                       
      do jj=1,ii-1                                                      
         if (iorb(trou(nd(i)+ii)).eq.iorb(trou(nd(i)+jj)))              
     &      nosh=nosh-2                                                 
         if (iorb(part(nd(i)+ii)).eq.iorb(part(nd(i)+jj)))              
     &      nosh=nosh-2                                                 
      enddo                                                             
      enddo                                                             
      fact=(szmult*szmult)+dfloat(nosh)/2                               
c      write(6,8172) i,nek,nosh,fact                                    
c8172  format(1x,'i,nek,nosh,fact',i3,i3,i3,d15.6)                      
      do mm=1,nvr                                                     
         inda=(mm-1)*ncf+i                                              
         fmul(mm)=fmul(mm)+fact*c(inda)*c(inda)                         
      enddo                                                             
      if(i.eq.ncf) go to 350                                            
      ip=i+1                                                            
      do 304 j=ip,ncf                                                   
      ndj=nd(j)+1                                                       
      nej=ne(j)                                                         
      ndif=ne(i)-nej                                                    
      nej=nej+ndj-1                                                     
      do 300 k=ndj,nej                                                  
         ntp=trou(k)                                                    
         if(.not.yocd(ntp,i))ndif=ndif+1                                
         ntp=part(k)                                                    
         if(.not.yocd(ntp,i))ndif=ndif+1                                
300   continue                                                          
      if(ndif.le.2) then                                                
         do mm=1,nvr                                                  
            inda=(mm-1)*ncf+i                                           
            indb=(mm-1)*ncf+j                                           
            fmul(mm)=fmul(mm)+2.*sntd(i,j)*c(inda)*c(indb)              
c            write(6,*) 'mm,i,j,inda,indb',mm,i,j,inda,indb             
c            write(6,*) 'sntd,ca,cb',sntd(i,j),c(inda),c(indb)          
         enddo                                                          
      endif                                                             
c      endif                                                            
304   continue                                                          
310   continue                                                          
350   continue                                                          
      return                                                            
      end                                                               
