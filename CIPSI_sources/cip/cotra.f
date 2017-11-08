      subroutine cotra(k,nf)                                            
      implicit real*8 (a-h,o-x,z),logical*1 (y)                         
      include 'pshf.prm'
      common/spo/fdiag(2*doa),hdiag(ndetz),tocom(2*doa),coeff,      
     1ndef,ispin(2*doa),iwspin(2*doa),numero(metz),          
     2ne(ndetz+1),nd(ndetz+1),trou(ntpz),part(ntpz),iorb(2*doa),      
     3itsym(2*doa),its(nsymz,nsymz),yoc(2*doa),yocs(2*doa),       
     4ygen(ndetz),yocd(2*doa,ndetz),yspin,yef(ndetz,nsymz)                   
      integer*4 nt(10),np(10)                                           
      integer*4 ne,nd,trou,part,iorb,iwspin,itsym,its  
c      contracte les excitations de nec a nec-1 ou nec-2 quand le proces
c     k donne i implique une ou deux desexcitations                     
      nu=ne(k)-1                                                        
      nup=ne(k)                                                         
      ndf=nd(nf)                                                        
      do 336 l=1,nup                                                    
      nt(l)=trou(ndf+l)                                                 
336   np(l)=part(ndf+l)                                                 
      jt=0                                                              
      jp=0                                                              
      do 340 l=1,nup                                                    
      if(nt(l).eq.0) go to 337                                          
      jt=jt+1                                                           
      trou(ndf+jt)=nt(l)                                                
337   if(np(l).eq.0) go to 340                                          
      jp=jp+1                                                           
      part(ndf+jp)=np(l)                                                
340   continue                                                          
350   continue                                                          
      call pert1(k,nf)                                                  
      return                                                            
      end                                                               
