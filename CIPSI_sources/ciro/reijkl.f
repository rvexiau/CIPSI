      subroutine reijkl                                                 
      implicit real*8 (a-h,o-x,z),logical*1 (y)                         
      include 'pshf.prm'                
      integer*4 ne,nd,iorb,ispin,itsym,its                    
      integer*4 trou,part                  
      common/ic/e(metz2),emp(metz2),fdiag(2*doa),                     
     & fmpb(2*doa),trou(8*ndimh),part(8*ndimh),ne(2*ndimh),nd(2*ndimh), 
     & ispin(2*doa),                                                  
     & iorb(2*doa),itsym(2*doa),its(nsymz,nsymz),igels(2),icount,   
     & norb,noca,nocb,mnorb,ncf,nsym,isym,ntrsy,metat1,metat2,          
     & yoc(2*doa),yprt
      common/gel/nao                                                    
      dimension isydeg(nsymz,nsymz),itsyp(2*doa)                        
      num(i)=i*(i-1)/2                                                  
      rewind 40                                                         
      read(40) nsym,norb,noc,ntrans,(itsyp(i),i=1,norb),nao,ydp         
      write(6,*) 'dans reijkl 40 ns,nrb,nc,ntr,nao', nsym,norb,
     &noc,ntrans,nao
      read(40) ((isydeg(i,j),i=1,nsym),j=1,nsym)                        
       do 7720 i=1,nsym                                                 
        do 7720 j=1,nsym                                                
7720   its(i,j)=isydeg(i,j)                                             
       do 7730 i=1,norb                                                 
7730   itsym(i)=itsyp(i)                                                
      if(ntrans.ne.0) then                                              
      do 55 i=1,ntrans                                                  
         read(40) j,((isydeg(k,l),k=1,ntyp),l=1,j)                      
   55 continue                                                          
      end if                                                            
      return                                                            
      end                                                               
