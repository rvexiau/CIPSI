      subroutine ijkf                                                   
      implicit real*8 (a-h,o-x,z),logical*1 (y)                         
      include 'pshf.prm'
      common/jkf/ f(doa*(doa+1)/2),aj(doa*(doa+1)/2),           
     2ak(doa*(doa+1)/2)                                             
      dimension ita(20),itb(20),nv1(2),nv2(2),natu1(2),natu2(2)         
      equivalence(nv11,nv1(1)),(nv12,nv1(2)),(nv21,nv2(1)),(nv22,nv2(2))
      common/spo/fdiag(2*doa),hdiag(ndetz),tocom(2*doa),coeff,      
     1ndef,ispin(2*doa),iwspin(2*doa),numero(metz),          
     2ne(ndetz+1),nd(ndetz+1),trou(ntpz),part(ntpz),iorb(2*doa),      
     3itsym(2*doa),its(nsymz,nsymz),yoc(2*doa),yocs(2*doa),       
     4ygen(ndetz),yocd(2*doa,ndetz),yspin,yef(ndetz,nsymz)                   
      common/degene/nesp(ndetz),nespo(ndetz),isytr(ndetz,nsymz),kkbrd,  
     1ysaut(ndetz)                                                      
      integer*4 ne,nd,trou,part,iorb,iwspin,itsym,its,isytr,nesp,nespo  
      common/cipsi/escf,test,teseff,tdet,twdet,ndiag,norb,mnorb,
     1             noca,nocb,noc1,noc2,igela,igelb,ncf,nsym,isym,isz,
     2             ntrsy,metat,lecdet,maxdet,yprt,yion,ypertu,ybrd,yen,
     3             ywvi,ywvf,ytoul,ystkic,ywdeg
      common/ic/c(metz*ndetz),e(ndetz),eb(metz),e2b(metz),emp(metz),    
     1teste(metz),e2mp(metz),dpivot,ietat,ishif,e2vp(metz),             
     2fmpb(2*doa),hefmp(nhefz),hefvp(nhefz),hefb(nhefz),              
     3psimp(metz),psivp(metz)                                           
      logical iden,iprin1                                               
      num(i)=i*(i-1)/2                                                  
      nter=norb*(norb+1)/2      
      thr=1.d-10
      read(40)                                                          
      read(40) (aj(l),l=1,nter)                                         
c      write(6,*) (aj(l),l=1,nter)                                       
      read(40) (ak(l),l=1,nter)                                         
c      write(6,*) (ak(l),l=1,nter)                                       
      read(40) (f(l),l=1,nter)                                          
      do 3240 i=1,norb                                                  
      do 3240 j=1,i                                                     
      if(dabs(f(i*(i-1)/2+j)).lt.thr)f(i*(i-1)/2+j)=0.d0                
 3240 continue                                                          
c     write(6,*) (f (l),l=1,nter)                                       
      l=0                                                               
      do 3200 ll=1,norb                                                 
      l=l+ll                                                            
      fdiag(ll)=f(l)                                                    
3200  fdiag(ll+norb)=f(l)                                               
      if(.not.yion) go to 3600                                          
      nter=nter+1                                                       
      nn=nter+norb+norb+2                                               
      do 3500 l=nter,nn                                                 
      f(l)=0.                                                           
c                                                                       
c                                                                       
      aj(l)=0.                                                          
3500  ak(l)=0.                                                          
      fdiag(mnorb+1)=0.                                                 
      fdiag(mnorb+2)=0                                                  
3600  continue                                                          
      return                                                            
      end                                                               
