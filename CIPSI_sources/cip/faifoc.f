      subroutine faifoc                                                 
      implicit real*8 (a-h,o-x,z),logical*1 (y)                         
      include 'pshf.prm'
      common/jkf/ f(doa*(doa+1)/2),aj(doa*(doa+1)/2),           
     2ak(doa*(doa+1)/2)                                             
      dimension nv1(2),nv2(2)         
      equivalence(nv11,nv1(1)),(nv12,nv1(2)),(nv21,nv2(1)),(nv22,nv2(2))
      common/spo/fdiag(2*doa),hdiag(ndetz),tocom(2*doa),coeff,      
     1ndef,ispin(2*doa),iwspin(2*doa),numero(metz),          
     2ne(ndetz+1),nd(ndetz+1),trou(ntpz),part(ntpz),iorb(2*doa),      
     3itsym(2*doa),its(nsymz,nsymz),yoc(2*doa),yocs(2*doa),       
     4ygen(ndetz),yocd(2*doa,ndetz),yspin,yef(ndetz,nsymz)                   
      common/degene/nesp(ndetz),nespo(ndetz),isytr(ndetz,nsymz),kkbrd,  
     1ysaut(ndetz)                                                      
      common/cipsi/escf,test,teseff,tdet,twdet,ndiag,norb,mnorb,
     1             noca,nocb,noc1,noc2,igela,igelb,ncf,nsym,isym,isz,
     2             ntrsy,metat,lecdet,maxdet,yprt,yion,ypertu,ybrd,yen,
     3             ywvi,ywvf,ytoul,ystkic,ywdeg
      integer*4 ne,nd,trou,part,iorb,iwspin,itsym,its,isytr,nesp,nespo  
      common/ic/c(metz*ndetz),e(ndetz),eb(metz),e2b(metz),emp(metz),    
     1teste(metz),e2mp(metz),dpivot,ietat,ishif,e2vp(metz),             
     2fmpb(2*doa),hefmp(nhefz),hefvp(nhefz),hefb(nhefz),              
     3psimp(metz),psivp(metz)                                           
      num(i)=i*(i-1)/2                                                  
      do 4500 l=1,norb                                                  
      t=fdiag(l)                                                        
      lj=num(l)                                                         
      do 4400 j=1,norb                                                  
      lj=lj+1                                                           
      if(j.gt.l) lj=num(j)+l                                            
4400  t=t+(tocom(j)+tocom(j+norb))*aj(lj)                               
      tt=t                                                              
      lj=num(l)                                                         
      do 4450 j=1,norb                                                  
      lj=lj+1                                                           
      if(j.gt.l) lj=num(j)+l                                            
      t=t-tocom(j)*ak(lj)                                               
4450  tt=tt-tocom(j+norb)*ak(lj)                                        
      fmpb(l)=t                                                         
4500  fmpb(l+norb)=tt                                                   
c                                                                       
c                                                                       
c                                                                       
      return                                                            
      end                                                               
