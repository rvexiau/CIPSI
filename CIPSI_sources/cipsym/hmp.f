      function hmp(iiii)                                                
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
      common/cipsi/escf,test,teseff,tdet,twdet,ndiag,norb,mnorb,
     1             noca,nocb,noc1,noc2,igela,igelb,ncf,nsym,isym,isz,
     2             ntrsy,metat,lecdet,maxdet,yprt,yion,ypertu,ybrd,yen,
     3             ywvi,ywvf,ytoul,ystkic,ywdeg
      integer*4 ne,nd,trou,part,iorb,iwspin,itsym,its,isytr,nesp,nespo  
      common/ic/c(metz*ndetz),e(ndetz),eb(metz),e2b(metz),emp(metz),    
     1teste(metz),e2mp(metz),dpivot,ietat,ishif,e2vp(metz),             
     2fmpb(2*doa),hefmp(nhefz),hefvp(nhefz),hefb(nhefz),              
     3psimp(metz),psivp(metz)                                           
      common/truc/indic(doa*(doa+1)/2),jndic(doa*(doa+1)/2),    
     1ndeb(500),nad(2000),kt(2000),                                     
     2nbo(99)                                                           
      integer*4 indic,jndic,nad,kt                                      
      logical iden,iprin1                                               
      num(i)=i*(i-1)/2                                                  
      nei=ne(iiii)                                                      
c                                                                       
      if(nei.eq.0) go to 2300                                           
c                                                                       
      l1=nd(iiii)                                                       
      l2=l1+nei                                                         
      l1=l1+1                                                           
      ac=0.d0                                                           
      do 2200 l=l1,l2                                                   
      n1p=part(l)                                                       
      n1t=trou(l)                                                       
2200  ac=ac+fmpb(n1p)-fmpb(n1t)                                         
      hmp=ac                                                            
      return                                                            
2300  hmp=0.d0                                                          
      return                                                            
      end                                                               
