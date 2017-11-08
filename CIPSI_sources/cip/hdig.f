      function hdig(iii)                                                
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
c                  ce ssp fonction calcule l energie du determinant     
c                  numerotee  i  ,par rapport a l energie de l etat     
c                  fondamental                                          
c                   j    integrales coulombiennes  stockees             
c                        dans le tableau aj                             
c                   k    integrales d echanges     stockees             
c                        dans le tableau ak                             
c                   f    elements de l operateur de fock                
      nei=ne(iii)                                                       
c     write(6,*)(trou(k),k=1,50),(part(k),k=1,50)                       
c     write(6,'('' ne='',20i4)') ne                                     
      if(nei.eq.0) go to 1300                                           
      l1=nd(iii)                                                        
      l2=l1+nei                                                         
      l1=l1+1                                                           
      ab=0.                                                             
      do 1200 l=l1,l2                                                   
      n1p=part(l)                                                       
      n1t=trou(l)                                                       
      ab=ab+fdiag(n1p)-fdiag(n1t)                                       
      ns1t=ispin(n1t)                                                   
      ns1p=ispin(n1p)                                                   
      n1t=iorb(n1t)                                                     
      n1p=iorb(n1p)                                                     
c                                                                       
c                                                                       
      do 1100 j=l1,l                                                    
      n2t=trou(j)                                                       
      n2p=part(j)                                                       
      ns2t=ispin(n2t)                                                   
      ns2p=ispin(n2p)                                                   
      n2t=iorb(n2t)                                                     
      n2p=iorb(n2p)                                                     
      if(n2t.gt.n1t) go to 1010                                         
      n12t=num(n1t) +n2t                                                
      go to 1020                                                        
1010  n12t=num(n2t)+n1t                                                 
1020  if(n2p.gt.n1p) go to 1030                                         
      n12p=num(n1p)+n2p                                                 
      go to 1040                                                        
1030  n12p=num( n2p)+n1p                                                
1040  n1tp=num(n2p)+n1t                                                 
      n2tp=num(n1p)+n2t                                                 
      ca=aj(n12t)+aj(n12p)-aj(n1tp)-aj(n2tp)                            
      if(ns1t.eq.ns2t) ca=ca-ak(n12t)                                   
      if(ns1p.eq.ns2p) ca=ca-ak(n12p)                                   
      if(ns1t.eq.ns2p) ca=ca+ak(n1tp)                                   
      if(ns1p.eq.ns2t) ca=ca+ak(n2tp)                                   
      if(j.eq.l) ca=0.5*ca                                              
1100  ab=ab+ca                                                          
1200  continue                                                          
      hdig=ab                                                           
c                                                                       
      return                                                            
1300  hdig=0.d0                                                         
c                                                                       
      return                                                            
      end                                                               
