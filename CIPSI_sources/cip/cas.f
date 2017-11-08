      subroutine cas                                                    
     *(ntm,npm,ntd,npd,ntt,npt,ntq,npq,ntp,npp,nth,nph,nt7,np7,nto,npo) 
      implicit real*8 (a-h,p,r-x,z),logical*1(y)                        
      include 'pshf.prm'
      parameter (niiz=25000)
      integer*4 ne,nd,trou,part,iorb,iwspin,itsym,its  
      integer*4 noma,nomb,nta,npa,ntb,npb,nts,nps,mtp 
      integer*4 ntm(ndetz),npm(ndetz),ntd(2*ndetz),npd(2*ndetz),        
     *ntt(3*ndetz),npt(3*ndetz),ntq(4*ndetz),npq(4*ndetz),               
     *ntp(5*ndetz),npp(5*ndetz),nth(6*ndetz),nph(6*ndetz),              
     *nt7(7*ndetz),np7(7*ndetz),nto(8*ndetz),npo(8*ndetz)              
      dimension istr(doa),noma(niiz,10),nomb(niiz,10),                
     *yalfa(doa),ybeta(doa),nta(niiz),npa(niiz),                    
     *ntb(niiz),npb(niiz),mtp(10),nts(10),nps(10)                       
      common/spo/fdiag(2*doa),hdiag(ndetz),tocom(2*doa),coeff,      
     *ndef,ispin(2*doa),iwspin(2*doa),numero(metz),          
     *ne(ndetz+1),nd(ndetz+1),trou(ntpz),part(ntpz),iorb(2*doa),      
     *itsym(2*doa),its(nsymz,nsymz),yoc(2*doa),yocs(2*doa),       
     *ygen(ndetz),yocd(2*doa,ndetz),yspin                                    
      common/cipsi/escf,test,teseff,tdet,twdet,ndiag,norb,mnorb,
     *             noca,nocb,noc1,noc2,igela,igelb,ncf,nsym,isym,isz,
     *             ntrsy,metat,lecdet,maxdet,yprt,yion,ypertu,ybrd,yen,
     *             ywvi,ywvf,ytoul,ystkic,ywdeg
      common/activ/numac(doa),noac,nusym,nelac                        
      common kf,km,kd,kt,kq,kp,kh,k7,ko                                 
c                                                                       
c   generation des determinants correspondant a l ic complete de        
c   nalfa electrons alpha et nbeta electrons beta dans noac orbitales   
c                                                                       
Ccc   RG : modified on 02/11/07 : parameter niiz set to 3000. Was 100.
Ccc   RG : modified on 14/01/08 : parameter niiz set to 10000.
Ccc   RG : modified on 13/02/08 : netz, in psf.prm, set to 16000.
Ccc   RG : modified on 28/11/08 : parameter niiz set to 12000.
Ccc   RV : modified on 27/01/16 : parameter niiz set to 25000.
      write(6,*)                                                        
      write(6,*) ' espace variationnel: cas '                           
      write(6,*) ' electrons actifs : ',nelac                           
      write(6,*) ' orbitales actives: ',noac                            
      write(6,*) ' numeros des orbitales actives  '                     
      write(6,999) (numac(i),i=1,noac)                                  
999   format(30(1x,i3))                                                 
      write(6,*) ' symetrie des determinants  ',nusym                   
      write(6,*)                                                        
      kf=0                                                              
      km=0                                                              
      kd=0                                                              
      kt=0                                                              
      kq=0                                                              
      kp=0                                                              
      kh=0                                                              
      k7=0                                                              
      ko=0                                                              
c                                                                       
c   electrons alpha                                                     
c                                                                       
      nalfa=(nelac+1)/2                                                 
      nbeta=nalfa                                                       
      if(nelac.gt.20) then                                              
      write(6,*) ' trop d electrons actifs  maximum 20'                 
      stop                                                              
      endif                                                             
      noc=nalfa                                                         
      if(yion) nalfa=nalfa-1                                            
      nca=0                                                             
      ytest=.false.                                                     
c                                                                       
c  electrons  alfa                                                      
c                                                                       
      do 5 k=1,nalfa                                                    
5     istr(k)=k                                                         
      if(nalfa.ne.0) istr(nalfa)=nalfa-1                                
      ii=0                                                              
      i=nalfa+1                                                         
10    i=i-1                                                             
      if(i.eq.0) go to 18                                               
14    newoc=istr(i)+1                                                   
      if(newoc.gt.(noac-nalfa+i)) go to 10                              
      istr(i)=newoc                                                     
      do 16 k=i+1,nalfa                                                 
16    istr(k)=istr(k-1)+1                                               
      i=nalfa                                                           
      ii=ii+1                                                           
      if(ii.gt.niiz) then                                               
      write(6,*) '  trop de combinaisons alpha: max=',niiz
      stop                                                              
      endif                                                             
      do 12 k=1,nalfa                                                   
12    noma(ii,k)=istr(k)                                                
      go to 14                                                          
18    continue                                                          
      nca=ii                                                            
      if(ytest) write(6,*) ' nbre  cbn alfa',nca                        
c                                                                       
c   electrons beta                                                      
c                                                                       
      do 15 k=1,nbeta                                                   
15    istr(k)=k                                                         
      istr(nbeta)=nbeta-1                                               
      ii=0                                                              
      i=nbeta+1                                                         
20    i=i-1                                                             
      if(i.eq.0) go to 28                                               
24    newoc=istr(i)+1                                                   
      if(newoc.gt.(noac-nbeta+i)) go to 20                              
      istr(i)=newoc                                                     
      do 26 k=i+1,nbeta                                                 
26    istr(k)=istr(k-1)+1                                               
      i=nbeta                                                           
      ii=ii+1                                                           
      if(ii.gt.niiz) then                                               
      write(6,*) '  trop de combinaisons beta:max=10000'  ! RG 2 novembre 2007                  
      stop                                                              
      endif                                                             
      do 22 k=1,nbeta                                                   
22    nomb(ii,k)=istr(k)                                                
      go to 24                                                          
28    continue       
      ncb=ii                                                            
      if(ytest) write(6,*) ' nbre  cbn alfa',ncb                        
      if(ytest) then                                                    
      write(6,*) 'noma'                                                 
      do 776 i=1,nca                                                    
776   write(6,*) (noma(i,k),k=1,nalfa)                                  
      write(6,*) 'nomb'                                                 
      do 777 i=1,ncb                                                    
777   write(6,*) (nomb(i,k),k=1,nbeta)                                  
      endif                                                             
      do 300 ia=1,max(1,nca)                                            
      itp=0                                                             
      do 305 k=1,noac                                                   
305   yalfa(k)=.false.                                                  
      do 310 k=1,nalfa                                                  
      kk=noma(ia,k)                                                     
310   yalfa(kk)=.true.                                                  
      do 315 k=1,noc                                                    
      if(.not.yalfa(k)) then                                            
      itp=itp+1                                                         
      mtp(itp)=numac(k)                                                 
      endif                                                             
315   continue                                                          
      nta(ia)=itp                                                       
      do 320 k=noc+1,noac                                               
      if(yalfa(k)) then                                                 
      itp=itp+1                                                         
      mtp(itp)=numac(k)                                                 
      endif                                                             
320   continue                                                          
      npa(ia)=itp-nta(ia)                                               
      if(itp.eq.0) go to 300                                            
      do 325 k=1,itp                                                    
325   noma(ia,k)=mtp(k)                                                 
300   continue                                                          
      do 400 ib=1,ncb                                                   
      itp=0                                                             
      do 405 k=1,noac                                                   
405   ybeta(k)=.false.                                                  
      do 410 k=1,nbeta                                                  
      kk=nomb(ib,k)                                                     
410   ybeta(kk)=.true.                                                  
      do 415 k=1,noc                                                    
      if(.not.ybeta(k)) then                                            
      itp=itp+1                                                         
      mtp(itp)=numac(k)                                                 
      endif                                                             
415   continue                                                          
      ntb(ib)=itp                                                       
      do 420 k=noc+1,noac                                               
      if(ybeta(k)) then                                                 
      itp=itp+1                                                         
      mtp(itp)=numac(k)                                                 
      endif                                                             
420   continue                                                          
      npb(ib)=itp-ntb(ib)                                               
      if(itp.eq.0) go to 400                                            
      do 425 k=1,itp                                                    
425   nomb(ib,k)=mtp(k)                                                 
400   continue                                                          
      if(ytest) then                                                    
      write(6,*) (nta(i),i=1,nca)                                       
      write(6,*) (npa(i),i=1,nca)                                       
      do 34 i=1,nca                                                     
      nt0=nta(i)+npa(i)                                                 
34    write(6,*) (noma(i,j),j=1,nt0)                                    
      write(6,*) (ntb(i),i=1,ncb)                                       
      write(6,*) (npb(i),i=1,ncb)                                       
      do 35 i=1,ncb                                                     
      nt0=ntb(i)+npb(i)                                                 
35    write(6,*) (nomb(i,j),j=1,nt0)                                    
      endif                                                             
c                                                                       
c   construction des determinants et test de symetrie                   
c                                                                       
      ndet=0                                                            
      ndi=0                                                             
      do 500 ia=1,max(1,nca)                                            
      kta=nta(ia)                                                       
      kpa=npa(ia)                                                       
      if(kta.eq.0) go to 507                                            
      do 505 k=1,kta                                                    
505   nts(k)=noma(ia,k)                                                 
507   if(kpa.eq.0) go to 509                                            
      do 510 k=1,kpa                                                    
510   nps(k)=noma(ia,k+kta)                                             
509   continue                                                          
      do 600 ib=1,ncb                                                   
      ktb=ntb(ib)                                                       
      kpb=npb(ib)                                                       
      if(ktb.eq.0) go to 607                                            
      do 605 k=1,ktb                                                    
605   nts(k+kta)=nomb(ib,k)+norb                                        
607   if (kpb.eq.0) go to 609                                           
      do 610 k=1,kpb                                                    
610   nps(k+kpa)=nomb(ib,k+ktb)+norb                                    
609   continue                                                          
      if(yion) then                                                     
      kpb=kpb+1                                                         
      nps(kpb+kpa)=norb+norb+1                                          
      endif                                                             
      nect=kta+ktb                                                      
      necp=kpa+kpb                                                      
      if(nect.ne.necp) then                                             
      write(6,*) '  attention necp et nect sont differents  '           
      write(6,*) '  nect,necp   ',nect,necp                             
      stop                                                              
      endif                                                             
      nec=necp                                                          
      is=1                                                              
      do 700 i=1,nec                                                    
      if(nec.eq.0) then                                                 
      is=1                                                              
      go to 725                                                         
      endif                                                             
      npi=nps(i)                                                        
      jps=itsym(npi)                                                    
      if(npi.eq.norb+norb+1) jps=1                                      
      is=its(is,itsym(nts(i)))                                          
700   is=its(is,jps)                                                    
725   if(nusym.eq.0) nusym=is                                           
      if(is.ne.nusym) then                                              
      go to 600                                                         
      endif                                                             
      ndet=ndet+1                                                       
      if(ndet.gt.ndetz) then                                            
      write(6,*)                                                        
      write(6,*) ' espace actif trop grand '                            
      write(6,*) ' maximum=',ndetz                                         
      stop                                                              
cCc       write(6,*) ndet 
      endif                                   

      i=nec+1                                                           
c                                                                       
      go to (820,910,920,930,940,950,960,970,980) ,i                    
c                                                                       
c                                                                       
c fondamental                                                           
820   if(isz.ne.0) go to 800                                            
      kf=kf+1                                                           
      if(kf.eq.1) go to 800                                             
      kf=kf-1                                                           
825   write(6,826)                                                      
826   format(' fondamental deja present')                               
      go to 800                                                         
c                                                                       
 910  continue                                                          
      km=km+1                                                           
      ntm(km)=nts(1)                                                    
      npm(km)=nps(1)                                                    
      go to 800                                                         
 920  continue                                                          
      do 925 i=1,2                                                      
      kd=kd+1                                                           
      ntd(kd)=nts(i)                                                    
      npd(kd)=nps(i)                                                    
 925  continue                                                          
      go to 800                                                         
 930  continue                                                          
      do 935 i=1,3                                                      
      kt=kt+1                                                           
      ntt(kt)=nts(i)                                                    
      npt(kt)=nps(i)                                                    
 935  continue                                                          
      go to 800                                                         
 940  continue                                                          
      do 945 i=1,4                                                      
      kq=kq+1                                                           
      ntq(kq)=nts(i)                                                    
      npq(kq)=nps(i)                                                    
 945  continue                                                          
      go to 800                                                         
 950  continue                                                          
      do 955 i=1,5                                                      
      kp=kp+1                                                           
      ntp(kp)=nts(i)                                                    
      npp(kp)=nps(i)                                                    
 955  continue                                                          
      go to 800                                                         
 960  continue                                                          
      do 965 i=1,6                                                      
      kh=kh+1                                                           
      nth(kh)=nts(i)                                                    
      nph(kh)=nps(i)                                                    
 965  continue                                                          
      go to 800                                                         
 970  continue                                                          
      do 975 i=1,7                                                      
      k7=k7+1                                                           
      nt7(k7)=nts(i)                                                    
      np7(k7)=nps(i)                                                    
 975  continue                                                          
      go to 800                                                         
 980  continue                                                          
      do 985 i=1,8                                                      
      ko=ko+1                                                           
      nto(ko)=nts(i)                                                    
      npo(ko)=nps(i)                                                    
 985  continue                                                          
 800  continue                                                          
600   continue                                                          
500   continue                                                          
      write(6,*) '   determinants engendres ',ndet                      
      return                                                            
      end                                                               
