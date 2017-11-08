      subroutine deter                                                  
      implicit real*8(a-h,o-x,z),logical*1(y)                           
      include 'pshf.prm'
      character *20 title
      character*4 iwspin,iplus,imoins                                   
cCc       character*2 iwspin,iplus,imoins                                   
      character*8 appel,xoui,xnon,xt
      integer*4 ne,nd,trou,part,iorb,itsym,its,isytr,nesp,nespo         
      integer*4 imul,itm                                                
      integer*4 ntm(ndetz),npm(ndetz),ntd(2*ndetz),npd(2*ndetz),        
     *ntt(3*ndetz),npt(3*ndetz),ntq(4*ndetz),npq(4*ndetz),               
     *ntp(5*ndetz),npp(5*ndetz),nth(6*ndetz),nph(6*ndetz),              
     *nt7(7*ndetz),np7(7*ndetz),nto(8*ndetz),npo(8*ndetz)              
      integer*4 itab(72*ndetz)                                         
      dimension indtr(8),indpa(8)                                       
      dimension isydeg(nsymz,nsymz),ind(2*doa)              
      dimension kex(9),inv(2),itsyp(2*doa)                            
      dimension nonact(doa)
      dimension indis(doa,nsymz)                                      
      common  kf,km,kd,kt,kq,kp,kh,k7,ko,nt(9),np(9),mt(9),mp(9)        
      common/recla/iord(ndetz),mspin(metz)                              
      common/hist/nom(10),nam(10,metz),aval(10),tabmp(10,metz),         
     *tau,semp(metz),sevp(metz),tabvp(10,metz),nclass,ymoyen            
      common/det/my,ny,idet(400),nexst(100),vst(100),long(8)            
      common/spo/fdiag(2*doa),hdiag(ndetz),tocom(2*doa),coeff,      
     *ndef,ispin(2*doa),iwspin(2*doa),numero(metz),          
     *ne(ndetz+1),nd(ndetz+1),trou(ntpz),part(ntpz),iorb(2*doa),      
     *itsym(2*doa),its(nsymz,nsymz),yoc(2*doa),yocs(2*doa),       
     *ygen(ndetz),yocd(2*doa,ndetz),yspin                                    
      common/degene/nesp(ndetz),nespo(ndetz),isytr(ndetz,nsymz),kkbrd,  
     *ysaut(ndetz),imul(nmulz),yheff(metz*(metz+1)/2),yeff(ndetz,nsymz) 
      common/nature/itm(metz,nsymz),ityper                              
      common/activ/numac(doa),noac,nusym,nelac                        
      common/ic/c(metz*ndetz),e(ndetz),eb(metz),e2b(metz),emp(metz),    
     *teste(metz),e2mp(metz),dpivot,ietat,ishif,e2vp(metz),             
     *fmpb(2*doa),hefmp(nhefz),hefvp(nhefz),hefb(nhefz),              
     *psimp(metz),psivp(metz)                                           
      common/cipsi/escf,test,teseff,tdet,twdet,ndiag,norb,mnorb,
     *             noca,nocb,noc1,noc2,igela,igelb,ncf,nsym,isym,isz,
     *             ntrsy,metat,lecdet,maxdet,yprt,yion,ypertu,ybrd,yen,
     *             ywvi,ywvf,ytoul,ystkic,ywdeg
      common/tab/ ntm,npm,ntd,npd,ntt,npt,ntq,npq,ntp,npp,nth,nph       
     *,nt7,np7,nto,npo                                                  
      common/nom33/title
      common/tome/taum(metz),yselm,ysele,ydp                            
      common/infogv/nbdet(9),iord60(ndetz),ynoic                        
      common/orper/yreduc,yorper(doa)
      equivalence (itab(1),ntm(1))                                      
      equivalence(kex(1),kf)                                            
      data iplus,imoins/'+','-'/                                        
c                                                                       
      namelist/icinp/test,ndiag,isz,metat,tocom,yprt,yion,              
     *ybrd,ityper,numero,coeff,ygen,title,                              
     *tau,ymoyen,aval,nclass,ndiag,ytoul,lecdet,tdet,maxdet,teseff,     
     *nvr,numac,nelac,nusym,noac,                                       
     *ywvi,ywvf,ystkic,ywdeg,ytocom,iselec,                             
     *iord,twdet,mspin,ydp,ninact,nonact                                
c                                                                       
c     lecdet=1  lecture non automatique des determinants                
c            2  lecture automatique sur file 60 avec seuil tdet         
c            3  lecture automatique sur file 60 des                     
c               maxdet premiers determinants de la file                 
c            4  generation d un cas                                     
c            5  lecture de l ic ( vecteurs et energies sur file 59      
c               (determinants sur file 60 apres bdav ou autre)          
c                                                                       
c                                                                       
c     iselec=1     selection fonction d onde avec renormalisation       
c            2     selection fonction d onde sans renormalisation       
c            3     selection sur l energie                              
                                                                        
      data appel/'        '/                                            
      data xoui/' oui    '/                                             
      data xnon/' non    '/                                             
      lecref=0   
      yreduc=.false.
      ninact=0
      twdet=0.d0                                                        
      indtr(1)=0                                                        
      indtr(2)=2*ndetz                                                  
      indtr(3)=6*ndetz                                                  
      indtr(4)=12*ndetz                                                 
      indtr(5)=20*ndetz                                                 
      indtr(6)=30*ndetz                                                 
      indtr(7)=42*ndetz                                                 
      indtr(8)=56*ndetz                                                 
      indpa(1)=ndetz                                                    
      indpa(2)=4*ndetz                                                  
      indpa(3)=9*ndetz                                                  
      indpa(4)=16*ndetz                                                 
      indpa(5)=25*ndetz                                                 
      indpa(6)=36*ndetz                                                 
      indpa(7)=49*ndetz                                                 
      indpa(8)=64*ndetz                                                 
      igela=0                                                           
      igelb=0                                                           
      nsym=1                                                            
      title='                   '
      ywdeg=.true.                                                      
      norb=1                                                            
      noca=0                                                            
      nocb=0                                                            
      ntyp=0                                                            
      maxdet=ndetz                                                      
      ntrans=0                                                          
      noac=0                                                            
      nusym=1                                                           
      nelac=0                                                           
      ydp=.true.                                                       
cibm                                                                    
cibm      call errset(213,1,0,0,1)                                      
cibm                                                                    
cibm      call errset(217,1,0,0,1)                                      
      ntrsy=1    
      read(40,end=7788,err=7788) nsym,norb,noc,ntrans,                  
     1 (itsyp(i),i=1,norb),noa,ydp                                      
c     write(6,*) nsym,norb,noc,ntrans,(itsyp(i),i=1,norb)               
7788  ndegen=ntrans                                                     
c                                                                       
      inv(1)=norb                                                       
      inv(2)=-norb                                                      
      if(norb.eq.1) go to 99                                            
      ns=nsym                                                           
c                                                                       
610   format(' symetrie par orbitale ',(30i3))                          
c                                                                       
c                                                                       
      read(40) ((isydeg(i,j),i=1,ns),j=1,ns)                            
c     write(6,*) ((isydeg(i,j),i=1,ns),j=1,ns)                          
      do 1 i=1,ns                                                       
c                                                                       
      do 1 j=1,ns                                                       
      its(i,j)=isydeg(i,j)                                              
 612  format(20x,20i3)                                                  
1     continue
      do 11 i=1,norb                                                    
      numac(i)=i                                                        
      if(ntyp.ne.0) go to 11                                            
      itsym(i)=itsyp(i)                                                 
      itsyp(i+norb)=itsyp(i)                                            
11    itsym(i+norb)=itsyp(i)                                            
c                                                                       
c                                                                       
      if(ntyp.eq.0) ntyp=nsym                                           
      mnorb=norb+norb                                                   
c                                                                       
      do 4 i=1,ntyp                                                     
      k=0                                                               
      do 4 j=1,mnorb                                                    
      yocs(j)=.false.                                                   
      if(itsyp(j).ne.i) go to 4                                         
      k=k+1                                                             
      ind(j)=k                                                          
      indis(k,i)=j                                                      
4     continue                                                          
      if(ndegen.eq.0) go to 16                                          
      do 15 i=1,ndegen                                                  
      read(40) j,((isydeg(k,l),k=1,ntyp),l=1,j)                         
c     write(6,*) j,((isydeg(k,l),k=1,ntyp),l=1,j)                       
      do 14 l=1,j                                                       
      do 14 k=1,mnorb                                                   
      inskl=0                                                           
      nskl=isydeg(itsyp(k),l)                                           
      ysin=.false.                                                      
      if(nskl.ge.0) go to 13                                            
      nskl=-nskl                                                        
      ysin=.true.                                                       
13    continue                                                          
      if(nskl.ne.0) inskl=indis(ind(k),nskl)                            
      if(ysin) inskl=-inskl                                             
14    isytr(k,ntrsy+l)=inskl                                            
      ntrsy=ntrsy+j                                                     
15    continue                                                          
16    continue                                                          
      ytoul=.true.                                                      
      yion=.false.                                                      
      isym=0                                                            
      isz=0                                                             
      test=0.5                                                          
      teseff=-1.d0                                                      
      yprt=.false.                                                      
      ybrd=.false.                                                      
      do 4444 i=1,ndetz                                                 
      iord(i)=i                                                         
4444  ygen(i)=.false.                                                   
      do 4445 i=1,noc                                                   
4445  tocom(i)=2.d0                                                     
      escf=0.d0                                                         
      ystkic=.false.                                                    
      ityper=1                                                          
      long(1)=1                                                         
      long(2)=1                                                         
      long(3)=2                                                         
      long(4)=2                                                         
      long(5)=3                                                         
      long(6)=3                                                         
      long(7)=4                                                         
      long(8)=4                                                         
      aval(1)=0.128d-6                                                  
      ywvi=.false.                                                      
      ywvf=.true.                                                       
      ndiag=6                                                           
      aval(2)=0.64d-6                                                   
      aval(3)=0.32d-5                                                   
      aval(4)=0.16d-4                                                   
      aval(5)=0.8d-4                                                    
      aval(6)=0.4d-3                                                    
      aval(7)=0.2d-2                                                    
      aval(8)=0.1d-1                                                    
      aval(9)=0.5d-1                                                    
      aval(10)=1.d0                                                     
      tau=1.d-4                                                         
      ytocom=.false.                                                    
      ypertu=.true.                                                     
      ymoyen=.false.                                                    
      yspin=.false.                                                     
      iselec=1                                                          
      lecdet=1                                                          
      ynoic=.false.                                                     
      nclass=10                                                         
      metat=1                                                           
      do 500 i=1,metz                                                   
      mspin(i)=0                                                        
500   numero(i)=i                                                       
      yidg=.false.                                                      
      irest=0                                                           
      read(5,icinp,end=9500,err=9500)                                   
      if(teseff.lt.0.d0) teseff=0.05d0*tau                              
      if(lecdet.eq.5) then                                              
        ynoic=.true.                                                    
        tdet=0.d0                                                       
       endif              
      rewind 21                                                         
      if(ytocom) read(21) (tocom(i),i=1,norb)                           
      if(ityper.eq.0) ypertu=.false.                                    
      if(ityper.eq.1) then                                              
      yen=.false.                                                       
      endif                                                             
      if(ityper.eq.2) then                                              
      yen=.true.                                                        
      endif                                                             
      if(iselec.eq.1) then                                              
      yselm=.true.                                                      
      endif                                                             
      if(iselec.eq.2) then                                              
      yselm=.false.                                                     
      endif                                                             
      if(iselec.eq.3) then                                              
      ysele=.true.                                                      
      endif                                                             
      write(6,2100) title                                               
2100  format(//80(1h*),/,1x,a20,/,80(1h*))
3000  format(//)                                                        
      write(6,3010)                                                     
      write(6,3020)                                                     
3020  format(1x,'informations generales ')                              
      write(6,3010)                                                     
3010  format(1x,38('*'))                                                
      write(6,3000)                                                     
      do 3333 kk=1,10                                                   
      nom(kk)=0                                                         
      do 3333 ii=1,metat                                                
      tabmp(kk,ii)=0.d0                                                 
      tabvp(kk,ii)=0.d0                                                 
3333  nam(kk,ii)=0                                                      
      if(noca.eq.0) go to 2                                             
      if(nocb.eq.0) nocb=noca                                           
      nocb=nocb+norb                                                    
         go to 3                                                        
2     if(nocb.eq.0) nocb=noc                                            
      noca=nocb                                                         
      nocb=nocb+norb                                                    
3       continue                                                        
      igela=igela+1                                                     
      if(igelb.lt.norb) igelb=igela+norb                                
      if(noca.gt.norb) go to 99                                         
      if(.not.yion) go to 109                                           
      mnor1=mnorb+1                                                     
      itsym(mnor1)=nsym+1                                               
      nsym1=nsym+1                                                      
      its(nsym1,nsym1)=nsym1                                            
      isym=nsym+isym                                                    
      do 105 i=1,nsym                                                   
c                                                                       
      its(i,nsym1)=nsym1                                                
105   its(nsym1,i)=nsym1                                                
      iorb(mnor1)=norb+1                                                
      ispin(mnor1)=1                                                    
      iwspin(mnor1)=imoins                                              
      isytr(mnor1,1)=0                                                  
      if(isz.ne.0) isytr(mnor1,1)=0                                     
      do 107 l=2,ntrsy                                                  
      isytr(mnor1,l)=mnor1                                              
107   continue                                                          
109   continue                                                          
c                                                                       
      n=norb                                                            
      igele=igela                                                       
      do 150 i=1,norb                                                   
      iorb(i)=i                                                         
      iorb(i+norb)=i                                                    
      ispin(i)=0                                                        
      ispin(i+norb)=1                                                   
      iwspin(i)=iplus                                                   
150   iwspin(i+norb)=imoins                                             
      do 170 i=1,norb                                                   
      k=i+norb                                                          
      if(k.lt.igelb.and.i.gt.igela) go to 169                           
      if(k.gt.igelb.and.i.lt.igela) go to 169                           
      if(k.lt.nocb.and.i.gt.noca) go to 169                             
      if(k.gt.nocb.and.i.lt.noca) go to 169                             
      isytr(i,1)=k                                                      
      isytr(k,1)=i                                                      
      go to 170                                                         
169   isytr(i,1)=0                                                      
      isytr(k,1)=0                                                      
170   continue                                                          
      write(6,3500) norb                                                
3500  format(1x,'nombre d orbitales(norb):',i3)                         
      write(6,3510) mnorb                                               
3510  format(1x,'nombre de spinorbitales(mnorb):',i3)                   
      mocb=nocb-norb                                                    
      mgela=igela-1                                                     
      mgelb=igelb-norb-1                                                
      write(6,3520) noca                                                
3520  format(1x,'nombre de spinorbitales occuppees de spin alpha:',i3)  
3530  format(1x,'nombre de spinorbitales occuppees de spin beta:',i3)   
      write(6,3530) mocb                                                
      write(6,3560)                                                     
3560  format(1x,'nombre d occupation des orbitales(tocom):')            
      write(6,3570) (tocom(i),i=1,norb)                                 
3570  format(10(2x,f8.6))                                               
      write(6,3580) nsym                                                
3580  format(1x,'nombre de symetries pour les orbitales(nsym):',i2)     
3590  format(1x,'nombre d etats traites en perturbation(metat):',i3)    
      write(6,3600) isz                                                 
3600  format(1x,'valeur de sz:',i2)                                     
      xt=xnon                                                           
      if(yion) xt=xoui                                                  
      write(6,3700) xt                                                  
3700  format(1x,'excitation vers l infini:',a8)                         
      xt=xnon                                                           
      write(6,3620) coeff                                               
3620  format(1x,'seuls sont generateurs les determinants de poids',     
     *' superieurs a ',f9.6,'sur au moins un etat utile')               
      write(6,3590) metat                                               
      if(metat.gt.metz) goto 777                                        
      write(6,3650)                                                     
3650  format(1x,'numero des etats traites en perturbation',             
     *'(numero)')                                                       
      write(6,3660) (numero(ii),ii=1,metat)                             
3660  format(20(2x,i2))                                                 
      write(6,3630) test                                                
3630  format(1x,'test absolu (test):',f9.6)                             
      if(ypertu) xt=xoui                                                
      write(6,3670) xt                                                  
3670  format(1x,'perturbation moller plesset barycentrique:',a8)        
      xt=xnon                                                           
      if(yen) xt=xoui                                                   
      write(6,3680) xt                                                  
3680  format(1x,'perturbation epstein-nesbet valeur propre:',a8)        
      write(6,3690) xt                                                  
3690  format(1x,'perturbation epstein -nesbet barycentrique:',a8)       
      go to 3810                                                        
3800  write(6,3830) test                                                
3830  format(' perturbation epstein-nesbet barycentrique seulement',/,  
     *' test: ',f9.6)                                                   
3810  continue                                                          
      if(ymoyen) xt=xoui                                                
      write(6,4100) xt                                                  
4100  format(1x,'histogramme et ecriture des moyens sur la file 60:',   
     *a8)                                                               
      if(.not.ymoyen) go to 4500                                        
      write(6,4110) tau                                                 
4110  format(1x,'seuils de selections des determinants moyens: ',f8.6)  
      aa=0.d0                                                           
      write(6,4120) aa,(aval(ii),ii=1,nclass)                           
4120  format(1x,'bornes des classes de l histogramme:',                 
     */,1x,11(2x,e8.2))                                                 
4500  continue                                                          
      xt=xnon                                                           
      if(ybrd) xt=xoui                                                  
      write(6,3710) xt                                                  
3710  format(1x,'hamiltonien effectif qdpt(ybrd):',a8)                  
      xt=xnon                                                           
      xt=xoui                                                           
      write(6,3720) xt                                                  
3720  format(1x,'resultats variationnels et energies gardes sur',       
     *' la file 04:',a8)                                                
      write(6,3000)                                                     
c                                                                       
c   lecdet=1 lecture des determinants en donnee                         
c                                                                       
      if(lecdet.eq.1) then                                              
      write(6,*)                                                        
      write(6,*) ' espace variationnel fourni en donnees'               
      write(6,*)                                                        
      call       lecnew(ntm,npm,ntd,npd,ntt,npt,ntq,npq,ntp,npp,nth,nph,
     *nt7,np7,nto,npo)                                                  
      endif                                                             
c                                                                       
c    lecdet=2 lecture des determinants sur file 60 avec seuil tdet      
c                                                                       
c    lecdet=3 lecture des maxdet premiers determinants de la file60     
c                                                                       
                                                                        
      if((lecdet.eq.2).or.(lecdet.eq.3).or.(lecdet.eq.5)) then          
      write(6,*)                                                        
      write(6,*) ' espace variationnel issu de la file 60'              
      if((lecdet.eq.2).or.(lecdet.eq.5))                                
     1 write(6,*) ' selection par seuil tdet= ',tdet                    
      if(lecdet.eq.3) write(6,*) ' selection par nombre  maxdet=',      
     *maxdet                                                            
      write(6,*)                                                        
      call       lec60(ntm,npm,ntd,npd,ntt,npt,ntq,npq,ntp,npp,nth,nph, 
     *nt7,np7,nto,npo)                                                  
         write(6,*) ' type des excitations'                             
         write(6,3684) ((kl-1),kl=1,9)                                  
3684  format(1x,' degre d excitation ',9(1x,i6))                        
         do 198 kl=1,9                                                  
            if (kl.ne.1) then                                           
               nbdet(kl)=kex(kl)/(kl-1)                                 
            else                                                        
               nbdet(kl)=kex(kl)                                        
            endif                                                       
198      continue                                                       
            write(6,3682) (nbdet(kl),kl=1,9)                            
3682  format(1x,' determinants       ',9(1x,i6))                        
      endif                                                             
c                                                                       
c    lecdet=4 generation d'un espace  cas                               
c                                                                       
      if(lecdet.eq.4)then                                               
      call     cas(ntm,npm,ntd,npd,ntt,npt,ntq,npq,ntp,npp,nth,nph,     
     *nt7,np7,nto,npo)                                                  
      endif                                                             
      do i=1,norb
       yorper(i)=.true.
       enddo
      if(ninact.gt.0) then
       write(6,*) ' troncation de l espace actif perturbatif '
       write(6,*) ' les determinants ne contenant que les orbitales'
       write(6,3660) (nonact(i),i=1,ninact)
       write(6,*) ' ne sont pas pris en compte dans la perturbation'
       yreduc=.true.
       do i=1,ninact
        yorper(nonact(i))=.false.
        enddo
       endif
       write(6,*)'yorper' , (yorper(i),i=1,norb)
      ncf=0                                                             
      if(ntrsy.lt.2) go to 191                                          
      do 190 n =2,ntrsy                                                 
      do 189 k=igela,mnorb                                              
      iskn=isytr(k,n)                                                   
      iskn=abs(iskn)                                                    
      if(iskn.lt.igela.or.iskn.gt.mnorb) go to 188                      
      if(k.gt.noca) go to 185                                           
      if(iskn.gt.noca) go to 188                                        
      go to 189                                                         
185   if(k.gt.norb) go to 186                                           
      if(iskn.le.noca.or.iskn.gt.norb) go to 188                        
      go to 189                                                         
186   if(k.gt.nocb) go to 187                                           
      if(iskn.lt.igelb.or.iskn.gt.nocb) go to 188                       
      go to 189                                                         
187   if(iskn.gt.nocb) go to 189                                        
188   isytr(k,n)=0                                                      
189   continue                                                          
190   continue                                                          
191   continue                                                          
c                                                                       
      do 7766 kk=1,ndetz                                                
      do 7766 ii=1,ntrsy                                                
7766  yeff(kk,ii)=.false.                                               
      mesp=1                                                            
      nespo(1)=0                                                        
770   if(kf.ne.0) call transf(0)                                        
      do 80 i=2,9                                                       
      nec=i-1                                                           
      nex=kex(i)/nec                                                    
      if(nex.eq.0) go to 80                                             
      monot=indtr(nec)                                                  
      monop=indpa(nec)                                                  
      do 79 j=1,nex                                                     
      ncfo=ncf+1                                                        
      ysaut(ncfo)=.false.                                               
      do 20 l=1,nec                                                     
      monot=monot+1                                                     
      monop=monop+1                                                     
      nt(l)=itab(monot)                                                 
20    np(l)=itab(monop)                                                 
      call transf(nec)                                                  
      do 31 l=1,nec                                                     
      mt(l)=nt(l)                                                       
31    mp(l)=np(l)                                                       
      if(ncf.lt.ncfo) go to 79                                          
      if(isz.ne.0) yeff(ncfo,1)=.true.                                  
      if(isz.ne.0) go to 30                                             
      do 25 l=1,nec                                                     
      ntl=isytr(nt(l),1)                                                
      npl=isytr(np(l),1)                                                
      if(ntl.eq.0)go to 30                                              
      if(npl.eq.0)go to 30                                              
      nt(l)=ntl                                                         
25    np(l)=npl                                                         
      call transf(nec)                                                  
      if(ncf.le.ncfo) go to 30                                          
      yeff(ncfo,1)=.true.                                               
      yeff(ncf,1)=.true.                                                
      nespo(ncf)=1                                                      
30    if(ndegen.eq.0) go to 75                                          
      do 70 ll=2,ntrsy                                                  
      ncfa=ncf                                                          
      yeff(ncfo,ll)=.false.                                             
      isigne=1                                                          
      ldeg=0                                                            
      do 34 l=1,nec                                                     
      nisy=isytr(mt(l),ll)                                              
      if(nisy.lt.0) isigne=-isigne                                      
      n=abs(nisy)                                                       
      if(n.eq.0) go to 70                                               
      nt(l)=n                                                           
      if(n.eq.mt(l)) go to 33                                           
      ldeg=ldeg+1                                                       
33    nisy=isytr(mp(l),ll)                                              
      if(nisy.lt.0) isigne=-isigne                                      
      n=abs(nisy)                                                       
      if(n.eq.0) go to 70                                               
      if(n.ne.mp(l)) ldeg=ldeg+1                                        
      np(l)=n                                                           
34     continue                                                         
      if(ldeg.eq.0.or.((2*(ldeg/2)).ne.ldeg))go to 40                   
      call transf(nec)                                                  
      if(ncf.le.ncfa) go to 65                                          
      yeff(ncfo,ll)=.true.                                              
      yeff(ncf,ll)=.true.                                               
      nespo(ncf)=isigne                                                 
      if(isz.ne.0) go to 65                                             
      do 37 l=1,nec                                                     
      ntl=isytr(nt(l),1)                                                
      npl=isytr(np(l),1)                                                
      if(ntl.eq.0) go to 40                                             
      if(npl.eq.0) go to 40                                             
      nt(l)=ntl                                                         
37    np(l)=npl                                                         
      call transf(nec)                                                  
      if(ncf.le.(ncfa+1)) go to 65                                      
      yeff(ncf,1)=.true.                                                
      yeff(ncf,ll)=.true.                                               
      nespo(ncf)=isigne                                                 
40    continue                                                          
65    continue                                                          
70    continue                                                          
75    if(ncfo.ge.ncf) go to 79                                          
      if(ncf-ncfo) 79,76,77                                             
76    nesp(ncf)=0                                                       
      go to 79                                                          
77    nesp(ncfo)=ncf-ncfo                                               
      ncfo=ncfo+1                                                       
      do 78 l=ncfo,ncf                                                  
      ysaut(l)=.true.                                                   
c                                                                       
c                                                                       
78    nesp(l)=0                                                         
79    continue                                                          
80    continue                                                          
      if(ncf.eq.0) go to 90                                             
160   format('0 dimension du sous espace variationnel:',i4,/)           
      mncf=ncf*ndiag                                                    
      if(mncf.le.(metz*ndetz)) go to 4530                               
      write(6,4535) mncf                                                
4535  format(/,1x,'dimension necessaire pour la matrice des vecteurs: ' 
     *,i10)                                                             
      stop                                                              
4530  continue                                                          
      return                                                            
777   write(6,778)                                                      
778   format(1x,/,1x,'trop d etats demandes en perturbation',//)        
      stop                                                              
90    write(6,91)                                                       
91    format(' erreur dans les donnees:aucun determinant dans le sous', 
     *' espace variationnel',/)                                         
      stop                                                              
99    write(6,98)                                                       
98    format(' erreur dans les informations generales')                 
9500  write(6,9501)                                                     
9501  format(' erreur dans la namelist icinp',/)                        
      stop                                                              
      end                                                               
