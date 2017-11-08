      subroutine brdnee                                                 
      implicit real*8(a-h,o-x,z),logical*1(y)                           
      include 'pshf.prm'
c                                                                       
c  pour utiliser la dimension maximale, ce parametre doit etre egal a   
c  nmulz= nhefz*(nombre de generateurs non equivalents)                 
c                                                                       
      common/degene/nesp(ndetz),nespo(ndetz),isytr(ndetz,nsymz),kkbrd,  
     1ysaut(ndetz),imul(nmulz),yheff(metz*(metz+1)/2),yef(ndetz,nsymz)  
      integer*4 isytr,nesp,nespo,imul,itm,itk,nti,ntl,npl               
      common/cipsi/escf,test,teseff,tdet,twdet,ndiag,norb,mnorb,
     1             noca,nocb,noc1,noc2,igela,igelb,ncf,nsym,isym,isz,
     2             ntrsy,metat,lecdet,maxdet,yprt,yion,ypertu,ybrd,yen,
     3             ywvi,ywvf,ytoul,ystkic,ywdeg
      common/ic/c(metz*ndetz),e(ndetz),eb(metz),e2b(metz),emp(metz),    
     1teste(metz),e2mp(metz),dpivot,ietat,ishif,e2vp(metz),             
     2fmpb(2*doa),hefmp(nhefz),hefvp(nhefz),hefb(nhefz),              
     3psimp(metz),psivp(metz)                                           
      common/spo/fdiag(2*doa),hdiag(ndetz),tocom(2*doa),coeff,      
     1ndef,ispin(2*doa),iwspin(2*doa),numero(metz),          
     2ne(ndetz+1),nd(ndetz+1),trou(ntpz),part(ntpz),iorb(2*doa),      
     3itsym(2*doa),its(nsymz,nsymz),yoc(2*doa),yocs(2*doa),       
     4ygen(ndetz),yocd(2*doa,ndetz),yspin                                    
      common/nature/itm(metz,nsymz)                                     
      integer*4 ne,nd,trou,part,iorb,iwspin,itsym,its,npk,ntk           
      dimension nti(nsymz),itk(ndetz,10)                                
      dimension ntk(10),npk(10),ntl(10),npl(10)                         
      equivalence (itk(1,1),imul(1))                                    
      equivalence (itm(1,1),ntk(1)),(itm(11,1),npk(1))                  
      equivalence (ntl(1),nti(1)),(npl(1),nti(11))                      
c                                                                       
c   determination des transformes itk(k,l) des determinants k de s par  
c   les operations de symetrie l (numeros et signes)                    
c                                                                       
      l1=1                                                              
      yisz0=isz.eq.0                                                    
      if(.not.yisz0) l1=2                                               
      do 100 k=1,ncf                                                    
      if(ysaut(k)) go to 100                                            
      nespk=nesp(k)+1                                                   
      nec=ne(k)                                                         
      if(nec.ne.0) go to 105                                            
      do 195 l=l1,ntrsy                                                 
195   itk(k,l)=k                                                        
      go to 100                                                         
105   nec1=nec-1                                                        
      do 110 kg=1,nespk                                                 
      kk=k+kg-1                                                         
      ndk=nd(kk)                                                        
      do 175 i=1,nec                                                    
      ntk(i)=trou(i+ndk)                                                
175   npk(i)=part(i+ndk)                                                
      do 120 l=l1,ntrsy                                                 
      isigne=1                                                          
      ncr=1                                                             
      do 125 i=1,nec                                                    
      nl=ntk(i)                                                         
      np=npk(i)                                                         
      ntli=isytr(nl,l)                                                  
      npli=isytr(np,l)                                                  
      if((ntli.ne.0).and.(npli.ne.0)) go to 155                         
      itk(kk,l)=0                                                       
      go to 120                                                         
155   continue                                                          
      if(ntli.le.0) isigne=-isigne                                      
      if(npli.le.0) isigne=-isigne                                      
      ntl(i)=abs(ntli)                                                  
      npl(i)=abs(npli)                                                  
125   continue                                                          
      if(yef(k,l)) go to 145                                            
       ii=kk                                                            
      go to 135                                                         
145   do 130 ig=1,nespk                                                 
      ii=k+ig-1                                                         
      do 140 i=1,nec                                                    
      it=ntl(i)                                                         
      ip=npl(i)                                                         
      if(.not.yocd(it,ii)) go to 130                                    
      if(.not.yocd(ip,ii)) go to 130                                    
140   continue                                                          
      go to 135                                                         
130   continue                                                          
      if(nec.eq.1) go to 198                                            
135   do 160 i=1,nec1                                                   
      ip=i+1                                                            
      do 150 j=ip,nec                                                   
      if(ntl(j).gt.ntl(i)) go to 165                                    
      n=ntl(j)                                                          
      ntl(j)=ntl(i)                                                     
      ntl(i)=n                                                          
      ncr=-ncr                                                          
165   if(npl(j).gt.npl(i)) go to 150                                    
      n=npl(j)                                                          
      npl(j)=npl(i)                                                     
      npl(i)=n                                                          
      ncr=-ncr                                                          
150   continue                                                          
160   continue                                                          
198   if(isigne.lt.0) ii=-ii                                            
      if(ncr.lt.0) ii=-ii                                               
      itk(kk,l)=ii                                                      
120   continue                                                          
110   continue                                                          
100   continue                                                          
      if(.not.yprt) go to 170                                           
      write(6,*) ' k '                                                  
       write(6,180) (k,k=1,ncf)                                         
      write(6,*) ' itk'                                                 
      do 185 l=l1,ntrsy                                                 
      write(6,*) ' l=',l                                                
185   write(6,180) (itk(k,l),k=1,ncf)                                   
180   format(20(1x,i4))                                                 
170   continue                                                          
c                                                                       
c  determination du signe itm(m,l) des transformes des etats m          
c  par les operations de symetrie l                                     
c                                                                       
      do 200 m=1,metat                                                  
      mm=m                                                              
      mncf=(mm-1)*ncf                                                   
      do 210 l=l1,ntrsy                                                 
      x=0.d0                                                            
      do 220 k=1,ncf                                                    
      if(itk(k,l).eq.0) go to 220                                       
      kk=itk(k,l)                                                       
      ikk=abs(kk)                                                       
      ckk=c(mncf+k)*c(mncf+ikk)                                         
      if(kk.lt.0) ckk=-ckk                                              
      x=x+ckk                                                           
220   continue                                                          
      itm(m,l)=1                                                        
      if(x.lt.0.d0) itm(m,l)=-1                                         
210   continue                                                          
200   continue                                                          
      if(.not.yprt) go to 250                                           
      write(6,*) 'itm'                                                  
      do 230 l=l1,ntrsy                                                 
      write(6,*) 'l=',l                                                 
230   write(6,180) (itm(m,l),m=1,metat)                                 
250   continue                                                          
c                                                                       
c  determination du tableau logique indiquant s il existe ou            
c  non un couplage entre deux etats                                     
c                                                                       
      ind=0                                                             
      do 650 m=2,metat                                                  
      ind=ind+m-1                                                       
      do 670 mm=1,m-1                                                   
      ii=ind+mm                                                         
      yheff(ii)=.true.                                                  
      do 660 j=l1,ntrsy                                                 
      imn=itm(m,j)*itm(mm,j)                                            
660    if(imn.eq.-1) yheff(ii)=.false.                                  
670   continue                                                          
650   continue                                                          
c                                                                       
c  symetrie des etats                                                   
                                                                        
c                                                                       
c  determination du facteur multiplicatif (imul) des groupes de         
c  de determinants degeneres apparaissant dans s.pour le determinant    
c  initial k de chaque groupe,nespo(k) donne l adresse de debut dans    
c  imul caracterisant ce groupe                                         
c                                                                       
c  on ne tient pas compte de la symetrie d espace des excitations       
c  lorsque les generateurs sont invariants d espace                     
c                                                                       
c                                                                       
c  version pert1 : on tient compte des symetries des generateurs        
c  mais pas de la symetrie de spin des excitations lorsque les          
c  generateurs sont invariants de spin                                  
c                                                                       
      do 380 kk=1,ncf                                                   
380   nespo(kk)=0                                                       
      mc=metat*(metat+1)/2                                              
      do 6688 kk=1,nmulz                                                
6688  imul(kk)=0.d0                                                     
      mulmax=nmulz-mc                                                   
      icf=-mc                                                           
      ni=0                                                              
      do 300 k=1,ncf                                                    
      if(ysaut(k)) go to 300                                            
      if(k.eq.1) go to 350                                              
      do 305 i=1,ni                                                     
      ki=nti(i)                                                         
      do 310 l=l1,ntrsy                                                 
      yefk=yef(k,l)                                                     
      yefki=yef(ki,l)                                                   
      yequ=(yefk.and.yefki).or.((.not.yefk).and.(.not.yefki))           
      if(.not.yequ) go to 305                                           
310   continue                                                          
      nespo(k)=nespo(ki)                                                
      go to 300                                                         
305   continue                                                          
350   icf=icf+mc                                                        
      ni=ni+1                                                           
      nti(ni)=k                                                         
      nespo(k)=icf                                                      
      if(icf.gt.mulmax) then                                            
      write(6,*) ' attention attention attention'                       
      write(6,*) ' tableau imul trop petit'                             
      write(6,*) ' modifier parametre nmulz'                            
      stop                                                              
      endif                                                             
      ind=0                                                             
      do 315 m=2,metat                                                  
      ind=ind+m-1                                                       
      do 315 n=1,m-1                                                    
      ix=ind+n                                                          
      im=1                                                              
      if(yprt) write(6,*) 'm,n,ix,im',m,n,ix,im                         
      do 320 l=2,ntrsy                                                  
      if(.not.yef(k,l)) go to 320                                       
      im=im+itm(m,l)*itm(n,l)                                           
320   continue                                                          
      if(yef(k,1).and.yisz0) im=im*(1+itm(m,1)*itm(n,1))                
      if(yprt) write(6,*) 'm,n,ix,im ',m,n,ix,im                        
315   imul(icf+ix)=im                                                   
300   continue                                                          
      if(.not.yprt) go to 365                                           
      icf=-mc                                                           
      do 370 i=1,ni                                                     
      k=nti(i)                                                          
      icf=icf+mc                                                        
      write(6,*) 'i,nti,icf'                                            
      write(6,180) i,nti(i),icf                                         
      write(6,*) 'imul'                                                 
      write(6,180) (imul(icf+ix),ix=1,mc)                               
370   continue                                                          
      write(6,*) 'k,nespo'                                              
      write(6,180) (k,k=1,ncf)                                          
      write(6,180) (nespo(k),k=1,ncf)                                   
365   continue                                                          
      return                                                            
c                                                                       
c  version diexcit : on tient compte des symetries des generateurs      
c  et de la symetrie de spin des excitations meme lorsque les           
c  generateurs sont invariants de spin                                  
c                                                                       
      entry dimul                                                       
      l1=1                                                              
      yisz0=isz.eq.0                                                    
      if(.not.yisz0) l1=2                                               
      do 580 kk=1,ncf                                                   
580   nespo(kk)=0.d0                                                    
      mc=metat*(metat+1)/2                                              
      icf=-mc                                                           
      ni=0                                                              
      do 500 k=1,ncf                                                    
      if(ysaut(k)) go to 500                                            
      if(k.eq.1) go to 550                                              
      do 505 i=1,ni                                                     
      ki=nti(i)                                                         
      do 510 l=l1,ntrsy                                                 
      yefk=yef(k,l)                                                     
      yefki=yef(ki,l)                                                   
      yequ=(yefk.and.yefki).or.((.not.yefk).and.(.not.yefki))           
      if(.not.yequ) go to 505                                           
510   continue                                                          
      nespo(k)=nespo(ki)                                                
      go to 500                                                         
505   continue                                                          
550   icf=icf+mc                                                        
      ni=ni+1                                                           
      nti(ni)=k                                                         
      nespo(k)=icf                                                      
      ind=0                                                             
      do 515 m=2,metat                                                  
      ind=ind+m-1                                                       
      do 515 n=1,m-1                                                    
      ix=ind+n                                                          
      im=1                                                              
      do 520 l=2,ntrsy                                                  
      if(.not.yef(k,l)) go to 520                                       
      im=im+itm(m,l)*itm(n,l)                                           
520   continue                                                          
      if(.not.yef(k,1).or.(.not.yisz0)) go to 515                       
      ims=itm(m,1)*itm(n,1)                                             
      if(ims.eq.-1) im=0                                                
c  doublement eventuel de la contribution de spin dans diexcit          
515   imul(icf+ix)=im                                                   
500   continue                                                          
      if(.not.yprt) go to 565                                           
      write(6,*) ' dimul'                                               
      icf=-mc                                                           
      do 570 i=1,ni                                                     
      k=nti(i)                                                          
      icf=icf+mc                                                        
      write(6,*) 'i ntl icf'                                            
      write(6,180) i,ntl(i),icf                                         
      write(6,*) ' imul'                                                
      write(6,180) (imul(icf+ix),ix=1,mc)                               
570   continue                                                          
      write(6,*) 'k   nespo'                                            
      write(6,180) (k,k=1,ncf)                                          
      write(6,180) (nespo(k),k=1,ncf)                                   
565   continue                                                          
      return                                                            
      entry degen(letat)                                                
c                                                                       
c  rotation des etats degeneres apres diagonalisation                   
c                                                                       
      if(letat.eq.1) return                                             
      l1=1                                                              
      yisz0=isz.eq.0                                                    
      if(.not.yisz0) l1=2                                               
      mmax=letat-1                                                      
      do 400 m=1,mmax                                                   
      if(m.le.jj) go to 400                                             
      jj=m-1                                                            
      mn1=m+1                                                           
      do 435 n=mn1,letat                                                
      jj=jj+1                                                           
      if(dabs(e(n)-e(m)).gt.1.d-6) go to 445                            
435   continue                                                          
445   if(jj.eq.m) go to 400                                             
      if(jj.eq.mn1) go to 420                                           
      write(6,415)                                                      
415   format(1x,/,'0 attention degenerescence multiple',/,              
     &' symmetrisation non prevue dans ce cas.pour les energies',       
     &/,' on peut traiter simultanement tous les etats degeneres ',/,   
     &' et diviser la somme de leurs energies perturbees par la ',/,    
     &' degenerescence.pour les vecteurs modifier degen',//)            
      go to 400                                                         
420   mncf1=m*ncf-ncf                                                   
      mncf2=mncf1+ncf                                                   
416   format(1x,/,'0les etats ',i2,' et ',i2,' sont degeneres',         
     *': symmetrisation par rotation',/)                                
c  rotation des etats                                                   
      do 450 k=1,ncf                                                    
      if(ysaut(k)) go to 450                                            
      nespk=nesp(k)                                                     
      if(nesp(k).eq.0) go to 450                                        
      do 440 ll=1,nespk                                                 
      c11=c(mncf1+k)                                                    
      c12=c(mncf1+k+ll)                                                 
      c1=dabs(c11)                                                      
      c2=dabs(c12)                                                      
      if(dabs(c1-c2).lt.1.d-5) go to 440                                
      c21=c(mncf2+k)                                                    
      c22=c(mncf2+k+ll)                                                 
      cc=c21+c22                                                        
      if(dabs(cc).lt.1.d-8)go to 440                                    
      tgphi=-(c11+c12)/cc                                               
      phi=datan(tgphi)                                                  
      cs=dcos(phi)                                                      
      si=dsin(phi)                                                      
      c1=c11*cs+c21*si                                                  
      c2=c12*cs+c22*si                                                  
      c1=dabs(c1)                                                       
      c2=dabs(c2)                                                       
      c3=dabs(c1-c2)                                                    
      if(c3.ge.1d-4) go to 440                                          
      write(6,416) m,mn1                                                
      do 480 i=1,ncf                                                    
      a=c(mncf1+i)*cs+c(mncf2+i)*si                                     
      c(mncf2+i)=-c(mncf1+i)*si+c(mncf2+i)*cs                           
480   c(mncf1+i)=a                                                      
      go to 400                                                         
440   continue                                                          
450   continue                                                          
400   continue                                                          
410   continue                                                          
      return                                                            
      end                                                               
