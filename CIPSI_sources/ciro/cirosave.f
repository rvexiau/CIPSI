       program ciro                                                    
c-----------------------------------------------------------------------
c      programme de calcul de la matrice densite a l'ordre zero         
c      dans la base moleculaire et la base atomique                     
c      calcul des matrices densites d etat et de transition             
c      toulouse-pise-mexico 1983-1990                                   
c      revisite mai 1991 f spiegelmann                                  
c avant derniere modif : fevrier 92 (romuald) : plus de perturbation    
c                                              mais vecteur+grand       
c    
c version revisite et elargie 17.11.1993 x.gadea, i.paidarova :
c    (changements necessaires aussi dans rozero,rntd et reijkl)
c     
c                                                                       
c      calcul des proprietes apres ic:                                  
c                                                                       
c         moments dipolaires                                            
c         populations de mulliken                                       
c         calcul des orbitales naturelles                               
c                                                                       
c-----------------------------------------------------------------------
c      mode d emploi                                                    
c-----------------------------------------------------------------------
c     namelist roinp                                                    
c                                                                       
c  lecture des etats et des determinants de l'ic                        
c-------------------------------------------------                      
c     files 04 venant de cipsi:                                         
c       n1fic,n2fic numero des files en lecture                         
c                  (defaut n1fic=3, n2fic=4)                            
c                                                                       
c     ymoyen=t, vecteurs de moyen (givens ou davidson) alors:           
c       n1fim,n2fim files 59 (vecteurs de la diagonalisation)           
c                  (defaut ymoyen=f, n1fim=58, n2fim=59)                
c     et                                                                
c       n1fid,n2fid files 60 issues de moyen (determinants)             
c                  (defaut n1fid=61, n2fid=62)                          
c                                                                       
c     yuni=f,deux espaces d'ordre zero(symmetries differentes)          
c            sinon n1fic=n2fic                                          
c                  n1fid=n2fid                                          
c                  n1fim=n2fim                                          
c            (defaut f)                                                 
c                                                                       
c     ietat(i),jetat(j) numeros des premiers et second etats des couples
c           d'etats entre lesquels on calcule les matrices densites     
c                                                                       
c  calcul de la matrice densite                                         
c------------------------------                                         
c                                                                       
c     ityper=0      pas de perturbation (defaut)                        
c            1      moller-plesset barycentrique                        
c            2      epstein nesbet valeur propre                        
c            3      epstein nesbet barycentrique                        
c                                                                       
c     yat=t,matrices densites dans la base d oa(defaut: en base d om)   
c     yprt=t,on imprime les matrices densites en base atomiques         
c         (defaut f)                                                    
c                                                                       
c   proprietes calculees                                                
c-----------------------                                                
c                                                                       
c     yona=t,on fait des orbitales naturelles (defaut f)                
c     ydipol=t,moments dipolaires (defaut f)                            
c     ypopul=t population de mulliken (defaut f)                        
c     yauto=t selection automatique (avec la file 33 ) des
c     etats lab1 (et eventuellement) lab2 (ietat et jetat sont alors 
c     generes automatiquement)
c-----------------------------------------------------------------------
c      files en lecture                                                 
c                                                                       
c        obligatoires                                                   
c                                                                       
c      n1fic et n2fic resultats de cipsi(files 04)                      
c      10: om et s issues de ijkl (file 10 de ijkl)                     
c      40: file des integrales de fok                                   
c      02: om issues de pshf                                            
c                                                                       
c         optionnelles                                                  
c                                                                       
c      si moyen+diag                                                    
c                                                                       
c      n1fid et n2fid vecteurs de givens ou davidson(files 59)          
c      n1fim et n2fim determinants de moyen (files 61)                  
c                                                                       
c                                                                       
c      files en creation                                                
c                                                                       
c      1:  matrices densites perturbees ou non,en base d oa ou d om     
c          (obligatoire)                                                
c     12:  orbitales naturelles comme file 11 de hondo (si yona=t)      
c     18:  moments de transition entre etats electroniques              
c-----------------------------------------------------------------------
c  program compatible avec pshf120/cip120/moy120/bdav/diagbig/          
c-----------------------------------------------------------------------
       implicit real*8(a-h,o-x,z),logical*1(y)                          
      include 'pshf.prm'
      integer*4 nttt
      character*9 day,hour
      character*40 typ_sym1,typ_sym2
      character*3 lab1,lab2
      character*1 charnum
      character*40 titre
      character*40 typdip
      common/uni/yunic,ymsym
      common/polar/calfa(dc),ycv
      common/ic/e(metz2),emp(metz2),fdiag(2*doa),                     
     & fmpb(2*doa),trou(8*ndimh),part(8*ndimh),ne(2*ndimh),nd(2*ndimh), 
     & ispin(2*doa),                                                  
     & iorb(2*doa),itsym(2*doa),its(nsymz,nsymz),igels(2),icount,   
     & norb,noca,nocb,mnorb,ncf,nsym,isym,ntrsy,metat1,metat2,          
     & yoc(2*doa),yprt,ypertu,yuni,ymp                                         
      integer*4 ne,ndcip,iorb,ispin,itsym,its                 
      integer*4 trou,part                 
      common/fimmif/n1fic,n2fic,n1fim,n2fim,ymoyen                  
      common/gel/nao                                                    
       common/cd/r(ncouz*(doa**2)),ndd,ndt,ncf1,ncf2,metat,icf2,ncou, 
     &ibeg(ncouz),ietat(ncouz),jetat(ncouz),ycou(metz,metz), 
     &yiet(metz),yjet(metz)                                             
      integer*4 title(40),numero(21)                          
      dimension ndcip(2*ndimh)
      dimension cmo(doa,doa),rao(doa,doa),rc(doa,doa)       
      dimension ccc(doa*doa)                                        
      datatitle/'  ','  ','et','at','s ',6*'  ','or','dr','e ',26*'  '/ 
      data numero/' 1',' 2',' 3',' 4',' 5',' 6',' 7',' 8',' 9','10',    
     *'11','12','13','14','15','16','17','18','19','29',' 0'/           
c      data charnum/'1','2','3','4','5','6','7','8','9'/
      dimension eb(metz2),nttt(nexz)                                    
      namelist/roinp/yuni,yprt,ypopul,ymoyen,ycv,                       
     & ietat,jetat,titre,ydipol,yat,yona,calfa,ymsym,
     * yauto,lab1,lab2,typ_sym1,typ_sym2,mspin
c      data titre/40*' '/                                               
c      call date(day)
c      call time(hour)
      write(6,9999)day,hour
 9999 format(/,80(1h*),/,8(10h   ciro   ),/, 10x,a10,5x,a10,/,
     * 80(1h*))
      call openf
      yuni=.false.                                                      
      yona=.false.                                                      
      ymoyen=.false.                                                    
      yat=.false.
      ymsym=.true.
      yauto=.false.
      yprt=.false.                                                      
      ityper=0                                                          
      ypertu=.false.                                                    
      ymp=.true.                                                        
      yenb=.false.                                                      
      yenvp=.false.                                                     
      ypopul=.false.                                                    
      ydipol=.false.
      ycv=.true.
      icount=0                                                          
      do 1 i=1,metz                                                     
      yiet(i)=.true.                                                    
      yjet(i)=.true.                                                    
      do 1 j=1,metz                                                     
1     ycou(i,j)=.true.                                                  
      n1fic=3                                                           
      n2fic=4                                                           
      n1fim=58                                                          
      n2fim=59                                                          
      n1fid=61                                                          
      n2fid=62                                                          
      maxm=metz                                                         
      maxm2=maxm+maxm                                                   
      ncf1=0                                                            
      ncf2=0                                                            
      metat1=0                                                          
      metat2=0                                                          
      do 33 i=1,ncouz                                                   
      ietat(i)=0                                                        
33    jetat(i)=0                                                        
      read(5,roinp)                                                     
      if (ityper.ne.0) then
	 write(6,*)'!! -- attention -- !!'
	 write(6,*) 'cette version de ciro ne sait plus faire 
     &	 de perturbation'
      endif
      if(ityper.ne.0)ypertu=.true.                                      
      if(ityper.eq.3) yenb=.true.                                       
      if(ityper.eq.2) yenvp=.true.                                      
      yunic=yuni                                                        
      if(ydipol) yat=.true.                                             
      if(yuni)then                                                      
      n2fic=n1fic                                                       
      n2fim=n1fim                                                       
      n2fid=n1fid
      endif                                                             
      typdip=titre//lab1//lab2
      write(6,950) titre                                                
950   format(///3x,90('*')/3x,'*',88x,'*'/3x,'    calcul de la ',       
     * 'matrice densite perturbee a l ordre zero',26x,'*'/              
     * 3x,'*',4x,a20,4x,'*'/3x,'*',88x,'*'/3x,90('*')///)              
      rewind n1fic                                                      
      if(.not.ymoyen) then                                              
      read(n1fic)norb,noca,nocb,metat1,ncf1,isz,yion,ii,i
     &,(ne(j),ndcip(j),j=1,ncf1),(trou(j),part(j),j=1,ii),(cj,j=1,i)   
c     &,(e(j),j=1,metat1),(emp(j),j=1,metat1),(eb(j),j=1,metat1),        
c     &(jj,j=1,metat1),(x,j=1,ncf1),(fmpb(j),j=1,norb)                   
      do i=1,ncf1
	 nd(i)=ndcip(i)
       enddo
       write(6,*) ' cipsi '                                             
       write(6,*) '   ncf1',ncf1                                        
       else                                                             
       read(n1fic)norb,noca,nocb,metat1,ncf1,isz,yion                   
       write(6,*) ' ciro  1'
       read(n1fim) ncf1,metat1,(cc,i=1,ncf1*metat1),(e(i),i=1,metat1)   
       write(6,*) ' determinants de moyen: ',ncf1                       
       ndf=0                                                            
       do 6300idet=1,ncf1                                               
       nd(idet)=ndf                                                     
       read(n1fid,end=6341)nec                                          
       ne(idet)=nec                                                     
       if(nec.ne.0)read(n1fid,end=6341)(nttt(j),j=1,nec+nec)            
       do 6301j=1,nec                                                   
       ndf=ndf+1                                                        
       trou(ndf)=nttt(j)                                                
6301   part(ndf)=nttt(j+nec)                                            
       ndi=nd(idet)                                                     
cwr       write(6,*) nec,nd(idet),' t ',(trou(ndi+kk),kk=1,nec),        
cwr     2 'p ', (part(ndi+kk),kk=1,nec)                                 
6300   continue                                                         
cwr       write(6,*) ' moyen '                                          
cwr       write(6,*) '   ncf1',ncf1                                     
cwr       write(6,7876) (cc,i=1,ncf1)                                   
      endif                                                             
420   continue                                                          
      write(6,*)' ncf1=',ncf1                                           
      mnorb=norb+norb                                                   
      rewind n2fic                                                      
      if(.not.ymoyen) then                                              
      read(n2fic)norb,noca,nocb,metat2,ncf2,isz,yion,iiii,iii           
     &,(ne(j+ncf1),ndcip(j+ncf1),j=1,ncf2)                              
     &,(trou(j+ii),part(j+ii),j=1,iiii)                                 
     &,(cj,j=1,iii)                                                     
     &,(e(j+metz),j=1,metat2)                                           
     &,(emp(j+metz),j=1,metat2),(eb(j+metz),j=1,metat2),(jj,j=1,metat2) 
     &,(x,j=1,ncf2),(fmpb(j+norb),j=1,norb)                             
c                                                                       
       do i=1,ncf2
	nd(i+ncf1)=ndcip(i+ncf1)
       enddo
       write(6,*) ' cipsi '                                             
       write(6,*) '   ncf2',ncf2                                        
7876   format(10(1x,f10.6))                                             
       else                                                             
       rewind n2fim                                                     
       rewind n2fid                                                     
       read(n2fic)norb,noca,nocb,metat2,ncf2,isz,yion                   
       read(n2fim) ncf2,metat2,(cc,i=1,ncf2*metat2),                    
     2 (e(metz+i),i=1,metat2)                                           
       write(6,*) ' determinants de moyen: ',ncf2                       
          write(6,*) ' ndf ', ndf                                       
       do 6500 idet=1,ncf2                                              
       nd(ncf1+idet)=ndf                                                
       read(n2fid,end=6341)nec                                          
       ne(ncf1+idet)=nec                                                
       if(nec.ne.0)read(n2fid,end=6341)(nttt(j),j=1,nec+nec)            
       do 6501 j=1,nec                                                  
       ndf=ndf+1                                                        
       trou(ndf)=nttt(j)                                                
6501   part(ndf)=nttt(j+nec)                                            
       ndi=nd(idet+ncf1)                                                
cwr       write(6,*) idet,nec,ndi,' t ',(trou(ndi+kk),kk=1,nec),        
cwr     2 'p ', (part(ndi+kk),kk=1,nec)                                 
6500   continue                                                         
       if(ndf.gt.8*ndimh) write(6,*) ' ndf velike', ndf
cwr       write(6,*) ' moyen '                                          
cwr       write(6,*) '   ncf2',ncf2                                     
cwr       write(6,7876) (cc,i=1,ncf2)                                   
      endif                                                             
      write(6,*)' ncf2=',ncf2                                           
      ncf=ncf1+ncf2                                                     
      if(yenb.or.yenvp) go to 24                                        
      do 22 j=1,maxm2                                                   
22    e(j)=emp(j)                                                       
      go to 21                                                          
24    ymp=.false.                                                       
      if(yenvp) go to 21                                                
      do 280 j=1,maxm2                                                  
280   e(j)=eb(j)                                                        
21    continue                                                          
      metat=metat1+metat2                                               
      call reijkl                                                       
      ndd=nao*nao                                                       
      do 440 k=ncf1+1,ncf                                               
440   nd(k)=nd(k)+ii                                                    
      ncou=0                                                            
      do 25 i=1,ncouz                                                   
25    if(ietat(i).ne.0) ncou=ncou+1                                     
      ndt=ndd*ncou                                                      
       do40 i=1,norb                                                    
       ispin(i)=0                                                       
       ispin(i+norb)=1                                                  
       iorb(i)=i                                                        
40     iorb(i+norb)=i                                                   
       iorb(norb+norb+1)=0
       igels(1)=1                                                       
       igels(2)=norb+1                                                  
      if(yauto) call selec(yuni,mspin,typ_sym1,lab1,metat1,
     *				      typ_sym2,lab2,metat2,
     *  				      ietat,jetat,ncou)
      ibegin=-ndd                                                       
      do 26 i=1,ncouz                                                   
      if(ietat(i).eq.0) go to 26                                        
      ibegin=ibegin+ndd                                                 
c      ibeg(ietat(i),jetat(i))=ibegin
      ibeg(i)=ibegin                                    
      yiet(ietat(i))=.false.                                            
      yjet(jetat(i))=.false.                                            
      ycou(ietat(i),jetat(i))=.false.                                   
26    continue                                                          
      do 800 i=1,ndt                                                    
800   r(i)=0.d0                                                         
      if(.not.yuni) go to 830                                           
      do 820 i=1,noca                                                   
      id=(i-1)*norb+i                                                   
      ind=-ndd                                                          
      do 810 icou=1,ncou                                                
      ind=ind+ndd                                                       
810   if(ietat(icou).eq.jetat(icou)) r(ind+id)=2.d0                     
820   continue                                                          
830   continue                                                          
      call rozero                                                       
      if(yprt) then                                                     
      ind=-ndd                                                          
      do 10 icou=1,ncou                                                 
      write(6,777) ietat(icou),jetat(icou)                              
777   format(/,2x,'matrice densite a l ordre zero entre les etats',     
     *i3,' et ',i2)                                                     
      ind=ind+ndd                                                       
10    call scrivi(ind,norb)                                             
      endif                                                             
c                                                                       
      if (yona) then                                                    
      kon=0                                                             
      do 881 i=1,norb                                                   
      do 881 j=i,norb                                                   
      kon=kon+1                                                         
      ij=(j-1)*norb + i                                                 
 881  ccc(kon)=-r(ij)                                                   
      call given(ccc,ccc(nao*nao+1),ccc(nao*nao+nao+1),norb,norb,norb)  
      end if                                                            
                                                                        
                                                                        
      write(6,*) ' nao,norb ',nao,norb                                  
         rewind 10                                                      
         read(10) ((cmo(i,j),i=1,nao),j=1,norb)                         
	 if(nao.ne.norb) then
	 read(11) ((rao(i,j),i=1,nao),j=1,nao)
	 do j=norb+1,nao
	  cmo(i,j)=rao(i,j)
         enddo
	 endif
      write(6,*) ' om'                                                  
      if(yprt) then                                                     
      call escriu(cmo,nao,norb,doa)                                   
      endif
      if(yona) then                                                     
      ncdep=nao*(nao+1)                                                 
      do 882 j=1,norb                                                   
        ccc(nao*nao+j)=-ccc(nao*nao+j)                                  
      do 882 k=1,nao                                                    
      kjna= (j-1)*nao + k                                               
      ccc(kjna)=0.d0                                                    
        do 882 i=1,norb                                                 
        ijnm=ncdep + (j-1)*norb + i                                     
        ccc(kjna)=ccc(kjna) + ccc(ijnm)*cmo(k,i)                        
 882  continue                                                          
      write(6,*) '   occupation des orbitales naturelles '              
      write(6,883)   (ccc(i),i=1+nao*nao,norb+nao*nao)                  
 883  format ((2x,5(f10.6,2x)))                                         
      write(6,*)                                                        
      write(6,*) '   orbitales naturelles '                             
      kon=norb/5                                                        
      if (kon*5.ne.norb) kon=kon+1                                      
      do 884 k=1,kon                                                    
      maxk=k*5                                                          
      maxk=min(norb,maxk)                                               
       do 885 i=1,5                                                     
 885   trou(i)=(k-1)*5 + i                                              
      write(6,886)(trou(i),i=1,5)                                       
 886  format (/,6x,5(i2,10x),/)                                         
      do 884 i=1,nao                                                    
      write(6,883) (ccc((j-1)*nao + i),j=1+(k-1)*5,maxk)                
 884  continue                                                          
      write(12) (ccc(i),i=1,nao*norb),(ccc(i),i=nao*nao+1,nao*nao+norb) 
      end if                                                            
c                                                                       
c  matrice densite dans la base atomique
c
      title(16)=numero(21)                                              
c      ipa=1
c      if(ipa.eq.1) go to 201
      do 200 icou=1,ncou                                                
      m1=ietat(icou)                                                    
      m2=jetat(icou)                                                    
      title(7)=numero(m1)                                               
      title(9)=numero(m2)                                               
      do 230 j=1,nao                                                    
      ind=ibeg(icou)                                                   
      do 230 k=1,norb                                                   
      rckj=0.d0                                                         
      do 220 l=1,norb                                                   
      ind=ind+1                                                         
220   rckj=rckj+r(ind)*cmo(j,l)                                         
230   rc(k,j)=rckj                                                      
c  calcul de la trace de la matrice densite dans la base moleculaire
      traro=0.d0
      do io=1,norb
      traro=traro+r(ibeg(icou)+1+(io-1)*(norb+1))
      enddo
      write(6,*)' trace de ro mol pour ietat jetat ',m1,m2,traro

      if(yuni.and.m1.eq.m2) go to 100                                   
      do 250 i=1,nao                                                    
      do 250 j=1,nao                                                    
      rij=0.d0                                                          
      do 240 k=1,norb                                                   
240   rij=rij+cmo(i,k)*rc(k,j)                                          
250   rao(i,j)=rij                                                      
      go to 190                                                         
100   continue                                                          
      do 150 i=1,nao                                                    
      do 150 j=1,i                                                      
      rij=0.d0                                                          
      do 140 k=1,norb                                                   
140   rij=rij+cmo(i,k)*rc(k,j)                                          
      rao(j,i)=rij                                                      
150   rao(i,j)=rij                                                      
190   continue                                                          
      if(yprt) then                                                     
      write(6,260) title                                                
      call escriu(rao,nao,nao,doa)                                    
      endif                                                             
260   format(//5x,'matrice densite en base atomique'/40a2/)             
      ind=ibeg(icou)                                                   
      do 180 j=1,nao                                                    
      do 180 i=1,nao                                                    
      ind=ind+1                                                         
180   r(ind)=rao(i,j)                                                   
      write(6,*)' trace de ro atomique pour ietat jetat ',m1,m2,traro
      traro=0.d0
      do io=1,norb
      traro=traro+r(ibeg(icou)+1+(io-1)*(norb+1))
      enddo
      write(6,*)' trace de ro atomique pour ietat jetat ',m1,m2,traro

200   continue
201   continue
c                                                                       
c  calcul des proprietes                                                
c                                                                       
      call prop(ydipol,ypopul,typdip)                                   
      close (1)                                                         
c      call date(day)
c      call time(hour)
      write(6,9998)day,hour
 9998 format(/,80(1h*),/,8(10h fin ciro ),/, 10x,a10,5x,a10,/,
     * 80(1h*))
      stop                                                              
6341  write(6,*) ' erreur lecture sur file moyen ', n1fim               
6342  write(6,*) ' erreur lecture sur file moyen ', n2fim               
      end                                                               
