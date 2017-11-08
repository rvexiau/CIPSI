c       Ecriture dans les fichiers de sortie desactive 
      subroutine pert1(kkmod,nf)                                        
      implicit real*8 (a-h,o-x,z),logical*1 (y)                         
      include 'pshf.prm'
      integer*4 ne,nd,trou,part,iorb,iwspin,itsym,its,isytr,nesp,nespo  
      integer*4 imul                                                    
      character*4 ibla2                                                 
      character*8 deg                                                   
      common/spo/fdiag(2*doa),hdiag(ndetz),tocom(2*doa),coeff,      
     1ndef,ispin(2*doa),iwspin(2*doa),numero(metz),          
     2ne(ndetz+1),nd(ndetz+1),trou(ntpz),part(ntpz),iorb(2*doa),      
     3itsym(2*doa),its(nsymz,nsymz),yoc(2*doa),yocs(2*doa),       
     4ygen(ndetz),yocd(2*doa,ndetz),yspin                                    
      common/degene/nesp(ndetz),nespo(ndetz),isytr(ndetz,nsymz),kkbrd,  
     1ysaut(ndetz),imul(nmulz),yheff(metz*(metz+1)/2),yef(ndetz,nsymz)  
      common/cipsi/escf,test,teseff,tdet,twdet,ndiag,norb,mnorb,
     1             noca,nocb,noc1,noc2,igela,igelb,ncf,nsym,isym,isz,
     2             ntrsy,metat,lecdet,maxdet,yprt,yion,ypertu,ybrd,yen,
     3             ywvi,ywvf,ytoul,ystkic,ywdeg
      common/ic/c(metz*ndetz),e(ndetz),eb(metz),e2b(metz),emp(metz),    
     1teste(metz),e2mp(metz),dpivot,ietat,ishif,e2vp(metz),   
     2fmpb(2*doa),hefmp(nhefz),hefvp(nhefz),hefb(nhefz),              
     3psimp(metz),psivp(metz)         
      common/hist/nom(10),nam(10,metz),aval(10),tabmp(10,metz),         
     1tau,semp(metz),sevp(metz),tabvp(10,metz),nclass,ymoyen  
      common/orper/yreduc,yorper(doa)
      common/tome/taum(metz),yselm,ysele                                
      dimension stok(ndetz),nts(8),nps(8),vil(metz),venvp(metz)                     
      dimension ikk(ndetz)                                              
      dimension him(metz),dmp(metz),dvp(metz),db(metz)                  
      data ibla1,ibla2/0,'=>'/,deg/'degenr'/                            
c     test sur la nullite
      if(kkmod.lt.0) kkmod=0
      kk=kkmod
      nec=ne(nf)                                                        
      nec1=nec-1                                                        
      yf=nec.eq.0                                                       
      ndf=nd(nf)                                                        
      if(yf) go to 17                                                   
      do 18 i=1,nec                                                     
         nt=trou(ndf+i)    
         nts(i)=nt                                                         
         if(yoc(nt)) go to 200                                             
         yoc(nt)=.true.                                                    
         np=part(ndf+i)                                                    
         nps(i)=np                                                         
         if(yoc(np)) go to 200                                             
 18   yoc(np)=.true.                                                    
 17   continue                                                          
c
c    exclusion des orbitales non actives
c
      if(yreduc) then
       if(nec.eq.0) go to 200
       do i=1,nec
	nt=iorb(trou(ndf+i))
	np=iorb(part(ndf+i))
	if(yorper(np)) go to 717
        enddo
	go to 200
717     continue
	endif
c                                                                       
c   test de non repetition et stockage utile                            
c                                                                       
      kk=kkmod                                                          
      mesp=1                                                            
c
c     si |i> interagit avec  l'etat |j> de S  pour  j < kk      
c     alors |i> a deja ete cree                                         
c
      j=1                                                               
30    if(j.ge.kk) goto 31                                               
      if(.not.ygen(j)) go to 26                                         
      mesp=mesp+nesp(j)                                                 
      ndif=nec-ne(j)                                                    
      if(ndif.gt.2) go to 26                                            
      if(ne(j).eq.0) go to 200                                          
      nej=ne(j)+nd(j)                                                   
      ndj=nd(j)+1                                                       
       do 25 i=ndj,nej                                                  
      if(.not.yoc(trou(i)))ndif=ndif+1                                  
25    if(.not.yoc(part(i)))ndif=ndif+1                                  
      if(ndif.le.2) go to 200                                           
26    j=j+1                                                             
      goto 30                                                           
31    ijune=0                                                           
      if(kk.eq.0) kk=1                                                  
      lesp=0                                                            
      nespf=nesp(kk)                                                    
      if (kkmod.eq.0) nespf=0                                           
      kkmax=kk+nespf                                                    
c                                                                       
c     si |i> efficace, calcul et stockage des hji dans stock                             
c     j appartient a S, j=1,ncf
c                                                                       
      do 40 j=kk,ncf
      ndif=nec-ne(j)                                                    
      if(ndif.gt.2) go to 40                                            
      if(ne(j).eq.0) go to 34                                           
      nej=ne(j)+nd(j)                                                   
      ndj=nd(j)+1                                                       
       do 32 i=ndj,nej                                                  
      if(.not.yoc(trou(i)))ndif=ndif+1                                  
32    if(.not.yoc(part(i)))ndif=ndif+1                                  
      if(ndif.gt.2) go to 40                                            
34    if(ndif.eq.0) go to 200                                           
      if(.not.ysaut(j)) go to 35                                        
      if(j.le.kkmax) lesp=lesp+1                                        
35    ijune=ijune+1                                                     
      ikk(ijune)=j                                                      
40    continue                                                          
      if(ijune.eq.0) go to 200                                          
      do 45 j=1,ijune                                                   
      k=ikk(j)                                                          
45    stok(j)=hntd(nf,k)                                                
c                                                                       
c   facteurs multiplicatifs tenant compte des equivalences              
c                                                                       
      fmul=0.                                                           
      nmul=0                                                            
      if(nespf.eq.0) go to 65                                           
      ybis=.true.                                                       
      fmulbr=1.d0/dfloat(lesp+1)                                        
      fmul=fmulbr*dfloat(nespf+1)-1.d0                                  
      nmul=nmul+idint(fmul)                                             
c                                                                       
c   calcul et stockage des hmi dans him                                 
c                                                                       
      go to 66                                                          
65    ybis=.false.                                                      
66    continue                                                          
      yet=.false.                                                       
      mncf = -ncf                                                       
      xam=0.d0                                                          
      ydeg=.false.                                                      
      hmpi=hmp(nf)                                                      
      if(yen) hvpi=hdig(nf)                                             
      do 60 m=1 ,metat                                                  
      bcalc=0.0                                                         
      mncf = mncf + ncf                                                 
c**   mncf=ncf*(m-1)                                                    
       do 61 j=1,ijune                                                  
      k=ikk(j)                                                          
      km=mncf+k                                                         
61    bcalc=bcalc+(stok(j)*c(km))                                       
      him(m)=bcalc                                                      
c                                                                       
c  partition moller plesset                                             
c                                                                       
      dmp(m)=emp(m)-hmpi                                                
      cimp=bcalc/dmp(m)                                                 
      car=cimp*cimp                                                     
      psimp(m)=psimp(m)+car                                             
      eimp=bcalc*cimp                                                   
      semp(m)=eimp                                                      
      e2mp(m)=e2mp(m)+eimp                                              
      xim=dabs(cimp)                                                    
      if(ysele)xim=dabs(eimp)                                           
      deno=dmp(m)                                                       
      eim=eimp                                                          
      cim=cimp                                                          
      hii=hmpi                                                          
c                                                                       
c   partitions epstein-nesbet valeur propre et barycentrique            
c                                                                       
      if(yen) then                                                      
c                                                                       
c   valeur propre                                                       
c                                                                       
         dvp(m)=e(m)-hvpi                                               
         civp=bcalc/(dvp(m))                                            
         carvp=civp*civp                                                
         psivp(m)=psivp(m)+carvp                                        
         eivp=bcalc*civp                                                
         sevp(m)=eivp                                                   
         e2vp(m)=e2vp(m)+eivp                                           
         xim=dabs(civp)                                                 
         if(ysele)xim=dabs(eivp)                                        
         deno=dvp(m)                                                    
         eim=eivp                                                       
         cim=civp                                                       
         hii=hvpi                                                       
c                                                                       
c   barycentrique                                                       
c                                                                       
         db(m)=eb(m)-hvpi                                               
         cib=bcalc/(db(m))                                              
         eib=bcalc*cib                                                  
         e2b(m)=e2b(m)+eib                                              
       endif                                                            
c                                                                       
c  test degenerescence precoce                                          
c                                                                       
        if(dabs(deno).gt.0.05) go to 80                                 
        if(dabs(eim).lt.1.d-6) go to 80                                 
c les degeneres sont inclus dans la selection                           
      if(ywdeg)ydeg=.true.                                              
      nex=ne(nf)                                                        
      mono1=nd(nf)                                                      
      mono2=mono1+nex                                                   
      mono1=mono1+1                                                     
      if(.not.yet) then                                                 
c      write(6,2810) numero(m),deg,eim,bcalc,hii,                        
c     * (iorb(trou(i)),iwspin(trou(i)),i=mono1,mono2),                   
c     *ibla1,ibla2,(iorb(part(i)),iwspin(part(i)),i=mono1,mono2)         
      else                                                              
c      write(6,2811) numero(m),deg,eim,bcalc,hii                         
      end if                                                            
 2810 format(1x,i2,2x,a6,2x,f9.6,2x,f8.5,2x,f6.3,11x,13(i3,a2,' '))     
 2811 format(1x,i2,2x,a6,2x,f9.6,2x,f8.5,2x,f6.3)                       
c                                                                       
c                                                                       
c   tests selection                                                     
c                                                                       
c   selection ecriture listing                                          
c                                                                       
80    continue                                                          
      if(xim.lt.teste(m)) go to 67                                      
      if(.not.yet) then                                                 
      nex=ne(nf)                                                        
      mono1=nd(nf)                                                      
      mono2=mono1+nex                                                   
      mono1=mono1+1                                                     
c       write(6,1060) numero(m),cim,eim,bcalc,hii,                       
c     * (iorb(trou(i)),iwspin(trou(i)),i=mono1,mono2),                   
c     * ibla1,ibla2, (iorb(part(i)),iwspin(part(i)),i=mono1,mono2)       
       yet=.true.                                                       
      else                                                              
c       write(6,1061) numero(m),cim,eim,bcalc,hii                        
      end if                                                            
1060  format(1x,i2,2x,f6.3,2x,f9.6,2x,f6.3,2x,f9.6,2x,10(i3,a2,' '))    
1061  format(1x,i2,2x,f6.3,2x,f9.6,2x,f6.3,2x,f9.6)                     
      go to 68                                                          
67    if(dpivot.lt.xim) dpivot=xim                                      
   68 continue                                                          
c                                                                       
c  selection file 60                                                    
c                                                                       
      if(.not.ymoyen) go to 3267                                        
      xim=xim*taum(m)                                                   
      vil(m)=xim                                                        
      xam=dmax1(xam,xim)               
      venvp(m)=eivp*(fmul+1.)     
 3267 continue                                                          
c                                                                       
c  calcul des contributions a l energie                                 
c                                                                       
      if(.not.ybis) go to 70                                            
        e2mp(m)=e2mp(m)+eimp*fmul                                       
        psimp(m)=psimp(m)+car*fmul                                      
      if(yen) then                                                      
         e2vp(m)=e2vp(m)+eivp*fmul          
         e2b(m)=e2b(m)+eib*fmul                                         
         psivp(m)=psivp(m)+carvp*fmul                                   
      endif                                                             
70    continue                                                          
   60 continue                                                          
c                                                                       
c  histogramme et ecriture file 60                                      
c                                                                       
      call isto(xam,vil,nmul,fmul,metat,teseff)                         
      if(ydeg) xam=1.d0                                                 
      if((xam.gt.tau).and.ymoyen)call pie(nec,nts,nps,xam,venvp,metat)              
c                                                                       
c  contributions a heff                                                 
c                                                                       
      if(.not.ybrd.or.metat.le.1) go to 171                             
      ind=0                                                             
      do 170 m=2,metat                                                  
      ind=ind+m-1                                                       
      do 170 mm=1,m-1                                                   
      mmm=ind+mm                                                        
      if(.not.yheff(mmm)) go to 170                                     
      rmul=imul(kkbrd+mmm)                                              
      facmul=rmul*fmulbr                                                
      hnum=him(m)*him(mm)                                               
      hhmp=hnum/(dmp(m)+dmp(mm))                                        
      hhmp=hhmp+hhmp                                                    
      if(ybis) hhmp=hhmp*facmul                                         
      hefmp(mmm)=hefmp(mmm)+hhmp                                        
        if(yen) then                                                    
          hhvp=hnum/(dvp(m)+dvp(mm))                                    
          hhb=hnum/(db(m)+db(mm))                                       
          hhvp=hhvp+hhvp                                                
          hhb=hhb+hhb                                                   
          if(ybis) then                                                 
             hhvp=hhvp*facmul                                           
             hhb=hhb*facmul                                             
          endif                                                         
         hefvp(mmm)=hefvp(mmm)+hhvp                                     
         hefb(mmm)=hefb(mmm)+hhb                                        
        endif                                                           
170   continue                                                          
171   continue                                                          
      ietat=ietat+nmul+1                                                
200   continue
      do  j=1,nec                                                    
         itt=trou(ndf+j)
         ipp=part(ndf+j)
         yoc(itt)=.false.                                          
         yoc(ipp)=.false.                                          
      enddo
      if(nf.le.ncf) go to 503                                           
      if(ne(kk).lt.nec)  then
          return                                          
      endif
      nec=ne(kk)                                                        
      ndk=nd(kk)                                                        
      if(yf) go to 503                                                  
  501 if(ne(nf).gt.nec) return                                          
      do 220 i=1,nec                                                    
      trou(ndf+i)=trou(ndk+i)                                           
220   part(ndf+i)=part(ndk+i)                                           
  503 continue                                                          
      return                                                            
      end                                                               
