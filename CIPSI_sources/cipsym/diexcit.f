c       Ecriture dans les fichiers de sortie desactive 
      subroutine diexcit                                                
      implicit real*8 (a-h,o-x,z),logical*1 (y)                         
      include 'pshf.prm'
      common/spo/fdiag(2*doa),hdiag(ndetz),tocom(2*doa),coeff,      
     1ndef,ispin(2*doa),iwspin(2*doa),numero(metz),          
     2ne(ndetz+1),nd(ndetz+1),trou(ntpz),part(ntpz),iorb(2*doa),      
     3itsym(2*doa),its(nsymz,nsymz),yoc(2*doa),yocs(2*doa),       
     4ygen(ndetz),yocd(2*doa,ndetz),yspin                                    
      common/degene/nesp(ndetz),nespo(ndetz),isytr(ndetz,nsymz),kkbrd,  
     1ysaut(ndetz),imul(nmulz),yheff(metz*(metz+1)/2),yef(ndetz,nsymz)  
      integer*4 ne,nd,trou,part,iorb,iwspin,itsym,its,isytr,nesp,nespo  
      integer*4 imul                                                    
      character*2 ibla2                                                 
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
      common/tome/taum(metz),yselm,ysele                                
      common/orper/yreduc,yorper(doa)
      dimension isp(4,3),stok(ndetz),ikk(ndetz),                        
     1 irapid(ndetz),yefsp(ndetz),                                      
     2 nts(8),nps(8),vil(metz),venvp(metz)                                          
      dimension him(metz),dmp(metz),dvp(metz),db(metz)                  
      character*8 deg                                                   
      integer*4 isptr(ndetz)                                            
      equivalence (yefsp(1),yef(1,1)),(isptr(1),isytr(1,1))             
      data ibla2/'=>'/,deg/'degenr'/                                    
      data isp/0,1,0,1,4*0,0,1,1,0/                                     
      data ibla1/0/                                                     
      nf=ncf+1                                                          
      write(6,*) 'diexcit yreduc=', yreduc
      if(ybrd) call dimul                                               
      ndf=nd(nf)                                                        
      n=norb                                                            
c                                                                       
c   generation des diexcitations sur s                                  
c                                                                       
      write(6,6)                                                        
6     format(' diexcitations')                                          
      igel=igela                                                        
      nbor=noc2                                                         
      nbor1=noc1+1                                                      
      nm1=n-1                                                           
      iva=1                                                             
 2002 continue                                                          
      nbi=nbor-1                                                        
      yij=.true.                                                        
      ykl=.true.                                                        
      yijkl=.true.                                                      
      nspin=1                                                           
      ii=igel-1                                                         
 10   ii=ii+1                                                           
      if(ii.gt.nbor) goto 13                                            
c     write(6,*) igel,nbor,ii                                           
      iis=itsym(ii)                                                     
c     write(6,*) iis,its(iis,iis)                                       
      if(its(iis,iis).ne.1) go to 10                                    
      jj=ii                                                             
      kk=nbor1-1                                                        
 11   kk=kk+1                                                           
      if(kk.gt.n) goto 14                                               
c     write(6,*) nbor1,n,kk                                             
      kks=itsym(kk)                                                     
c     write(6,*) kks,its(kks,kks)                                       
      if(its(kks,kks).ne.1) go to 11                                    
      ll=kk                                                             
      goto 100                                                          
12    continue                                                          
c     write(6,*) ii,kk,kk,ii,ii                                         
      goto 11                                                           
 14   continue                                                          
      goto 10                                                           
 13   continue                                                          
      iva=2                                                             
c     write(6,*) ii,kk,ii,kk,ii,kk                                      
      ykl=.false.                                                       
      yijkl=.false.                                                     
      if(nm1.eq.nbor) go to 2230                                        
      ii=igel-1                                                         
 20   ii=ii+1                                                           
      if(ii.gt.nbor) goto 24                                            
      iis=itsym(ii)                                                     
      if(its(iis,iis).ne.1) go to 20                                    
      jj=ii                                                             
      kk=nbor1-1                                                        
 21   kk=kk+1                                                           
      if(kk.gt.nm1) goto 25                                             
      kks=itsym(kk)                                                     
      kk1=kk+1                                                          
      ll=kk1-1                                                          
 22   ll=ll+1                                                           
      if(ll.gt.n) goto 26                                               
      lls=itsym(ll)                                                     
      kls=its(kks,lls)                                                  
      if(kls.eq.1) goto 100                                             
23    continue                                                          
      goto 22                                                           
 26   continue                                                          
      goto 21                                                           
 25   continue                                                          
      goto 20                                                           
 24   continue                                                          
 2230 continue                                                          
      if(nbi.le.0) go to 5000                                           
      iva=3                                                             
      yij=.false.                                                       
      ykl=.true.                                                        
      ii=igel-1                                                         
 30   ii=ii+1                                                           
      if(ii.gt.nbi) goto 34                                             
      iis=itsym(ii)                                                     
      ii1=ii+1                                                          
      jj=ii1-1                                                          
 31   jj=jj+1                                                           
      if(jj.gt.nbor) goto 35                                            
      jjs=itsym(jj)                                                     
      ijs=its(iis,jjs)                                                  
      if(ijs.ne.1) goto 31                                              
      kk=nbor1-1                                                        
 32   kk=kk+1                                                           
      if(kk.gt.n) goto 36                                               
      kks=itsym(kk)                                                     
      if(its(kks,kks).ne.1) go to 32                                    
      ll=kk                                                             
      goto 100                                                          
 33   continue                                                          
      goto 32                                                           
 36   continue                                                          
      goto 31                                                           
 35   continue                                                          
      goto 30                                                           
 34   continue                                                          
      iva=4                                                             
      if(nm1.eq.nbor) go to 2240                                        
      ykl=.false.                                                       
      nspin=3                                                           
      ii=igel-1                                                         
 40   ii=ii+1                                                           
      if(ii.gt.nbi) goto 45                                             
      iis=itsym(ii)                                                     
      ii1=ii+1                                                          
      jj=ii1-1                                                          
 41   jj=jj+1                                                           
      if(jj.gt.nbor) goto 46                                            
      jjs=itsym(jj)                                                     
      ijs=its(iis,jjs)                                                  
      kk=nbor1-1                                                        
 42   kk=kk+1                                                           
      if(kk.gt.nm1) goto 47                                             
      kks=itsym(kk)                                                     
      kk1=kk+1                                                          
      ll=kk1-1                                                          
 43   ll=ll+1                                                           
      if(ll.gt.n) goto 48                                               
      lls=itsym(ll)                                                     
      kls=its(kks,lls)                                                  
      ijkls=its(ijs,kls)                                                
      if(ijkls.eq.1)  goto 100                                          
44    continue                                                          
      goto 43                                                           
 48   continue                                                          
      goto 42                                                           
 47   continue                                                          
      goto 41                                                           
 46   continue                                                          
      goto 40                                                           
 45   continue                                                          
 2240 continue                                                          
c                                                                       
      goto 5000                                                         
c                                                                       
100   continue    
      ysyspi=.not.yijkl                                                 
      haj=ai(ii,kk,jj,ll)                                               
      if(nspin.gt.1) hak=ai(ii,ll,jj,kk)                                
      do 2000 isa=1,nspin                                               
      iso=ii+isp(1,isa)*norb                                            
      jso=jj+isp(2,isa)*norb                                            
      kso=kk+isp(3,isa)*norb                                            
      lso=ll+isp(4,isa)*norb                                            
      if(.not.yorper(iorb(iso)).or.
     * .not.yorper(iorb(jso)).or.
     * .not.yorper(iorb(kso)).or.
     * .not.yorper(iorb(lso)))go to 1005
      hijkl=0.d0                                                        
      if(isa.le.2) hijkl=haj                                            
      if(isa.ge.2) hijkl=hijkl-hak                                      
      yspi=.false.                                                      
1800   continue                                                         
      yoc(iso)=.true.                                                   
      yoc(jso)=.true.                                                   
      yoc(kso)=.true.                                                   
      yoc(lso)=.true.                                                   
      ndf1=ndf+1                                                        
      ndf2=ndf+2                                                        
      part(ndf1)=kso                                                    
      part(ndf2)=lso                                                    
      trou(ndf2)=jso                                                    
      trou(ndf1)=iso                                                    
      do 230 k=1,ncf                                                    
      i=0                                                               
      if(yocd(iso,k))i=i+1                                              
      if(yocd(jso,k))i=i+1                                              
      if(yocd(kso,k))i=i+1                                              
      if(yocd(lso,k))i=i+1                                              
230   irapid(k)=i                                                       
      yquick=.not.(yocs(iso).or.yocs(jso).or.yocs(kso).or.yocs(lso))    
      if(yquick.and.hijkl.eq.0.d0) go to 2000                           
      mesp=1                                                            
      do 1000 kkk=1,ncf                                                 
      if(.not.ygen(kkk)) go to 1000                                     
      nec =ne(kkk)                                                      
      yf=nec.eq.0                                                       
      if(ysaut(kkk)) go to 1000                                         
      nespf=nesp(kkk)                                                   
      ybis=nespf.gt.0                                                  0
      if(yspi.and..not.yefsp(kkk)) go to 1000                           
c                                                                       
      kkkmax=kkk+nespf                                                  
      lesp=0                                                            
      ndk=nd(kkk)                                                       
 503  nec1=nec+1                                                        
      nec2=nec+2                                                        
      ne(nf)=nec2                                                       
       ijune=1                                                          
      ikk(1)=kkk                                                        
      stok(1)=hijkl                                                     
      if(yquick) go to 310                                              
250   nef=nec+2                                                         
      if(nec.eq.0) go to 260                                            
      if(yocd(iso,kkk).or.yocd(jso,kkk).or.yocd(kso,kkk).or.yocd(lso,   
     *kkk)) go to 1000                                                  
c                                                                       
c
c   tests de non-repetition                                             
c                                                                       
c     si !i> interagit avec  l'etat !j> de -s-  pour  j )= (kk-1 )      
c     alors !i> a deja ete cree                                         
260   j=1                                                               
 300  if(j.ge.kkk) goto 310                                             
      if(.not.ygen(j)) go to 302                                        
      if(irapid(j).eq.0) go to 302                                      
      ndif=ne(j)-nec+2-irapid(j)                                        
      if(ndif.ge.2) go to 302                                           
      if(yf) go to 1000                                                 
      do 301 l=1,nec                                                    
      if(.not.yocd(trou(ndk+l),j))ndif=ndif+1                           
301   if(.not.yocd(part(ndk+l),j))ndif=ndif+1                           
      if(ndif.gt.2) go to 302                                           
      go to 1000                                                        
  302 j=j+1                                                             
      goto 300                                                          
c                                                                       
c     calcul et stockage de hki dans le tableau stok                    
c                                                                       
 310  if(yf) go to 311                                                  
      do 308 i=1,nec                                                    
      part(ndf2+i)=part(ndk+i)                                          
 308  trou(ndf2+i)=trou(ndk+i)                                          
  311 j=kkk+1                                                           
      if(yquick) go to 450                                              
401    if(j.gt.ncf) go to 450                                           
c                                                                       
      if(irapid(j).eq.0) go to 400                                      
      ndif=ne(j)-nec+2-irapid(j)                                        
      if(ndif.gt.2) go to 400                                           
      if(yf) go to 342                                                  
      do 340 l=1,nec                                                    
      if(.not.yocd(trou(ndk+l),j))ndif=ndif+1                           
340   if(.not.yocd(part(ndk+l),j))ndif=ndif+1                           
      if(ndif.gt.2) go to 400                                           
342   if(ndif.eq.0) go to 900                                           
      if(j.le.kkkmax)lesp=lesp+1                                        
640   continue                                                          
      ijune=ijune+1                                                     
      ikk(ijune)=j                                                      
  400 j=j+1                                                             
      goto 401                                                          
450   continue                                                          
      i=0                                                               
      do 455 m=1,ijune                                                  
      j=ikk(m)                                                          
      if(m.ne.1) go to 453                                              
      bcalc=stok(1)                                                     
      go to 454                                                         
453   bcalc=hntd(nf,j)                                                  
 454  if(bcalc.eq.0.d0) go to 455                                       
      i=i+1                                                             
      stok(i)=bcalc                                                     
      ikk(i)=j                                                          
455   continue                                                          
      ijune=i                                                           
      if(ijune.eq.0) go to 900                                          
c                                                                       
c  facteurs multiplicatifs tenant compte des equivalences               
c                                                                       
      yefsk=yefsp(kkk)                                                  
      ysbrd=yefsk.or.ysyspi                                             
      if(.not.ysbrd) go to 457                                          
      ybis=.true.                                                       
      if(.not.yefsk) nespf=nespf+nespf+1                                
457   continue                                                          
      if(isz.ne.0) ysbrd=.false.                                        
      fmul=0.                                                           
      nmul=0                                                            
      kkbrd=nespo(kkk)                                                  
      fmulbr=1.d0/dfloat(lesp+1)                                        
      fmul=fmulbr*dfloat(nespf+1)-1.d0                                  
      nmul=nmul+idint(fmul)                                             
c                                                                       
c   calcul et stockage de hmi dans le tyableau him                      
c                                                                       
650   continue                                                          
      mncf = -ncf                                                       
      yet=.false.                                                       
      xam=0.d0                                                          
      ydeg=.false.                                                      
      hmpi=hmp(nf)                                                      
      if(yen) hvpi=hdig(nf)                                             
      do 460 m=1,metat                                                  
      bcalc=0.0                                                         
      mncf = mncf + ncf                                                 
      do 61 j=1,ijune                                                   
      km=mncf+ikk(j)                                                    
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
c   selection et ecriture listing                                       
c                                                                       
80    continue                                                          
                                                                        
      if(xim.lt.teste(m)) go to 67                                      
      if(.not.yet) then                                                 
      yet=.true.                                                        
      if(yquick) then                                                   
      if(nec.ne.0) then                                                 
c       write(6,1060) numero(m),cim,eim,bcalc,hii,                       
c     * ii,iwspin(iso),jj,iwspin(jso),(iorb(trou(i+ndk)),                
c     *iwspin(trou(i+ndk)),i=1,nec),ibla1,ibla2,                         
c     *kk,iwspin(kso),ll,iwspin(lso),(iorb(part(i+ndk)),                 
c     *iwspin(part(i+ndk)),i=1,nec)                                      
      else                                                              
c       write(6,1060) numero(m),cim,eim,bcalc,hii,                       
c     * ii,iwspin(iso),jj,iwspin(jso),ibla1,ibla2,kk,                    
c     *iwspin(kso),ll,iwspin(lso)                                        
      endif                                                             
      else                                                              
      nex=ne(nf)                                                        
       mono1=ndf                                                        
      mono2=mono1+nex                                                   
      mono1=mono1+1                                                     
c       write(6,1060) numero(m),cim,eim,bcalc,hii,                       
c     * (iorb(trou(i)),iwspin(trou(i)),i=mono1,mono2),                   
c     * ibla1,ibla2, (iorb(part(i)),iwspin(part(i)),i=mono1,mono2)       
      end if                                                            
      else                                                              
c       write(6,1061) numero(m),cim,eim,bcalc,hii                        
      endif                                                             
1060  format(1x,i2,2x,f6.3,2x,f9.6,2x,f6.3,2x,f9.6,2x,10(i3,a2,' '))    
1061  format(1x,i2,2x,f6.3,2x,f9.6,2x,f6.3,2x,f9.6)                     
      go to 68                                                          
67    if(dpivot.lt.xim) dpivot=xim                                      
68    continue                                                          
c                                                                       
c   selection file 60                                                   
c                                                                       
                                                                        
      if(.not.ymoyen) go to 2356                                        
      xim=xim*taum(m)                                                   
      vil(m)=xim              
      xam=dmax1(xam,xim)     
      venvp(m)=eivp*(fmul+1.)   
2356  continue                                                          
c                                                                       
c   contributions a l energie                                           
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
460   continue                                                          
      ietat=ietat+nmul+1                                                
c                                                                       
c   histograme et ecriture file 60                                      
c                                                                       
      if(.not.ymoyen) go to 2358                                        
      call isto(xam,vil,nmul,fmul,metat,teseff)                         
      nex=ne(nf)                                                        
      mono1=ndf                                                         
      mono2=mono1+nex                                                   
      mono1=mono1+1                                                     
      istf=0                                                            
      do 1 ji=mono1,mono2                                               
      istf=istf+1                                                       
      nts(istf)=trou(ji)                                                
      nps(istf)=part(ji)                                                
 1    continue                                                          
      if(ydeg) xam=1.d0        
      if((xam.gt.tau).and.ymoyen) call pie(nex,nts,nps,xam,venvp,metat)             
2358  continue                                                          
c                                                                       
c   contributions a heff                                                
c                                                                       
      if(metat.eq.1.or..not.ybrd) go to 900                             
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
      if(ysbrd) hhmp=hhmp+hhmp                                          
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
          if(ysbrd) then                                                
              hhvp=hhvp+hhvp                                            
              hhb=hhb+hhb                                               
          endif                                                         
         hefvp(mmm)=hefvp(mmm)+hhvp                                     
         hefb(mmm)=hefb(mmm)+hhb                                        
        endif                                                           
170   continue                                                          
900   continue                                                          
c                                                                       
1000  mesp=mesp+nesp(kkk)                                               
1005  continue
       yoc(iso)=.false.                                                 
       yoc(jso)=.false.                                                 
       yoc(kso)=.false.                                                 
       yoc(lso)=.false.                                                 
1900  if(.not.ysyspi) go to 2000                                        
      if(yspi) go to 2000                                               
      yspi=.true.                                                       
      iso=isptr(iso)                                                    
      jso=isptr(jso)                                                    
      kso=isptr(kso)                                                    
      lso=isptr(lso)                                                    
      go to 1800                                                        
 2000 continue                                                          
c     write(6,*)ii,kk                                                   
      goto (12,23,33,44) iva                                            
c     write(6,*) ii,kk,ii                                               
      if(iva.eq.1) goto 12                                              
c     write(6,*) ii,kk,ii,kk                                            
      if(iva.eq.2) goto 23                                              
      if(iva.eq.3) goto 33                                              
      if(iva.eq.4) goto 44                                              
5000  write(6,2810)                                                     
      return                                                            
      end                                                               
