      subroutine igene(k,yf,ym,yd)                                      
      implicit real*8(a-h,o-x,z),logical*1(y)                           
      include 'pshf.prm'
c                                                                       
c     generation des configurations i obtenues a partir des configurat  
c     ions k de l'i.c                                                   
c                                                                       
      common/spo/fdiag(2*doa),hdiag(ndetz),tocom(2*doa),coeff,      
     1ndef,ispin(2*doa),iwspin(2*doa),numero(metz),          
     2ne(ndetz+1),nd(ndetz+1),trou(ntpz),part(ntpz),iorb(2*doa),      
     3itsym(2*doa),its(nsymz,nsymz),yoc(2*doa),yocs(2*doa),       
     4ygen(ndetz),yocd(2*doa,ndetz),yspin,yef(ndetz,nsymz)                   
      integer*4 ne,nd,trou,part,iorb,iwspin,itsym,its,isytr,nesp,nespo  
      dimension noc(2),igels(2),norbs(2)                                
      common/cipsi/escf,test,teseff,tdet,twdet,ndiag,norb,mnorb,
     1             noca,nocb,noc1,noc2,igela,igelb,ncf,nsym,isym,isz,
     2             ntrsy,metat,lecdet,maxdet,yprt,yion,ypertu,ybrd,yen,
     3             ywvi,ywvf,ytoul,ystkic,ywdeg
      equivalence (noca,noc(1)),(igela,igels(1)),(norb,norbs(1))        
      integer*4 ntps(ntpz,2)                                            
      equivalence (ntps(1,2),part(1)),(ntps(1,1),trou(1))               
      nf=ncf+1                                                          
      ndk=nd(k)                                                         
      ndf=nd(nf)                                                        
      k0=k                                                              
      nec =ne(k)                                                        
      if (nec) 501,503,501                                              
  501 continue                                                          
      do 502 ir=1,nec                                                   
      trou(ndf+ir)=trou(ndk+ir)                                         
      part(ndf+ir)=part(ndk+ir)                                         
  502 continue                                                          
  503 continue                                                          
      if(nec.eq.0.and.ym)go to 999                                      
      if(nec.eq.1.and.yd) go to 999                                     
c     monoexcitations                                                   
      nec1=nec+1                                                        
      ne(nf) = nec+1                                                    
      do 522 ks=1,2                                                     
      igel=igels(ks)                                                    
      n=norbs(ks)                                                       
      nbor=noc(ks)                                                      
      nbor1=nbor+1                                                      
      do 524 ii=igel,nbor                                               
      trou(ndf+nec1)=ii                                                 
      iis=itsym(ii)                                                     
      do 526 jj=nbor1,n                                                 
c                                                                       
      if(yocd(ii,k).or.yocd(jj,k)) go to 526                            
c                                                                       
      part(ndf+nec1)=jj                                                 
      jjs=itsym(jj)                                                     
      ijs=its(iis,jjs)                                                  
      if(ijs.eq.1)call pert1(k,nf)                                      
      if(nec.eq.0) go to 526                                            
  530 do 532 kr=1,2                                                     
c     changement detrou et de particule                                 
      do 534 lk=1,nec                                                   
      la=ntps(ndk+lk,kr)                                                
      ls=ispin(la)+1                                                    
      las=itsym(la)                                                     
      if (kr-1) 536,536,538                                             
536   na1=igels(ls)                                                     
      if(ls.eq.ks) na1=ii+1                                             
      na2= noc(ls)                                                      
      go to 540                                                         
  538 na1 =noc(ls)+1                                                    
      if(ls.eq.ks) na1=jj+1                                             
      na2=norbs(ls)                                                     
540   if(na1.gt.na2) goto 534                                           
      do 542 ll=na1,na2                                                 
      if(ll.eq.la.or.ll.eq.ii.or.ll.eq.jj) go to 542                    
544   ntps(ndf+lk,kr)=ll                                                
      lls=itsym(ll)                                                     
      kls=its(las,lls)                                                  
      ijkls=its(ijs,kls)                                                
      if(ijkls.eq.1) call pert1(k,nf)                                   
      ntps(ndf+lk,kr)=la                                                
  542 continue                                                          
  534 continue                                                          
  532 continue                                                          
  526 continue                                                          
 524  continue                                                          
  522 continue                                                          
      if(nec.lt.1) go to 999                                            
c modif jp .daudey dec 89 .and.not.yion                                 
c dans le cas des ions,les monoexcites ne sont pas generes dans pertu   
c donc ils doivent etre inclus ici                                      
      if(ym.and.nec.eq.1.and..not.yion) go to 999                       
      if(nec.eq.2.and.yd) go to 999                                     
c      meme degre d'excitation                                          
c      changement d'un trou ou d'une particule                          
  548 ne(nf)=nec                                                        
      do 550 lk=1,nec                                                   
      lk1=lk+1                                                          
      do 552 kr=1,2                                                     
      la=ntps(ndk+lk,kr)                                                
      ls=ispin(la)+1                                                    
      las=itsym(la)                                                     
      if(kr.gt.1) go to 556                                             
554   na1=igels(ls)                                                     
      na2=noc(ls)                                                       
      go to 558                                                         
  556 na1=noc(ls)+1                                                     
      na2=norbs(ls)                                                     
558   do 569 ll=na1,na2                                                 
c                                                                       
      ll1=ll+1                                                          
      if(ll.eq.la) go to 569                                            
562   ntps(ndf+lk,kr)=ll                                                
      lls=itsym(ll)                                                     
      lalls=its(las,lls)                                                
      if(lalls.eq.1) call pert1(k,nf)                                   
c       changement d'un deuxieme trou ou particule                      
564   if(lk.eq.nec) go to 569                                           
  566 do 568 lkp=lk1,nec                                                
      lb=ntps(ndk+lkp,kr)                                               
      lbs=ispin(lb)+1                                                   
      if(lbs.eq.ls.and.ll.ge.na2) go to 568                             
      lbsy=itsym(lb)                                                    
      if(kr.eq.2)go to 572                                              
c      je change un 2ieme trou                                          
      nb2=noc(lbs)                                                      
      nb1=igels(lbs)                                                    
      go to 578                                                         
572   nb2=norbs(lbs)                                                    
      nb1=noc(lbs)+1                                                    
578   if(ll.eq.lb) go to 587                                            
      if(ls.eq.lbs) nb1=ll1                                             
      do 585 ii=nb1,nb2                                                 
      if(ii.eq.lb) go to 585                                            
      if(ii.eq.la) go to 585                                            
580   ntps(ndf+lkp,kr)=ii                                               
      ntps(ndf+lk ,kr)=ll                                               
      iis=itsym(ii)                                                     
      ibs=its(lbsy,iis)                                                 
      if(its(ibs,lalls).eq.1)  call pert1(k,nf)                         
  585 continue                                                          
587   continue                                                          
      ntps(ndf+lkp,kr)=lb                                               
  568 continue                                                          
569   continue                                                          
      ntps(ndf+lk,kr)=la                                                
  552 continue                                                          
  550 continue                                                          
c      changement d,une paire trou -particule                           
      do 586 lk=1,nec                                                   
      la=trou(ndk+lk)                                                   
      ls=ispin(la)+1                                                    
      las=itsym(la)                                                     
      na1=igels(ls)                                                     
      na2=noc(ls)                                                       
      do 588 lkp=1,nec                                                  
      lb=part(ndk+lkp)                                                  
      lbs=ispin(lb)+1                                                   
      nb1=noc(lbs)+1                                                    
      nb2=norbs(lbs)                                                    
      lbsy=itsym(lb)                                                    
      lalbs=its(las,lbsy)                                               
c     si le trou t la particule annihiles sont de meme spin , je peux   
c     creer un trou et une particule de spin  identique mais quelconque 
c      dans le cas contraire , le nouveau trou doit etre du meme spin   
c      que l,ancien ,la nouvelle particule du meme spin que l'ancienne  
      ntsp=ls-lbs                                                       
      do 590 ll=na1,na2                                                 
      if (ll.eq.la) goto 590                                            
      lls=itsym(ll)                                                     
592   do 594 ii=nb1,nb2                                                 
      if (ii.eq.lb) goto  594                                           
596   part(ndf+lkp)=ii                                                  
      trou(ndf+lk )=ll                                                  
      iis=itsym(ii)                                                     
      lis=its(lls,iis)                                                  
      liabs=its(lis,lalbs)                                              
      if(liabs.eq.1) call pert1(k,nf)                                   
  594 continue                                                          
  590 continue                                                          
      if (ntsp) 588,598,588                                             
  598 lsp=-ls+3                                                         
c        changement de spin                                             
      nocu =noc(lsp)                                                    
      nocu1=nocu+1                                                      
      igelu=igels(lsp)                                                  
      norbu=norbs(lsp)                                                  
      do 600 ll=igelu,nocu                                              
      lls=itsym(ll)                                                     
      do 602 ii=nocu1,norbu                                             
      trou(ndf+lk)=ll                                                   
      part(ndf+lkp)=ii                                                  
      iis=itsym(ii)                                                     
      lis=its(lls,iis)                                                  
      liabs=its(lis,lalbs)                                              
      if(liabs.eq.1)  call pert1(k,nf)                                  
 602  continue                                                          
  600 continue                                                          
588   part(ndf+lkp)=lb                                                  
      trou(ndf+lk )=la                                                  
  586 continue                                                          
      if(nec.eq.3.and.yd)go to 999                                      
      if(nec.eq.2.and.ym) go to 999                                     
c      etats a une excitation de moins                                  
      ne(nf)=nec-1                                                      
      do 604 lk=1,nec                                                   
      la=trou(ndk+lk)                                                   
      ls=ispin(la)+1                                                    
      las=itsym(la)                                                     
      do 606 lkp=1,nec                                                  
      lb=part(ndk+lkp)                                                  
      lbs=ispin(lb)+1                                                   
      if (ls.ne.lbs) goto 606                                           
      lbsy=itsym(lb)                                                    
c      simple desexcitation                                             
      trou(ndf+lk)=0                                                    
      part(ndf+lkp)=0                                                   
      lalbs=its(las,lbsy)                                               
      if(lalbs.eq.1) call cotra(k,nf)                                   
c      changeons en plus un trou ou une particule                       
c      sale prof , va                                                   
      do 616 lkt=1,nec                                                  
c       si les deux trous ont le meme spin (ls=lts) , ils jouent un role
c       identique ,donc lkt>lk                                          
c      si les 2 trous ont des spins differents , ils jouent des roles   
c      differents et les boucles sont sans restrictions                 
      lat=trou(ndk+lkt)                                                 
      lts=ispin(lat)+1                                                  
      lbt=part(ndk+lkt)                                                 
c      chgt de trou                                                     
      if(lkt-lk) 619,624,619                                            
  619 kr=1                                                              
      ici=itsym(lat)                                                    
      na1=igels(lts)                                                    
      na2=noc(lts)                                                      
      nevi2=lat                                                         
      if (ls-lts) 618,620,618                                           
  618 nevi1=0                                                           
      do 1010 ii=na1,na2                                                
      if(ii.eq.nevi1.or.ii.eq.nevi2) go to 1010                         
      ntps(ndf+lkt,kr)=ii                                               
      trou(ndf+lk)=0                                                    
      part(ndf+lkp)=0                                                   
      iici=its(ici,itsym(ii))                                           
      if(its(iici,lalbs).eq.1) call cotra(k,nf)                         
1010  continue                                                          
      ntps(ndf+lkt,kr)=ntps(ndk+lkt,kr)                                 
      go to 624                                                         
  620 if (lkt-lk) 624,624,628                                           
  628 nevi1=la                                                          
      do 1020 ii=na1,na2                                                
      if(ii.eq.nevi1.or.ii.eq.nevi2) go to 1020                         
      ntps(ndf+lkt,kr)=ii                                               
      trou(ndf+lk)=0                                                    
      part(ndf+lkp)=0                                                   
      iici=its(ici,itsym(ii))                                           
      if(its(iici,lalbs).eq.1) call cotra(k,nf)                         
1020  continue                                                          
      ntps(ndf+lkt,kr)=ntps(ndk+lkt,kr)                                 
c      changement de particule                                          
624   lps=ispin(lbt)+1                                                  
      kr=2                                                              
      ici=itsym(lbt)                                                    
      na1=noc(lps)+1                                                    
      na2=norbs(lps)                                                    
      nevi2=lbt                                                         
      if (lkt-lkp) 621,616,621                                          
  621 if(ls-lps) 623,625,623                                            
  623 nevi1=0                                                           
      do 1030 ii=na1,na2                                                
      if(ii.eq.nevi1.or.ii.eq.nevi2) go to 1030                         
      ntps(ndf+lkt,kr)=ii                                               
      trou(ndf+lk)=0                                                    
      part(ndf+lkp)=0                                                   
      iici=its(ici,itsym(ii))                                           
      if(its(iici,lalbs).eq.1) call cotra(k,nf)                         
1030  continue                                                          
      ntps(ndf+lkt,kr)=ntps(ndk+lkt,kr)                                 
      go to 616                                                         
  625 if  (lkt-lkp) 616,616,627                                         
  627 nevi1=lb                                                          
      do 1040 ii=na1,na2                                                
      if(ii.eq.nevi1.or.ii.eq.nevi2) go to 1040                         
      ntps(ndf+lkt,kr)=ii                                               
      trou(ndf+lk)=0                                                    
      part(ndf+lkp)=0                                                   
      iici=its(ici,itsym(ii))                                           
      if(its(iici,lalbs).eq.1) call cotra(k,nf)                         
1040  continue                                                          
      ntps(ndf+lkt,kr)=ntps(ndk+lkt,kr)                                 
  616 continue                                                          
606   part(ndf+lkp)=lb                                                  
604   trou(ndf+lk )=la                                                  
c      deux desexcitations                                              
c modif jp daudey juillet 89 yd.or.yf                                   
c erreur ne se produisant que dans le cas ou l'on traite un etat de meme
c symetrie que le fondamental sans que celui-ci soit dans s.            
c il est alors engendre deux fois                                       
      if(nec.eq.2.and.(yd.or.yf))go to 999                              
      if(nec.eq.3.and.ym)go to 999                                      
      if(nec.eq.4.and.yd)go to 999                                      
      necm=nec-1                                                        
      ne(nf)=nec-2                                                      
c      il faut annihiler 2 trous et 2 particules , pour ca , prendre    
c      toutes les paires de trous (i.e. lt1) lt2), et toutes les paires 
c      de particules (lp1)lp2). la compatibilite de spin exige que      
c      si les spins de lt1 et lt2 sont les memes , les spins de lp1 et  
c      lp2 soient aussi egaux a celui de lt1 et lt2                     
c      si les spisn de lt1 et lt2 sont differnts , il faut que las      
c       spins  de lp1et lp2 le soient aussi                             
  634 do 636 lt1 =1,necm                                                
      lt1p=lt1+1                                                        
      lat=trou(ndk+lt1)                                                 
c                                                                       
      lats=ispin(lat)+1                                                 
c                                                                       
      ilat=itsym(lat)                                                   
      do 638 lt2=lt1p,nec                                               
      lbt=trou(ndk+lt2)                                                 
c                                                                       
      lbts=ispin(lbt)+1                                                 
c                                                                       
      nts=lats-lbts                                                     
      ilbt=itsym(lbt)                                                   
      ilab=its(ilat,ilbt)                                               
      do 640 lp1=1,necm                                                 
      lp1p=lp1+1                                                        
c                                                                       
      lcp=part(ndk+lp1)                                                 
c                                                                       
      lcps=ispin(lcp)+1                                                 
      ilcp=itsym(lcp)                                                   
      do 642 lp2=lp1p,nec                                               
c                                                                       
      ldp=part(ndk+lp2)                                                 
c                                                                       
      ldps=ispin(ldp)+1                                                 
      if(nts.ne.0) go to 644                                            
646   if(lcps.ne.lats.or.ldps.ne.lats) go to 642                        
      go to 650                                                         
644   if(lcps.eq.ldps) go to 642                                        
650   trou(ndf+lt1)=0                                                   
      trou(ndf+lt2)=0                                                   
      part(ndf+lp1)=0                                                   
      part(ndf+lp2)=0                                                   
      ildp=itsym(ldp)                                                   
      ilcd=its(ilcp,ildp)                                               
      ilabcd=its(ilab,ilcd)                                             
      if(ilabcd.eq.1) call cotra(k,nf)                                  
642   part(ndf+lp2)=ldp                                                 
640   part(ndf+lp1)=lcp                                                 
638   trou(ndf+lt2)=lbt                                                 
636   trou(ndf+lt1)=lat                                                 
  999 continue                                                          
      return                                                            
      end                                                               
