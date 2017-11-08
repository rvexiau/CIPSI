      subroutine pertu                                                  
      implicit real *8(a-h,o-x,z),logical*1(y)                          
      include 'pshf.prm'
      common/spo/fdiag(2*doa),hdiag(ndetz),tocom(2*doa),coeff,      
     1ndef,ispin(2*doa),iwspin(2*doa),numero(metz),          
     2ne(ndetz+1),nd(ndetz+1),trou(ntpz),part(ntpz),iorb(2*doa),      
     3itsym(2*doa),its(nsymz,nsymz),yoc(2*doa),yocs(2*doa),       
     4ygen(ndetz),yocd(2*doa,ndetz),yspin                                    
      common/degene/nesp(ndetz),nespo(ndetz),isytr(ndetz,nsymz),kkbrd,  
     *ysaut(ndetz),imul(nmulz),yef(ndetz,nsymz)                         
      integer*4 ne,nd,trou,part,iorb,iwspin,itsym,its,isytr,nesp,nespo  
      integer*4 imul                                                    
      common/cipsi/escf,test,teseff,tdet,twdet,ndiag,norb,mnorb,
     1             noca,nocb,noc1,noc2,igela,igelb,ncf,nsym,isym,isz,
     2             ntrsy,metat,lecdet,maxdet,yprt,yion,ypertu,ybrd,yen,
     3             ywvi,ywvf,ytoul,ystkic,ywdeg
      common/ic/c(metz*ndetz),e(ndetz),eb(metz),e2b(metz),emp(metz),    
     1teste(metz),e2mp(metz),dpivot,ietat,ishif,e2vp(metz),             
     2fmpb(2*doa),hefmp(nhefz),hefvp(nhefz),hefb(nhefz),              
     3psimp(metz),psivp(metz)                                           
      common/tome/taum(metz),yselm,ysele,ydp                            
      dimension noc(2),igels(2),norbs(2)                                
      equivalence (noca,noc(1)),(igela,igels(1)),(norb,norbs(1))        
      write(6,1888)                                                     
      write(6,999)                                                      
      write(6,1888)                                                     
1888  format(1x,77('*'),/)                                              
      if(yen) then                                                      
      write(6,*) ' selection envp '                                     
      else                                                              
      write(6,*) ' selection mpb'                                       
      endif                                                             
      if(ysele) then                                                    
         write(6,*) ' selection sur l energie'                          
       else                                                             
         write(6,*) ' selection sur la fonction'                        
       endif                                                            
      if(yselm) write(6,*) ' seuil renormalise par etat'                
      write(6,*) 

      write(6,1999)                                                     
      write(6,2001)                                                     
      write(6,2002)                                                     
      write(6,2003)                                                     
      write(6,2004)                                                     
      write(6,*)                                                        
      write(6,1059)                                                     
1059  format(1x,' m       cmi     emi      hmi    ei',                  
     * 10x,'determinant i')                                             
      nsy=isym                                                          
      nf=ncf+1                                                          
      yd=.false.                                                        
      ym=.false.                                                        
      yf=.false.                                                        
      do400k=1,ncf                                                      
      if(ne(k).eq.1)ym=.true.                                           
      if(ne(k).eq.0)yd=.true.                                           
      if(ne(k).eq.2)yf=.true.                                           
400   continue                                                          
      k=0                                                               
      if(ym)yf=ym                                                       
      if(yion) go to 450                                                
      if(.not.yf) go to 450                                             
      if(yd.or.isz.ne.0) go to 410                                      
c                                                                       
      if(nsy.ne.1) go to 410                                            
      ne(nf)=0                                                          
999   format(1x,'selection des determinants i ameliorant l etat m',     
     *' et de leur contributions')                                      
1999  format(1x,'emi est la contribution du determinant i a l ene',     
     *'rgie de l etat m')                                               
2001  format(1x,'cmi est le coefficient du determinant i dans',         
     *'la correction d ordre 1 a la fonction de l etat m')              
2002  format(1x,'hmi est l element de matrice <m/h/i>')                 
2003  format(1x,'ei est l energie du determinant i')                    
2004  format(1x,'l energie de reference est celle de l etat du',        
     *' vide',/)                                                        
      write(6,403)                                                      
403   format(1x,'reference  ')                                          
      call pert1(0,nf)                                                  
410   ne(nf)=1                                                          
      if(.not.ym) go to 450                                             
      write(6,404)                                                      
404   format('0monoexcitations')                                        
      ndf=nd(nf)+1                                                      
      is=0                                                              
       if(isz.ne.0) go to 426                                           
415   is=is+1                                                           
      igel=igels(is)                                                    
      itm=noc(is)                                                       
      ip1=itm+1                                                         
      n=norbs(is)                                                       
416   do 420 it=igel,itm                                                
      ist=itsym(it)                                                     
      do 420 ip=ip1,n                                                   
      isp=itsym(ip)                                                     
      if(its(ist,isp).ne.nsy) go to 420                                 
      part(ndf)=ip                                                      
      trou(ndf)=it                                                      
      call pert1(0,nf)                                                  
420   continue                                                          
      if(is.lt.2.and.isz.eq.0) go to 415                                
      go to 450                                                         
426   igel=igela                                                        
      n=mnorb                                                           
      ip1=nocb+1                                                        
      itm=noca                                                          
      go to 416                                                         
450   continue                                                          
      if((.not.ym).or.(yion)) write(6,404)                              
      do 500 k=1,ncf                                                    
      if(.not.ygen(k)) go to 500                                        
      if(ysaut(k)) go to 500                                            
      kkbrd=nespo(k)                                                    
      call igene(k,yf,ym,yd)                                            
500   continue                                                          
      return                                                            
      end                                                               
