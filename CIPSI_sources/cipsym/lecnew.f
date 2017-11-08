      subroutine lecnew(ntm,npm,ntd,npd,ntt,npt,ntq,npq,ntp,npp,nth,nph,
     * nt7,np7,nto,npo)                                                 
      implicit real*8(a-h,o-x,z),logical*1(y)                           
      include 'pshf.prm'
      integer*4 ntm(1),npm(1),ntd(1),npd(1),ntt(1),                     
     1npt(1),nth(1),nph(1),                                             
     2ntq(1),npq(1),ntp(1),npp(1)                                       
      integer*4 nt7(1),np7(1),nto(1),npo(1)                             
      integer*4 diex3(4,2),diex4(4,2),masqd1(4),masqd2(4)               
      character*1 inex_min,inex_maj,irec,itasp,nex                     
      character*8 xmsg                                                  
      common kf,km,kd,kt,kq,kp,kh,k7,ko                                 
      common/cipsi/escf,test,teseff,tdet,twdet,ndiag,norb,mnorb,
     1             noca,nocb,noc1,noc2,igela,igelb,ncf,nsym,isym,isz,
     2             ntrsy,metat,lecdet,maxdet,yprt,yion,ypertu,ybrd,yen,
     3             ywvi,ywvf,ytoul,ystkic,ywdeg
      dimension itabs(40),itab(40)                                      
      dimension itri1(20,3),itri2(20,3)                                 
      dimension itri3(15,3),itri4(15,3),itri5(6,3),itri6(6,3)           
      dimension masq1(20),masq2(20),masq3(15),masq4(15)                 
      dimension mas1(20),mas2(20),mas3(15),mas4(15)                     
      dimension inex_min(9),inex_maj(9),irec(14),itasp(80)             
      dimension xmsg(8)                                                 
      data masq1,masq2/4*1,3*0,6*1,3*0,4*1,1,1,0,1,1,0,1,1,0,1,1,0,     
     &1,1,0,1,1,0,1,1/                                                  
      data itri1/1,3*2,6*1,6*2,3*1,2,4*1,3*2,3*1,3*2,3*1,4*2,7*1,3*2,   
     &3*1,7*2/                                                          
      data itri2/1,2,1,1,2,1,1,2,1,1,2,2,1,2,2,1,2,2,1,2,               
     &1,1,2,1,1,2,1,1,2,1,2,1,2,2,1,2,2,1,2,2,                          
     *1,1,1,2,1,1,2,1,1,2,1,2,2,1,2,2,1,2,2,2/                          
      data masq3,masq4/6*1,3*0,4*1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,3*1/     
      data itri3,itri4/3*1,3*2,6*1,2*2,7*1,3*2,3*1,2,1,2,9*1,3*2,       
     &1,2,2,2,1,1,2,2,1,2,2,1,2,2,1,2,   2,2,1,2,1,2,1,2,2,1,2,2,1,2,   
     *2,2,2,1,1,2,1,2,2,1,2,2,1,2,2,2,2,2/                              
      data itri5,itri6/3*1,2,6*1,2,6*1,2,1,6*2,1,6*2,1,3*2/             
      data mas1,mas2/1,3*0,12*1,3*0,1,1,0,4*1,0,2*1,2*0,2*1,0,4*1,0,1/  
      data mas3,mas4/3*1,2*0,4*1,0,2*1,0,2*1,0, 2*1,0,2*1,0,4*1,0,3*1/  
      data diex3,diex4/3*1,2,1,1,2,1,1,3*2,2,1,2,2/                     
      data masqd1,masqd2/3*1,0,0,3*1/                                   
      data inex_min/'f','m','d','t','q','p','h','7','o' /               
      data inex_maj/'F','M','D','T','Q','P','H','7','O' /               
      data irec/'0','1','2','3','4','5','6','7','8','9',' ',',','+','-'/
      data xmsg/'    mono','      di','     tri','  quadri','   penta', 
     *'    hexa','   hepta','    octo'/                                 
 5000 format(/5x,'le nombre de determinants',a8,'excites est limite a ',
     * i5,' programme stoppe')                                          
      maxm=ndetz                                                        
      maxd=ndetz                                                        
      maxt=ndetz                                                        
      maxq=ndetz                                                        
      maxp=ndetz                                                        
      maxh=ndetz                                                        
      max7=ndetz                                                        
      max8=ndetz                                                        
      ncart=0                                                           
      kf=0                                                              
      km=0                                                              
      kd=0                                                              
      kt=0                                                              
      kq=0                                                              
      kp=0                                                              
      kh=0                                                              
      ki=0                                                              
      k7=0                                                              
      ko=0                                                              
700   continue                                                          
      k=kf+km+kd/2+kt/3+kq/4+kp/5+kh/6                                  
      k=k +k7/7+k8/8                                                    
      if(k.eq.0) go to 5557                                             
      l=k-ki                                                            
c                                                                       
c                                                                       
5557  continue                                                          
      read(5,705,end=11,err=10000) itasp                                
705   format(80a1)                                                      
      ncart=ncart+1                                                     
c                                                                       
c                                                                       
5558   ki=k                                                             
      nex=itasp(1)                                                      
      do 20 i=1,9                                                       
      if(nex.eq.inex_min(i).or.nex.eq.inex_maj(i)) goto 21              
20    continue                                                          
      go to 11                                                          
21    nec=i-1                                                           
      yb=.false.                                                        
      l=0                                                               
      nur=0                                                             
      il=10                                                             
      ll=0                                                              
      do 40 i=2,80                                                      
      do 25 j=1,14                                                      
      if(itasp(i).eq.irec(j)) go to 26                                  
25    continue                                                          
      write(6,60)                                                       
60    format(' caractere non permis')                                   
      stop                                                              
26    if(j.gt.10) go to 30                                              
      nur=nur*il+j-1                                                    
c                                                                       
      yb=.false.                                                        
      go to 40                                                          
30    if(nur.eq.0) go to 32                                             
      l=l+1                                                             
      itab(l)=nur                                                       
      nur=0                                                             
32    jj=j-11                                                           
c                                                                       
      if(jj.eq.0) go to 34                                              
      jj=jj-1                                                           
      if(jj.ne.0) go to 35                                              
34    if(yb) go to 710                                                  
      yb=.true.                                                         
      go to 40                                                          
35    ll=ll+1                                                           
      itabs(ll)=(jj-1)*norb                                             
      yb=.false.                                                        
40    continue                                                          
710   if(ll.eq.0.or.ll.eq.l) go to 719                                  
      write(6,711)                                                      
711   format(' pas assez d informations en spin')                       
      go to 10000                                                       
719   ncos=ll                                                           
      i=nec+1                                                           
      maxe=l                                                            
c                                                                       
c                                                                       
      go to (720,730,750,900,2000,3000,4000,7000,8000),i                
c                                                                       
c                                                                       
c fondamental                                                           
720   if(isz.ne.0) go to 700                                            
      kf=kf+1                                                           
      if(kf.eq.1) go to 700                                             
      kf=kf-1                                                           
725   write(6,726)                                                      
726   format(' fondamental deja present')                               
      go to 700                                                         
c                                                                       
c                                                                       
c     mono                                                              
730   if(isz.gt.1) go to 700                                            
      ii=maxe                                                           
737   continue                                                          
      it1=itab(1)                                                       
      do 740 j=2,ii                                                     
c si isz=1 1spin alpha--> 1spin  beta                                   
      km=km+1                                                           
      ip1=itab(j)                                                       
      ntm(km)=it1                                                       
      npm(km)=ip1+isz*norb                                              
740   continue                                                          
      if(km.le.maxm) go to 700                                          
      write(6,5000) xmsg(1),maxm                                        
      stop                                                              
c                                                                       
c                                                                       
c                                                                       
c diexcitees                                                            
750   continue                                                          
      if(isz.gt.2) go to 700                                            
      ii=maxe                                                           
      if(ii.le.2) go to 700                                             
      np1=itab(1)                                                       
      np2=itab(2)                                                       
      do 760 j=3,ii,2                                                   
      nt1=itab(j)                                                       
      nt2=itab(j+1)                                                     
      if(ncos.ne.0) go to801                                            
      if(isz-1) 7559,770,795                                            
7559  is1=1                                                             
      is2=2                                                             
      is3=1                                                             
      is4=2                                                             
756   kd=kd+2                                                           
      kd1=kd-1                                                          
      ntd(kd1)=np1+(is1-1)*norb                                         
      ntd(kd )=np2+(is2-1)*norb                                         
      npd(kd1)=nt1+(is3-1)*norb                                         
      npd(kd )=nt2+(is4-1)*norb                                         
757   continue                                                          
      if(nt1.eq.nt2.or.np1.eq.np2) go to 760                            
      if(is4.eq.1) go to 758                                            
      is2=1                                                             
      is4=1                                                             
      go to 756                                                         
758   if(is3.eq.2) go to 760                                            
      is2=2                                                             
      is3=2                                                             
      go to 756                                                         
770   continue                                                          
      do 790 l=1,4                                                      
      if(np1.ne.np2) go to 771                                          
      if(masqd1(l).ne.0) go to 790                                      
771   if(nt1.ne.nt2) go to 772                                          
      if(masqd2(l).ne.0) go to 790                                      
772   kd=kd+1                                                           
      ntd(kd)=np1+   (diex3(l,1)-1)*norb                                
      npd(kd)=nt1+   (diex4(l,1)-1)*norb                                
      kd=kd+1                                                           
      ntd(kd)=np2+   (diex3(l,2)-1)*norb                                
      npd(kd)=nt2+   (diex4(l,2)-1)*norb                                
790   continue                                                          
       go to 760                                                        
795   if(np1.eq.np2.or.nt1.eq.nt2) go to 760                            
      kd=kd+1                                                           
      ntd(kd)=np1                                                       
      npd(kd)=nt1+norb                                                  
      kd=kd+1                                                           
      ntd(kd)=np2                                                       
      npd(kd)=nt2+norb                                                  
      go to 760                                                         
801   kd=kd+1                                                           
      ntd(kd)=np1+itabs(1)                                              
      npd(kd)=nt1+itabs(j)                                              
      kd=kd+1                                                           
      ntd(kd)=np2+itabs(2)                                              
      npd(kd)=nt2+itabs(j+1)                                            
760   continue                                                          
       if(kd/2.le.maxd) go to 700                                       
      write(6,5000) xmsg(2),maxd                                        
      stop                                                              
c                                                                       
c                                                                       
c triexcitees                                                           
900   continue                                                          
      if(isz.gt.3) go to 700                                            
      ii=maxe                                                           
      if(ii.le.3) go to 700                                             
      np1=itab(1)                                                       
      np2=itab(2)                                                       
      np3=itab(3)                                                       
903   ipo=0                                                             
      if(np1.eq.np2.or.np2.eq.np3) ipo=1                                
      if(np1.eq.np3) ipo=-1                                             
      do 1000 j=4,ii,3                                                  
      jp1=j+1                                                           
      jp2=j+2                                                           
      nt1=itab(j)                                                       
      nt2=itab(jp1)                                                     
      nt3=itab(jp2)                                                     
913   ito=0                                                             
      if(nt1.eq.nt2.or.nt2.eq.nt3) ito=1                                
      if(nt1.eq.nt3) ito=-1                                             
      if(ncos.eq.0) go to 928                                           
      do 910 l=1,3                                                      
      kt=kt+1                                                           
      ntt(kt)=itab(l)+itabs(l)                                          
910   npt(kt)=itab(j+l-1)+itabs(j+l-1)                                  
      go to 1000                                                        
928   iszp=isz+1                                                        
      go to (929,950,990,980),iszp                                      
929   continue                                                          
      do 945 l=1,20                                                     
      if(ipo) 930,933,932                                               
930   if(mas1(l).ne.0) go to 945                                        
      go to 933                                                         
932   if(masq1(l).ne.0) go to 945                                       
933   if(ito) 934,936,935                                               
934   if(mas2(l).ne.0) go to 945                                        
      go to 936                                                         
935   if(masq2(l).ne.0) go to 945                                       
936   continue                                                          
      do 938 ll=1,3                                                     
      kt=kt+1                                                           
      ntt(kt)=itab(ll)    +(itri1(l,ll)-1)*norb                         
938   npt(kt)=itab(j+ll-1)+(itri2(l,ll)-1)*norb                         
945   continue                                                          
      go to 1000                                                        
950   continue                                                          
c                      sz=1   spin+    -->   spin-                      
      do 975 l=1,15                                                     
      if(ipo) 960,962,961                                               
960   if(mas3(l).ne.0) go to 975                                        
      go to 962                                                         
961   if(masq3(l).ne.0) go to 975                                       
962   if(ito) 963,966,964                                               
963   if(mas4(l).ne.0) go to 975                                        
      go to 966                                                         
964   if(masq4(l).ne.0) go to 975                                       
966   continue                                                          
      do 968 ll=1,3                                                     
      kt=kt+1                                                           
      ntt(kt)=itab(ll)    +(itri3(l,ll)-1)*norb                         
968   npt(kt)=itab(j+ll-1)+(itri4(l,ll)-1)*norb                         
975   continue                                                          
      go to 1000                                                        
c                                                                       
c                    sz=3     +++  -->  ---                             
980   continue                                                          
      if(ipo.ne.0.or.ito.ne.0) go to 700                                
      do 982 ll=1,3                                                     
      kt=kt+1                                                           
      ntt(kt)=itab(ll)                                                  
982   npt(kt)=itab(j+ll-1)+norb                                         
      go to 1000                                                        
c                                                                       
c                         sz=2 +,+ --> -,-                              
990   continue                                                          
      if(ipo.ne.0.and.ito.ne.0) go to700                                
      if(ipo) 991,993,992                                               
991   l=4                                                               
      go to 997                                                         
992   l=5                                                               
      go to 997                                                         
993   if(ito) 994,996,995                                               
994   l=1                                                               
      go to 997                                                         
995   l=2                                                               
997   do 989 ll=1,3                                                     
      kt=kt+1                                                           
      ntt(kt)=itab(ll)+(itri5(l,ll)-1)*norb                             
989    npt(kt)=itab(j+ll-1)+(itri6(l,ll)-1)*norb                        
      go to 1000                                                        
996   do 999 l=1,6                                                      
      do 998 ll=1,3                                                     
      kt=kt+1                                                           
      ntt(kt)=itab(ll)+(itri5(l,ll)-1)*norb                             
998    npt(kt)=itab(j+ll-1)+(itri6(l,ll)-1)*norb                        
999   continue                                                          
1000  continue                                                          
       if(kt/3.le.maxt) go to 700                                       
      write(6,5000) xmsg(3),maxt                                        
      stop                                                              
c                                                                       
c                                                                       
c                                                                       
c     quadri (pas de generation auto)                                   
2000  ii=maxe/4                                                         
      ii=maxe                                                           
      if(ii.le.7) go to 700                                             
      do 2100 jj=5,ii,4                                                 
      j=jj-1                                                            
      do 2050 i=1,4                                                     
      kq=kq+1                                                           
      ntq(kq)=  itab(i)+itabs(i)                                        
2050  npq(kq)=  itab(i+j)+itabs(i+j)                                    
2100  continue                                                          
       if(kq/4.le.maxq) go to 700                                       
      write(6,5000) xmsg(4),maxq                                        
      stop                                                              
c                                                                       
c     penta (pas de generation auto)                                    
3000   ii=maxe/5                                                        
      ii=maxe                                                           
      if(ii.le.9) go to 700                                             
      do 3100 jj=6,ii,5                                                 
      j=jj-1                                                            
      do 3050 i=1,5                                                     
      kp=kp+1                                                           
      ntp(kp)=itab(i)+itabs(i)                                          
3050  npp(kp)=itab(i+j)+itabs(i+j)                                      
c                                                                       
3100  continue                                                          
       if(kp/5.le.maxp) go to 700                                       
      write(6,5000) xmsg(5),maxp                                        
      stop                                                              
c                                                                       
c     hexa (pas de generation auto)                                     
4000  ii=maxe/6                                                         
      ii=maxe                                                           
      if(ii.le.11) go to 700                                            
      do 4100 jj=7,ii,6                                                 
      j=jj-1                                                            
      do 4050 i=1,6                                                     
      kh=kh+1                                                           
      nth(kh)=itab(i)+itabs(i)                                          
4050  nph(kh)=itab(i+j)+itabs(i+j)                                      
4100   continue                                                         
       if(kh/6.le.maxh) go to 700                                       
      write(6,5000) xmsg(6),maxh                                        
      stop                                                              
c                                                                       
c     hepta (pas de generation auto)                                    
c                                                                       
 7000 ii=maxe/7                                                         
      ii=maxe                                                           
      if(ii.le.13)go to 700                                             
      do 7100 jj=8,ii,7                                                 
      j=jj-1                                                            
      do 7050 i=1,7                                                     
      k7=k7+1                                                           
      nt7(k7)=itab(i)+itabs(i)                                          
 7050 np7(k7)=itab(i+j)+itabs(i+j)                                      
 7100 continue
      if(k7/7.le.max7) go to 700                                        
      write(6,5000) xmsg(7),max7                                        
      stop                                                              
c                                                                       
c     octo (pas de generation auto)                                     
c                                                                       
 8000 ii=maxe                                                           
      if(ii.le.15) go to 700                                            
      do 8100 jj=9,ii,8                                                 
      j=jj-1                                                            
      do 8050 i=1,8                                                     
      ko=ko+1                                                           
      nto(ko)=itab(i)+itabs(i)                                          
 8050 npo(ko)=itab(i+j)+itabs(i+j)                                      
 8100 continue                                                          
      if(ko/8.le.maxo) go to 700                                        
      write(6,5000) xmsg(8),maxo                                        
      stop                                                              
11    continue                                                          
c                                                                       
      return                                                            
10000   write(6,10001)                                                  
10001  format(' erreur en lecture determinant')                         
      stop                                                              
      end                                                               
