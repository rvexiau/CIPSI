      subroutine dipole(typdip)
      implicit real*8 (a-h,o-z)                                         
      character*40 typener
      logical*1 yuni,ycv,yiet,yjet,ycou
      include 'pshf.prm' 
      character*40 typdip
      common/uni/yuni
      common/polar/calfa(dc),ycv
      common/cd/r(ncouz*(doa**2)),ndd,ndt,ncf1,ncf2,metat,icf2,ncou,    
     &ibeg(ncouz),ietat(ncouz),jetat(ncouz),ycou(metz,metz),
     &yiet(metz),yjet(metz)                                             
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,zan(dc),c(3,dc)          
      common/iofile/ir,iw,ip,is,iq,ih,iv,if2,if3,if4,ioda(19)           
      common/dipxyz/xs(doas),ys(doas),zs(doas)                          
      common /inout/lec,imp                                             
      common/ic/e(metz2)
      dimension dmuel(ncouz,3),dmunuc(ncouz,3)                          
      dimension dmxyz(ncouz),rob(doa**2)
cdudu
      dimension ff(ncouz)
cdudu
      dimension epol(doa*(doa+1)/2,3)
      data fac /2.54158059d+00/                                         
c      if(ich.ne.0) return                                              
c                                                                       
c     ----- calculate c                                                 
      call xyzdip                                                       
c      write(6,*) jetat
c                                                                       
c     ----- electronic contribution to dipole moment -----              
c                                                                       
      write(6,*)                                                        
      write(6,*) ' 1 ua =2.54158059 debyes'                             
      write(6,*)                                                        
      nusq=num*num                                                      
      do 2000 icou=1,ncou                                               
      dmx=0                                                             
      dmy=0                                                             
      dmz=0                                                             
      dmxnuc=0                                                          
      dmynuc=0                                                          
      dmznuc=0                                                          
      m1=ietat(icou)                                                    
      m2=jetat(icou)                                                    
      ind=ibeg(icou)  
c calcul de la trace de la matrice densite  
c      traro=0.d0
c      do io=1,num
c      traro=traro+r(ind+(io-1)*num+io)
c      enddo
c      write(6,*)' trace de ro atom pour ietat jetat ',m1,m2,
c     *   traro
c
      indm1=ibeg(icou)
      ilim=num*num
      do i=1,ilim
      rob(i)=r(indm1+i)
      enddo
      dmuel(icou,1)=-tracp(rob,xs,num)                                  
      dmuel(icou,2)=-tracp(rob,ys,num)                                  
      dmuel(icou,3)=-tracp(rob,zs,num)                                  
      dmunuc(icou,1)=0.d0                                               
      dmunuc(icou,2)=0.d0                                               
      dmunuc(icou,3)=0.d0                                               
c                                                                       
c      ------  nuclear contribution to dipole moment    ----------      
c                                                                       
      if(yuni.and.(ietat(icou).eq.jetat(icou)))then                    
        do 1800 j=1,nat                                                 
           dmxnuc=dmxnuc+zan(j)*c(1,j)                                  
           dmynuc=dmynuc+zan(j)*c(2,j)                                  
1800       dmznuc=dmznuc+zan(j)*c(3,j)
c----------------------------------------------------
c
c     correction nucleaire (coeur valence)
c
      if (.not.ycv) go to 950
      naa=num*(num+1)/2
      rewind 17
      do 830 ic=1,nat
      alfad=calfa(ic)
      if (dabs(alfad).lt.1.0d-06) go to 830
      read (17) (epol(i,1),i=1,naa)
      read (17) (epol(i,2),i=1,naa)
      read (17) (epol(i,3),i=1,naa)
      read (17)
  830 continue
      dcorx=0.0d0
      dcory=0.0d0
      dcorz=0.0d0
      do 840 i=1,nat
      alfad=calfa(i)
      if (alfad.lt.1.0d-6) go to 840
      px=0.0d0
      py=0.0d0
      pz=0.0d0
      do 835 j=1,nat
      if (j.eq.i) go to 835
      rd=(c(1,i)-c(1,j))**2 + (c(2,i)-c(2,j))**2 + (c(3,i)-c(3,j))**2
      rd=zan(j)/(rd*dsqrt(rd))
      px=px+rd*(c(1,j)-c(1,i))
      py=py+rd*(c(2,j)-c(2,i))
      pz=pz+rd*(c(3,j)-c(3,i))
  835 continue
      dcorx=dcorx+alfad*px
      dcory=dcory+alfad*py
      dcorz=dcorz+alfad*pz
  840 continue
      dmxnuc=dmxnuc-dcorx
      dmynuc=dmynuc-dcory
      dmznuc=dmznuc-dcorz
c      print *,icou,'dmnuc   ',dmxnuc+dcorx,dmynuc+dcory,dmznuc+dcorz
c      print *,icou,'dmnuc cv',dmxnuc,dmynuc,dmznuc
950   continue
c------------------------------------------------------
           dmunuc(icou,1)=dmxnuc                                        
           dmunuc(icou,2)=dmynuc                                        
           dmunuc(icou,3)=dmznuc                                        
      endif                                                             
2000  continue                                                          
      write(6,*)                                                        
      write(6,*) ' contribution electronique en ua'               
      write(6,9999)                                                     
      do 2020 i=1,ncou                                                  
      write(6,9998)ietat(i),jetat(i),                                   
     * dmuel(i,1),dmuel(i,2),dmuel(i,3)                                 
2020  continue                                                          
      write(6,*)                                                        
      write(6,*) ' contribution nucleaire en ua'                        
      write(6,9999)                                                     
      do 2030 i=1,ncou                                                  
      write(6,9998)ietat(i),jetat(i),                                   
     * dmunuc(i,1),dmunuc(i,2),dmunuc(i,3)                              
2030  continue                                                          
      write(6,*)                                                        
      write(6,*) ' nombre d etats traites'
      write(6,*) ncou
      write(6,*)                                                        
      write(6,*) ' contribution totale en ua'                           
      write(6,*) ' forces d oscillateurs en absorption'
      write(6,9997)                                                     
      do 2040 i=1,ncou                                                  
      dmx=dmuel(i,1)+dmunuc(i,1)                                        
      dmy=dmuel(i,2)+dmunuc(i,2)                                        
      dmz=dmuel(i,3)+dmunuc(i,3)                                        
      dmxyz(i)=dmx*dmx+dmy*dmy+dmz*dmz
      tua=e(metz+jetat(i))-e(ietat(i))
      tev=tua*27.2103d0
      tcm=tua*219474.63137d0
      ff(i)=(2.d0/3.d0)*dmxyz(i)*tua
      write(6,9998)ietat(i),jetat(i),                                   
cdudu* dmx,dmy,dmz,ff,tev,tcm
     * dmx,dmy,dmz,tcm,tev,ff(i)
2040  continue      
cdudu
      typener='energie de polarisation des coeurs'
      call rdener(typener,epolnuc,1)
      write(6,*)
      write(6,*) ' energie de polarisation des coeurs ',epolnuc
      write(6,*) ' ietat(1) energie totale du fondamental'
      write(6,9996)ietat(1), e(ietat(1))+epolnuc
      write(6,*)
      write(6,*) ' jetat(i) energies totales des etats excites'
      do 2050 i=1,ncou                                                  
      write(6,9996)jetat(i), e(metz+jetat(i))+epolnuc
2050  continue
c  ecriture sur ft99 de ncou,ecor(i),f si calcul du fondamental, 
c                                      si e(ietat(1))=e(metz+jetat(1))     
        if(e(ietat(1)).eq.e(metz+jetat(1))) then
        open(unit=99,status='UNKNOWN',form='FORMATTED',file='ft99')
        write(99,*) ncou
          do i=1,ncou
          write(99,*) e(metz+jetat(i))+epolnuc
          end do
          do i=1,ncou
          write(99,*) ff(i)
          end do
        end if
cdudu

c stockage sur la file 33
      call stkener(typdip,dmxyz,ncou)
c      ires=23
c      open(unit=ires,file='/disc/romu/enertemp',status='unknown',
c     &  form='unformatted')
c      rewind ires
c      read(ires) mmetat,(e,i=1,mmetat)
c1067  write(ires) ncou,(dmxyz(i),i=1,ncou)
      write(1) ncou,(ietat(i),jetat(i),                                 
     *dmuel(i,1),dmuel(i,2),dmuel(i,3),i=1,ncou)                        
 9999 format(/,1x,'etats',9x,'dmx',7x,'dmy',7x,'dmz')
 9997 format(/,1x,'etats',9x,'dmx',7x,'dmy',7x,'dmz',
cdudu*          7x,'fe',6x,'te(eV)',3x,'te(cm-1)')                  
     *          5x,'te(cm-1)',3x,'te(eV)',5x,'fe')                  
cdudu9998 format(2i4,1x,5f10.4,1x,f10.1)
cCc  9998 format(2i4,1x,3f10.4,1x,f10.1,2f10.4)
 9998 format(2i4,1x,3f16.10,1x,f10.1,2f10.4)
 9996 format(i4,1x,f18.12,1x,f12.6)
      return                                                            
      end                                                               
