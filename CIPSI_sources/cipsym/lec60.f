      subroutine lec60                                                  
     *(ntm,npm,ntd,npd,ntt,npt,ntq,npq,ntp,npp,nth,nph,nt7,np7,nto,npo) 
c                                                                       
c  lecture des determinants sur la file 60                              
c                                                                       
c    tdet :seuil de selection                                           
c    maxdet : nombre maximal selectionne pour s                         
c                                                                       
      implicit real*8(a-h,o-p,r-x,z),logical*1(y)                       
      include 'pshf.prm'
      integer*4 ntm(ndetz),npm(ndetz),ntd(2*ndetz),npd(2*ndetz),        
     *ntt(3*ndetz),npt(3*ndetz),ntq(4*ndetz),npq(4*ndetz)               
     *,ntp(5*ndetz),npp(5*ndetz),nth(6*ndetz),nph(6*ndetz)              
     *,nt7(7*ndetz),np7(7*ndetz),nto(8*ndetz),npo(8*ndetz)              
      common kf,km,kd,kt,kq,kp,kh,k7,ko                                 
      common/cipsi/escf,test,teseff,tdet,twdet,ndiag,norb,mnorb,
     1             noca,nocb,noc1,noc2,igela,igelb,ncf,nsym,isym,isz,
     2             ntrsy,metat,lecdet,maxdet,yprt,yion,ypertu,ybrd,yen,
     3             ywvi,ywvf,ytoul,ystkic,ywdeg
      common/det/my,ny,idet(16*lungo),nexst(lungo),vst(lungo)
     * ,evp(lungo*metz),ntot
      common/infogv/nbdet(9),iord60(ndetz),ynoic                        
      dimension nexs(ndetz)                             
      dimension itab(40)                         
c                                                                       
      ncart=0                                                           
      kf=0                                                              
      km=0                                                              
      kd=0                                                              
      kt=0                                                              
      kq=0                                                              
      kp=0                                                              
      kh=0                                                              
      k7=0                                                              
      ko=0                                                              
      ki=0                                                              
      ndet=0                                                            
      rewind 60
c-- ou est la premiere mono, la premiere di...?                         
      do 344 ij=1,9                                                     
        nbdet(ij)=0                                                     
344   continue                                                          
600   read(60,end=3487) mw,nw,(idet(i),i=1,mw),(nexst(i),i=1,nw),       
     &(vst(i),i=1,nw)                                                       
c                                                                       
c                                                                       
c                                                                       
c                                                                       
c             
      myy=0
      do  lstf=1,nw
        nec=nexst(lstf)
        if(nec.gt.8) then
          write(6,*) 'erreur degre d excitation',nec,' depasse les ',
     *    'possibilites du programme'
          stop
        end if
        if(vst(lstf).lt.tdet) then
	  myy=myy+2*nec
	  else
          ndet=ndet+1
          if(ndet.gt.maxdet) then      
            if(lecdet.eq.3)return   
            if((lecdet.eq.2).or.(lecdet.eq.5)) then  
              write(6,799)                         
799   format(/,1x,' trop de determinants selectionnes sur la file60',   
     */,' augmenter le seuil de selection  tdet ou le nombre maxdet',   
     */,' maxdet maximum=2000')                                         
              stop                                
            endif                                
          endif                                 
          if(nec.ne.0)then
	    do  kstf=1,nec
	      myy=myy+1
              itab(kstf)=idet(myy)
	      myy=myy+1
              itab(kstf+nec)=idet(myy)
            end do
          end if
      i=nec+1
c                                                                       
c                                                                       
c                                                                       
      go to (720,110,120,130,140,150,160,170,180) ,i                    
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
 110  continue                                                          
      km=km+1                                                           
      ntm(km)=itab(1)                                                   
      npm(km)=itab(2)                                                   
      go to 700                                                         
 120  continue                                                          
      do 220 i=1,2                                                      
      kd=kd+1                                                           
      ntd(kd)=itab(i)                                                   
      npd(kd)=itab(i+2)                                                 
 220  continue                                                          
      go to 700                                                         
 130  continue                                                          
      do 230 i=1,3                                                      
      kt=kt+1                                                           
      ntt(kt)=itab(i)                                                   
      npt(kt)=itab(i+3)                                                 
 230  continue                                                          
      go to 700                                                         
 140  continue                                                          
      do 240 i=1,4                                                      
      kq=kq+1                                                           
      ntq(kq)=itab(i)                                                   
      npq(kq)=itab(i+4)                                                 
 240  continue                                                          
      go to 700                                                         
 150  continue                                                          
      do 250 i=1,5                                                      
      kp=kp+1                                                           
      ntp(kp)=itab(i)                                                   
      npp(kp)=itab(i+5)                                                 
 250  continue                                                          
      go to 700                                                         
 160  continue                                                          
      do 260 i=1,6                                                      
      kh=kh+1                                                           
      nth(kh)=itab(i)                                                   
      nph(kh)=itab(i+6)                                                 
 260  continue                                                          
      go to 700                                                         
 170  continue                                                          
      do 270 i=1,7                                                      
      k7=k7+1                                                           
      nt7(k7)=itab(i)                                                   
      np7(k7)=itab(i+7)                                                 
 270  continue                                                          
      go to 700                                                         
 180  continue                                                          
      do 280 i=1,8                                                      
      ko=ko+1                                                           
      nto(ko)=itab(i)                                                   
      npo(ko)=itab(i+8)                                                 
 280  continue                                                          
 700  continue                                                          
       nbdet(nec+1)=nbdet(nec+1)+1                                      
        end if
      end do
      go to 600                                                         
 3487 continue                                                          
      rewind 60                                                         
      write(6,*) ' determinants lus: ',ndet                             
      return                                                            
      end                                                               
