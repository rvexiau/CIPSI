      subroutine lec60                                                  
     *(ntm,npm,ntd,npd,ntt,npt,ntq,npq,ntp,npp,nth,nph) 
c                                                                       
c  lecture des determinants sur la file 60                              
c                                                                       
c    tau :seuil de selection                                                                 
c                                                                       
      implicit real*8(a-h,o-p,r-x,z),logical*1(y)                       
      include 'pshf.prm'
      integer*4 ntm(ndimh),npm(ndimh),ntd(2*ndimh),npd(2*ndimh),        
     *ntt(3*ndimh),npt(3*ndimh),ntq(4*ndimh),npq(4*ndimh)               
     *,ntp(5*ndimh),npp(5*ndimh),nth(6*ndimh),nph(6*ndimh)                          
      common/kexc/ kf,km,kd,kt,kq,kp,kh                   
      common/e0cip/ebmp(metz),eben(metz),tau      
      dimension idot(16*lungo),necst(lungo),vs(lungo)                                           
      dimension itab(40)                         
c                                                                                                                               
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
                                                     
600   read(60,end=3487) mw,nw,(idot(i),i=1,mw),(necst(i),i=1,nw),       
     &(vs(i),i=1,nw)                                                       
c                                                                       
c                                                                       
c                                                                       
c                                                                       
c             
      myy=0
      do  lstf=1,nw
        nec=necst(lstf)
        if(nec.gt.6) then
          write(6,*) 'erreur degre d excitation',nec,' depasse les ',
     *    'possibilites du programme'
          stop
        end if
        if(vs(lstf).lt.tau) then
	  myy=myy+2*nec
	  else
          ndet=ndet+1                              
          if(nec.ne.0)then
	    do  kstf=1,nec
	      myy=myy+1
              itab(kstf)=idot(myy)
	      myy=myy+1
              itab(kstf+nec)=idot(myy)
            end do
          end if
      i=nec+1
c                                                                       
c                                                                       
c                                                                       
      go to (720,110,120,130,140,150,160) ,i                    
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
 700  continue                                                          
                                  
        end if
      end do
      go to 600                                                         
 3487 continue                                                          
      rewind 60                                                         
      write(6,*) ' determinants lus: ',ndet                             
      return                                                            
      end                                                               
