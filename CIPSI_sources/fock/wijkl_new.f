      subroutine wijkl_new(bijkls,bijkld,ydp)                           
      implicit real*8 (a-h,o-x,z),logical*1 (y)                         
      include 'pshf.prm'
      integer*2 ij0,ipair                                              
      common/jkf/ f(doa*(doa+1)/2),aj(doa*(doa+1)/2),           
     1 ak(doa*(doa+1)/2)                                            
      dimension bufd(lsize)
      real*4 bufs(lsize)
      dimension bijkld(*)
      real*4 bijkls(*)
c *** table de multiplication des groupes                               
c ***   its(nsymz,nsymz)                                                
c ***                                                                   
      common/sy/itsym(doa),its(nsymz,nsymz)         
      common/nodim/norb,noca,nocb,igel,ncf,isym,nsym,isz,ntrsy
      dimension ipair(doa*doa,nsymz,2)           
      dimension ntt(nsymz),ntp(nsymz),npp(nsymz)                        
      dimension ntttt(nsymz)
      dimension nt(nsymz),np(nsymz),ij0(doa,doa)         
      noc=noca
      do is=1,nsym                                                  
      nt(is)=0                                                          
      np(is)=0                                                          
      end do
      do 150 i=1,noc                                                    
      is=itsym(i)                                                       
      nt(is)=nt(is)+1                                                   
      ij0(is,nt(is))=i
 150  continue
      do 250 i=noc+1,norb                                               
      is=itsym(i)                                                       
      np(is)=np(is)+1                                                   
      ij0(is,np(is)+noc)=i
 250  continue
c ***                                                                   
c *** trous et particules donnant une distribution de symetrie is       
c ***                                                                   
      do 440 is=1,nsym                                                  
      ntt(is)=0                                                         
      npp(is)=0                                                         
      ntp(is)=0                                                         
      nn=0                                                              
      do 360 ks=1,nsym                                                  
      do 360 ls=ks,nsym                                                 
      if(its(ks,ls).ne.is) go to 360                                    
c *** trou-trou                                                         
      if(nt(ks).eq.0.or.nt(ls).eq.0) go to 360
      do 350 ixx=1,nt(ks)                                          
      i=ij0(ks,ixx)
      if(ks.eq.ls) then                                                 
         lim=ixx                                                       
         else                                                           
         lim=1                                                    
      endif                                                             
      do 340 jxx=lim,nt(ls)                                             
      j=ij0(ls,jxx)
      nn=nn+1                                                           
      ipair(nn,is,1)=i                                                  
 340  ipair(nn,is,2)=j                                                  
 350  continue                                                          
 360  continue                                                          
      ntt(is)=nn                                                        
c *** trou-part                                                         
      do 375 ks=1,nsym                                                  
      do 375 ls=1,nsym                                                  
      if(its(ks,ls).ne.is) go to 375                                    
      if(nt(ks).eq.0.or.np(ls).eq.0)go to 375
      do 370 ixx=1,nt(ks)                                          
      i=ij0(ks,ixx)
      do 365 jxx=noc+1,noc+np(ls)                                       
      j=ij0(ls,jxx)
      nn=nn+1                                                           
      ipair(nn,is,1)=i                                                  
      ipair(nn,is,2)=j                                                  
 365  continue                                                          
 370  continue                                                          
 375  continue                                                          
      ntp(is)=nn                                                        
c *** part-part                                                         
      do 395 ks=1,nsym                                                  
      do 395 ls=ks,nsym                                                 
      if(its(ks,ls).ne.is) go to 395                                    
      if(np(ks).eq.0.or.np(ls).eq.0)go to 395
      do 390 ixx=noc+1,noc+np(ks)                                       
      i=ij0(ks,ixx)
      if(ks.eq.ls) then                                                 
         lim=ixx                                                        
         else                                                           
         lim=noc+1                                                    
      endif                                                             
      do 380 jxx=lim,noc+np(ls)                                             
      j=ij0(ls,jxx)
      nn=nn+1                                                           
      ipair(nn,is,1)=i                                                  
 380  ipair(nn,is,2)=j                                                  
 390  continue                                                          
 395  continue                                                          
      npp(is)=nn                                                        
 400  continue                                                          
c      nt(is)=nn                                                        
 440  continue                                                          
      do 650 is=1,nsym                                                  
      do 650 ij=1,npp(is)                                               
      i=ipair(ij,is,1)                                                  
      j=ipair(ij,is,2)                                                  
      if(i.le.j) then                                                   
         ij0(i,j)=ij                                                    
         if(i.ne.j) ij0(j,i)=is                                         
         else                                                           
         ij0(j,i)=ij                                                    
         ij0(i,j)=is                                                    
      endif                                                             
 650  continue                                                          
      write(40) ((ij0(i,j),i=1,norb),j=1,norb)                          
c ***                                                                   
c *** integrales                                                        
c ***                                                                   
c *** 4 trous                                                           
c ***                                                                   
      nint=1                                                               
      do is=1,nsym                                                  
        ntttt(is)=nint                                                 
        do ij=1,npp(is)                                               
          i=ipair(ij,is,1)                                              
          j=ipair(ij,is,2)                                              
          do kl=1,ij                                                 
            k=ipair(kl,is,1)                                         
            l=ipair(kl,is,2)                                          
            nint=nint+1                                                 
          end do
        end do
      end do
      nint=nint-1
      write(40)nint,(ntttt(i),i=1,nsym)
      n=1                                                               
      do is=1,nsym                                                  
        do ij=1,npp(is)                                               
          i=ipair(ij,is,1)                                              
          j=ipair(ij,is,2)                                              
          do kl=1,ij                                                 
            k=ipair(kl,is,1)                                         
            l=ipair(kl,is,2)                                          
            if(ydp)then
	      bufd(n)=ai(i,j,k,l,bijkls,bijkld,ydp)                    
	      if(n.eq.lsize)then
		write(40)bufd
		n=0
              end if
            else
	      bufs(n)=ai(i,j,k,l,bijkls,bijkld,ydp)                    
	      if(n.eq.lsize)then
		write(40)bufs
		n=0
              end if
	    end if
            n=n+1                                                    
          end do
        end do
      end do
 500  continue                                                          
      if(ydp)then
	write(40)bufd
      else
	write(40)bufs
      end if
      write(6,5000)nint,nsym                                            
      do 6000 is=1,nsym                                                 
      write(6,5001)is,ntttt(is)                                         
 6000 continue                                                          
 5000 format(/,' nombre d integrales',i9,                               
     1       /,' nombre de symetries',i9)                               
 5001 format(i5,6x,i8)                                                  
      return                                                            
      end                                                               
