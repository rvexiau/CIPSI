      subroutine rozero                                                 
c      calcul de la matrice densite a l ordre zero                      
c      entres etats                                                     
      implicit real*8(a-h,o-x,z),logical*1(y)                      
 
      include 'pshf.prm'                  
      integer*4 ne,nd,iorb,ispin,itsym,its                    
      integer*4 trou,part                  
      integer*4 nedum,nddum
cRG
      integer*4 trdum,pardum                                
cRG      
! ### DB : vv( :, : ) contient les vecteurs propres pour impression en clair
!      real*8 vv( metz, ndetz )
! ###      
      common/comc/c(2*ndimh*metz)                                       
      common/ic/e(metz2),trou(8*ndimh),part(8*ndimh),
     & ne(2*ndimh),nd(2*ndimh),ispin(2*doa), 
     & iorb(2*doa),itsym(2*doa),its(nsymz,nsymz),igels(2),icount,
     & norb,noca,nocb,mnorb,ncf,nsym,isym,ntrsy,metat1,metat2,
     & yoc(2*doa),yprt
      common/uni/yuni
      common/cd/r(ncouz*(doa**2)),ndd,ndt,ncf1,ncf2,metat,icf2,ncou,  
     *ibeg(ncouz),ietat(ncouz),jetat(ncouz),ycou(metz,metz),      
     &yiet(metz),yjet(metz)                                             
      common/fimmif/n1fic,n2fic,n1fim,n2fim,ymoyen                      
      common/gel/nao      
     
!!$  >>> debug <<<
ccc        write(*,*) 'beginning rozero' 
ccc        write(*,*) 'nao = ',nao
!!$  >>> debug <<<      
c      ymoyen=.true.
      if (.not.ymoyen) then                                             
         rewind n1fic                                                   
         read(n1fic)nnn,nnn,nnn,mnn,nncf1,nnn,ydum,ii,i,                
     &   (nedum,nddum,j=1,nncf1),(trdum,pardum,j=1,ii),(c(j),j=1,i),
     *   (e(j),j=1,mnn)     
         rewind n2fic                                                   
         read(n2fic)nnn,nnn,nnn,mnn,nncf2,nnn,ydum,iiii,iii             
     &   ,(nedum,nddum,j=1,nncf2)                                       
     &   ,(trdum,pardum,j=1,iiii)                                       
     &   ,(c(j+metz*ndimh),j=1,iii),                                   
     *   (e(metz+j),j=1,mnn)     
      else                                                              
         rewind n1fim                                                   
         read(n1fim) nncf1,mmetat ,(c(i),i=1,nncf1*mmetat )
!#################         
! DB impression des vecteurs propres en clair
!          vv( 1 : mmetat, 1 : nncf1 ) =+
!      +reshape( c( 1:nncf1*mmetat ),+
!      +(/ mmetat, nncf1 /), ORDER=(/ 2, 1 /) )
!      
!           write( 236, '(<mmetat>(2X, f9.6))' )+
!      +vv( 1 : mmetat, 1 : nncf1 )
!##################         
         rewind n2fim                                                   
         read(n2fim) nncf2,mmetat ,(c(i+metz*ndimh),i=1,nncf2*mmetat )
      endif                                                             
!!$  >>> debug <<<
ccc        write(*,*) 'flag1' 
ccc        write(*,*) 'nao = ',nao
ccc           write(*,*) 'ideb =',ideb
ccc        write(*,*) 'ncf1 = ',ncf1,ideb,ncf
!!$  >>> debug <<<      
      ideb=ncf1+1                                                       
      kkk=0                                                             
      icf2=metz*ndimh
      do 5 k1=1,ncf1                                                    
         do ijk=1,norb*2                                                
            yoc(ijk)=.false.                                            
         enddo                                                          
         nek1=ne(k1)                                                    
         if (nek1.eq.0) goto 1067                                       
         ndi =nd(k1)                                                    
         do ijk=1,nek1                                                  
            nt=trou(ndi+ijk)                                            
            np=part(ndi+ijk)                                            
            yoc(nt)=.true.                                              
            yoc(np)=.true.                                              
         enddo                                                          
1067     continue                                                       
!!$  >>> debug <<<
ccc        write(*,*) 'inside the 5 loop' 
ccc        write(*,*) 'k1/nao = ',k1,nao,norb
ccc        if(nao/=99) stop
!!$  >>> debug <<<      
      do 6 k2=ideb,ncf                      
      nek2=ne(k2)                                                       
      kkk=kkk+1                                                         
      ndif=nek1-nek2           
c      if((nek1.eq.2).and.(nek2.eq.2)) goto 6
      if((ndif.gt.1).or.(ndif.lt.-1)) go to 6                           
      ndi=nd(k2)+1                                                      
      nei=nek2          
8     if(nei.eq.0) go to 88                                             
      nfi=ndi+nei-1           
      do 300 i=ndi,nfi     
      m=trou(i)                                                         
      if(.not.yoc(m)) ndif=ndif+1                                       
      m=part(i)                                                         
      if(.not.yoc(m)) ndif=ndif+1                                       
300   continue                             
      if(ndif.gt.1) go to 6                                             
88    if(ndif.eq.0) go to 800                                           
      call rntd(k1,k2,ndif,mu,mv,isign)      
      mu=iorb(mu)                                                       
      mv=iorb(mv)                                                       
      id=(mu-1)*norb+mv                                                 
      kk2=k2-ncf1                                                       
      do 65 icou=1,ncou           
      m1=ietat(icou)                                                    
      c1=c(ncf1*(m1-1)+k1)                                              
      m2=jetat(icou)                                                    
      c2=c(icf2+(m2-1)*ncf2+kk2)                                
      ind=ibeg(icou)                                                   
      indix=ind+id                                                      
65    r(indix)=r(indix)+c1*c2*isign                                     
      go to 6                                                           
800   if(nei.eq.0) go to 6                                              
      do 810 i=ndi,nfi                                                  
      mu=trou(i)                                                        
      mu=iorb(mu)                                                       
      mv=part(i)                                                        
      mv=iorb(mv)                                                       
      idmu=(mu-1)*norb+mu                                               
      idmv=(mv-1)*norb+mv                                               
      kk2=k2-ncf1                                                       
c      km1=-ncf1                                                         
       kk1=-ncf1                                                         
      do 865 icou=1,ncou   
      m1=ietat(icou)                                                    
      c1=c(ncf1*(m1-1)+k1)                                              
      m2=jetat(icou)                                                    
      c2=c(icf2+ncf2*(m2-1)+kk2)                                
      ind=ibeg(icou)                                                   
      indix=ind+idmu                                                    
      cccc=c1*c2                                                        
      r(indix)=r(indix)-cccc                                            
      if(mv.gt.0) then
        index=ind+idmv                                                    
        r(index)=r(index)+cccc                                            
      endif
865   continue                                                          
810   continue                                                          
!!$  >>> debug <<<
ccc        write(*,*) 'inside the 6 loop' 
ccc        write(*,*) 'k1/nao = ',k1,nao,norb
ccc        if(nao/=99) stop
!!$  >>> debug <<<      
6     continue                                                          
5     continue                                                          
c calcul de la trace de ro
!!$  >>> debug <<<
ccc        write(*,*) 'flag2' 
ccc        write(*,*) 'nao = ',nao
!!$  >>> debug <<<      
      do icou=1,ncou
       traro=0.d0
       do io=1,norb
        traro=traro+r(ibeg(icou)+1+(io-1)*(norb+1))
       enddo
c       write(6,*)'dans rozero,icou,ietat,jetat,tracede ro',icou,
c     *            ietat(icou),jetat(icou),traro
      enddo
!!$  >>> debug <<<
ccc        write(*,*) 'Ending rozero' 
ccc        write(*,*) 'nao = ',nao
!!$  >>> debug <<<      
      return                                                            
      end                                                               
