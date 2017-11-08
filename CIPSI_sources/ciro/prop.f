      subroutine prop(ydipol,ypopul,typdip)                            
c---------------------------------------------------------------------  
c   namelist propi                                                      
c                                                                       
c   ce programe calcule les proprietes moleculaires suivantes:          
c    apres pshf,cip120 et ro120                                         
c   1.- la population electronique de mulliken (ypopul=t)               
c   2.- les moments dipolaires                 (ydipol=t)               
c                                                                       
c  nprint:0 resultat seul (defaut)                                      
c         1 impresion  des matrices x y z dans la base atomique         
c         2 impresion  des matrices densite                             
c         3 impresion  de tout                                          
c------------------------------------------------------------------     
      implicit real*8 (a-h,o-x,z)                                       
      implicit logical*1 (y)                                            
      logical last,out                                                  
      character*40 typdip
      character*4 iflab                                                 
      include 'pshf.prm'
      common/cd/r(ncouz*(doa**2)),ndd,ndt,ncf1,ncf2,metat,icf2,ncou,  
     &ibeg(ncouz),ietat(ncouz),jetat(ncouz),ycou(metz,metz),      
     &yiet(metz),yjet(metz)                                             
      common/iofile/ir,iw,ip,is,iq,ih,iv,if2,if3,if4,ioda(19),title(10) 
      common/output/nprint,itol,icut,normf,normp                        
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,zan(dc),c(3,dc)          
      common/nshel/ex(dgp),cs(dgp),cp(dgp),cd(dgp),cf(dgp),cg(dgp),     
     1kstart(ds),katom(ds),ktype(ds),kng(ds),kloc(ds),        
     2 kmin(ds),kmax(ds),nshell,kmaz(ds)                          
      common/runlab/iflab(3,doa)                                      
      common/atlim/limlow(dc),limsup(dc)                                
      common/infpot/iatno(dc)                   
      dimension mmin(10),mmax(10),mmaz(10)                              
      dimension s(doa*(doa+1)/2)                                    
c      dimension e(dc*(dc+1)/2),gac(dc,2)                                        
      dimension a(dc),ns(dc),ks(dc)                                     
      dimension newsh(ds,48),ptr(3,144),dtr(6,288),ftr(10,480)          
      data mmin/1,2,5,11,21,1,1,5,11,21/                                
      data mmax/1,4,10,20,35,1,4,10,20,35/                              
      data mmaz/1,4,9,17,30,1,4,10,20,35/                               
      data pi32/5.5683279968317d0/                                      
      nprint=0                                                          
      ir=5                                                              
      iw=6                                                              
      is=8                                                              
      ih=10                                                             
      if2=2                                                             
      itol=9
      icut=9                                                            
      normf=0                                                           
      normp=0                                                           
      npsd=0                                                            
 9000 continue                                                          
      read(2) ioda,nt,nshell,num,nveca,nvecb,bname,ngauss,              
     1ne,ich,mul,na,nb,nat,group,naxis,title,                           
     2(ktype(k),kstart(k),kloc(k),katom(k),kng(k),k=1,nshell),          
     3(ex(k),cs(k),cp(k),cd(k),cf(k),cg(k),k=1,ngauss),                 
     4((c(i,k),i=1,3),a(k),ns(k),ks(k),iatno(k),zan(k),k=1,nat),        
     5((newsh(i,k),i=1,nshell),k=1,nt),((ptr(i,k),i=1,3),k=1,3*nt),     
     6((dtr(i,k),i=1,6),k=1,6*nt),((ftr(i,k),i=1,10),k=1,10*nt)         
c    7                ,((gtr(i,k),i=1,15),k=1,15*nt)                    
     8  ,((iflab(i,k),i=1,3),k=1,num)                                   
      nx=(num*(num+1))/2                                                
      do 1926 i=1,nshell                                                
         kmin(i)=mmin(ktype(i))                                         
         kmax(i)=mmax(ktype(i))                                         
         kmaz(i)=mmaz(ktype(i))                                         
1926  continue                                                          
      if(ypopul)then                                                    
      write(6,*)'*****************************'                         
      write(6,*) ' analyse de population'                               
      write(6,*)'*****************************'                         
      call mulken                                                       
      endif                                                             
      if(ydipol) then                                                   
      write(6,*)'*****************************'                         
      write(6,*) ' moments dipolaires '                                 
      write(6,*)'*****************************'                         
      call dipole(typdip)                                               
      endif                                                             
      return                                                              
      end                                                               
