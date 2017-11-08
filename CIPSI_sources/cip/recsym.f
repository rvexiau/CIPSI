      subroutine recsym(m,sym)                                          
      implicit real*8(a-h,o-x,z),logical*1 (y)                          
      include 'pshf.prm'
      integer*4 itm                                                
      common/nature/itm(metz,nsymz)                                     
      common/cipsi/escf,test,teseff,tdet,twdet,ndiag,norb,mnorb,
     1             noca,nocb,noc1,noc2,igela,igelb,ncf,nsym,isym,isz,
     2             ntrsy,metat,lecdet,maxdet,yprt,yion,ypertu,ybrd,yen,
     3             ywvi,ywvf,ytoul,ystkic,ywdeg
      character*3 sym                                                   
      character*7 tabsym                                                
      character*1 pm,tab                                                
      dimension tabsym(10),tab(9)                                       
      data tab/'1','2','3','4','5','6','7','8','9'/                     
      call initsm(tabsym,nsym)                                          
      itab=isym                                                         
      pm=' '                                                            
      if (nsym.eq.10) then                                              
         if (isym.eq.1.or.isym.eq.2) then                               
            itrnsf=itm(m,2)                                             
            pm='+'                                                      
            if (itrnsf.eq.-1) then                                      
               itab=isym+8                                              
               pm=' '                                                   
            endif                                                       
         else if (isym.eq.9.or.isym.eq.10) then                         
            itrnsf=itm(m,ntrsy+1)                                       
            if (itrnsf.eq.-1) then                                      
               itab=isym-8                                              
               pm='-'                                                   
            endif                                                       
         endif                                                          
         irs=2                                                          
      else if (nsym.eq.5) then                                          
         if (isym.eq.1) then                                            
            itrnsf=itm(m,2)                                             
            pm='+'                                                      
            if (itrnsf.eq.-1) then                                      
               itab=isym+4                                              
               pm=' '                                                   
            endif                                                       
         else if (isym.eq.5) then                                       
            itrnsf=itm(m,4)                                             
            if (itrnsf.eq.-1) then                                      
               itab=isym-4                                              
               pm='-'                                                   
            endif                                                       
         endif                                                          
         irs=1                                                          
      else if (nsym.eq.2.or.nsym.eq.4.or.nsym.eq.8) then                
         irs=3                                                          
      else                                                              
         irs=1                                                          
         itab=isym                                                      
         tabsym(itab)=tab(itab)                                         
      endif                                                             
      sym=tabsym(itab)(1:irs)//pm                                       
      return                                                            
      end                                                               
