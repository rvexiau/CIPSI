      subroutine openf
      logical*1 yuni
      character*80 prefix
      character*20 bdveca,bdvecb,f04a,f04b,ijkl,ijom,info,detcla,detclb,
     *detspina,detspinb,omg,ona,dip,dmat,cv,ener
      character*90 bdveca1,bdvecb1,f04a1,f04b1,ijkl1,ijom1,info1,
     *detcla1,detclb1,detspina1,detspinb1,cv1,
     *omg1,ona1,dip1,dmat1,blan90,ener1
      common/uni/yuni
      namelist/rofil/prefix,bdveca,bdvecb,f04a,f04b,ijkl,ijom,info,
     * detcla,detclb,detspina,detspinb,cv,ener,ener1,
     *omg,ona,dip,dmat,
     *bdveca1,bdvecb1,f04a1,f04b1,ijkl1,ijom1,info1,
     *detcla1,detclb1,cv1,
     *omg1,ona1,dip1,dmat1
      blan90=' '                                                        
      bdveca='bdvec'
      bdvecb='bdvec'
      f04a='f04'
      f04b='f04'
      ijkl='ijkl'
      ijom='ijom'
      info='info'
      detcla='det_cl'
      detclb='det_cl'
      detspina='det_spin'
      detspinb='det_spin'      
      omg='om'
      ona='ona'
      dip='dip'
      dmat='dmat'
      cv='vcpp'
      ener='ener'
      prefix= '          '                                              
      bdveca1=blan90
      bdvecb1=blan90
      f04a1=blan90
      f04b1=blan90
      ijkl1=blan90
      ijom1=blan90
      info1=blan90
      detcla1=blan90
      detclb1=blan90
      detspina1=blan90
      detspinb1=blan90      
      omg1=blan90
      ona1=blan90
      dip1=blan90
      dmat1=blan90
      cv1=blan90
      ener1=blan90
      read(5,rofil) 
      write(6,rofil)                                                   
      call nomfil(1,dmat,dmat1,prefix,'unknown   ')
      call nomfil(2,info,info1,prefix,'old       ')
      call nomfil(10,ijom,ijom1,prefix,'unknown   ')
      call nomfil(11,omg,omg1,prefix,'unknown   ')
      call nomfil(12,ona,ona1,prefix,'unknown   ')
      call nomfil(18,dip,dip1,prefix,'unknown   ')
      call nomfil(40,ijkl,ijkl1,prefix,'old       ')
      call nomfil(3,f04a,f04a1,prefix,'old       ')
      call nomfil(58,bdveca,bdveca1,prefix,'unknown   ')
      call nomfil(61,detcla,detcla1,prefix,'unknown   ')
      call nomfil(64,detspina,detspina1,prefix,'unknown   ')           
      call nomfil(17,cv,cv1,prefix,'unknown   ') 
      call nomfil(33,ener,ener1,prefix,'old       ') 
      if((f04a.eq.f04b).or.(bdveca.eq.bdvecb).or.(detcla.eq.detclb))then
         yuni=.true.
         write(6,*) 'f04b = f04a' 
         write(6,*) 'bdvecb = bdveca' 
         write(6,*) 'detclb = detcla' 
         write(6,*) 'detspinb = detspina'          
      else
         yuni=.false.
         call nomfil(4,f04b,f04b1,prefix,'old       ')
         call nomfil(59,bdvecb,bdvecb1,prefix,'unknown   ')
         call nomfil(62,detclb,detclb1,prefix,'unknown   ')
         call nomfil(65,detspinb,detspinb1,prefix,'unknown   ')            
      endif
      return                                                            
      end                                                               
