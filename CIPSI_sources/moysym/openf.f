      subroutine openf                                                  
      character*80 prefix
      character*20 det,ijkl,f04,det_cl,det_spin,csf,mat,enmp,ener
      character*20 mats1,mats2,mats3,mats4,mats5,mats6
      character*90 mats1b,mats2b,mats3b,mats4b,mats5b,mats6b      
      character*90 det1,f041,ijkl1,blan90,det_cl1,mat1,enmp1,
     * csf1,ener1,det_spin1            
      namelist/moyfil/prefix,det,ijkl,f04,det1,ijkl1,f041,
     * det_cl,det_spin,mat,enmp,csf,det_cl1,mat1,enmp1,
     * ener,ener1,csf1,det_spin1
      blan90=' '                                                        
      det='det       '                                                  
      ijkl='ijkl     '                                                  
      f04='f04       '
      mat='mat       '
      mats1='mat1    '
      mats2='mat2    '
      mats3='mat3    '
      mats4='mat4    '
      mats5='mat5    '
      mats6='mat6    '
      det_cl='det_cl    '
      det_spin='det_spin    '      
      csf='csf       '
      enmp='enmp      '
      ener='ener      '
      prefix= '          '                                              
      det1=blan90                                                       
      f041=blan90                                                       
      ijkl1=blan90                                               
      det_cl1=blan90
      csf1=blan90
      enmp1=blan90
      ener1=blan90
      mat1=blan90
      det_spin1=blan90
      mats1b=blan90
      mats2b=blan90
      mats3b=blan90
      mats4b=blan90
      mats5b=blan90
      mats6b=blan90      
      read(5,moyfil)                                                    
      write(6,moyfil)      
      call nomfil(70,mat,mat1,prefix,'unknown   ')
      call nomfil(4,f04,f041,prefix,'old       ')                          
      call nomfil(33,ener,ener1,prefix,'old       ')                           
      call nomfil(40,ijkl,ijkl1,prefix,'old       ')                           
      call nomfil(60,det,det1,prefix,'old       ')                         
      call nomfil(61,det_cl,det_cl1,prefix,'unknown   ')
      call nomfil(62,csf,csf1,prefix,'unknown   ')
      call nomfil(63,enmp,enmp1,prefix,'unknown   ')
      call nomfil(64,det_spin,det_spin1,prefix,'unknown   ')   

c     open a H matrix file for each value of 2S+1 
      call nomfil(71,mats1,mats1b,prefix,'unknown   ')
      call nomfil(72,mats2,mats2b,prefix,'unknown   ')
      call nomfil(73,mats3,mats3b,prefix,'unknown   ')
      call nomfil(74,mats4,mats4b,prefix,'unknown   ')
      call nomfil(75,mats5,mats5b,prefix,'unknown   ')
      call nomfil(76,mats6,mats6b,prefix,'unknown   ')
                                          
      return                                                            
      end                                                               
