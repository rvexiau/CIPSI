      subroutine openf                                                  
      character*80 prefix
      character*20 f04,det_cl,mat,bdvec,ener,det_spin
      character*90 f041,blan90,det_cl1,mat1,bdvec1,ener1,
     * det_spin1
      namelist /bdfil/ prefix,f04,f041,det_spin,det_spin1,
     * det_cl,mat,det_cl1,mat1,ener,ener1,bdvec,bdvec1
      blan90=' '                                                        
      f04='f04       '
      mat='mat       '
      det_cl='det_cl    '
      det_spin='dspin    '
      bdvec='bdvec     '
      ener='ener      '
      prefix= '          '                                              
      f041=blan90                                                       
      bdvec1=blan90
      det_cl1=blan90
      ener1=blan90
      mat1=blan90
      det_spin1=blan90
      read(5,bdfil)        
      write(6,bdfil)                  
      call nomfil(1,mat,mat1,prefix,'unknown   ')
      call nomfil(4,f04,f041,prefix,'old       ')                          
      call nomfil(33,ener,ener1,prefix,'unknown   ')
      call nomfil(59,bdvec,bdvec1,prefix,'unknown   ')
      call nomfil(61,det_cl,det_cl1,prefix,'unknown   ')
      return                                                            
      end                                                               
