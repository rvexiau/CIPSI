      subroutine openf                                                  
      character*80 prefix
      character*10 f04ref,som,f04,bdvec,ener,bdvecr,det,detr,diavec
      character*90 blan90,f041,f04ref1,bdvec1,ener1,som1,
     * det1,detr1,bdvecr1,diavec1            
      namelist/hefil/prefix,f04,f041,som,som1,f04ref,f04ref1,
     *bdvec,bdvec1,ener,ener1,det,det1,detr,detr1,bdvecr,bdvecr1,
     *diavec,diavec1
      blan90=' '                                                        
      f04='f04       '
      f04ref='f04ref    '
      som='som      '
      bdvec='bdvec    '
      bdvecr='bdvecr   '
      diavec='diavec   '
      det='det      '
      detr='detr     '
      ener='ener     '
      prefix= '          '                                              
      som1=blan90                                                       
      f041=blan90                                                       
      f04ref1=blan90                                               
      bdvec1=blan90
      bdvecr1=blan90
      diavec1=blan90
      det1=blan90
      detr1=blan90
      ener1=blan90
      write(6,*) 'op hefil'
      read(5,hefil)                                                    
      write(6,*) 'finlec hefil'
      write(6,hefil)                                                   
      call nomfil(12,som,som1,prefix,'unknown   ')                         
      call nomfil(4,f04,f041,prefix,'unknown   ')                          
      call nomfil(3,f04ref,f04ref1,prefix,'unknown   ')   
      call nomfil(58,bdvec,bdvec1,prefix,'unknown   ')
      call nomfil(56,diavec,diavec1,prefix,'new       ')
      call nomfil(60,det,det1,prefix,'unknown   ')
      call nomfil(59,bdvecr,bdvecr1,prefix,'unknown   ')
      call nomfil(61,detr,detr1,prefix,'unknown   ')
      call nomfil(33,ener,ener1,prefix,'old       ')
      return                                                            
      end                                                               
