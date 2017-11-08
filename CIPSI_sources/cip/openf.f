      subroutine openf                                                  
      character*80 prefix
      character*10 det,ijkl,f04,ocno,bdvec,ener
      character*90 det1,f041,ijkl1,blan90,ocno1,bdvec1,ener1            
      namelist /cipfil/ prefix,det,ijkl,f04,det1,ijkl1,f041,ocno1,
     *bdvec,bdvec1,ocno,ener,ener1
      blan90=' '                                                        
      det='det       '                                                  
      ijkl='ijkl     '                                                  
      f04='f04       '
      ocno='ocno     '
      bdvec='bdvec    '
      ener='ener     '
      prefix= '          '                                              
      det1=blan90                                                       
      f041=blan90                                                       
      ijkl1=blan90                                               
      ocno1=blan90
      bdvec1=blan90
      ener1=blan90
      read(5,cipfil)                                                    
      write(6,cipfil)                                                   
      call nomfil(60,det,det1,prefix,'unknown   ')                         
      call nomfil(40,ijkl,ijkl1,prefix,'old       ')                           
      call nomfil(4,f04,f041,prefix,'unknown   ')                          
      call nomfil(21,ocno,ocno1,prefix,'unknown   ')
      call nomfil(59,bdvec,bdvec1,prefix,'unknown   ')
      call nomfil(33,ener,ener1,prefix,'old       ')
      return                                                            
      end                                                               
