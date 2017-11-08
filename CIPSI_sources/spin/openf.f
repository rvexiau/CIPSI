      subroutine openf                                                  
      character*80 prefix
      character*20 f04,det_cl,mat,bdvec,ener,det_spin,ijkl
      character*90 f041,blan90,det_cl1,mat1,bdvec1,ener1,ijkl1
      character*90 det_spin1
      namelist /bdfil/ prefix,f04,f041,
     * det_cl,mat,det_cl1,mat1,ener,ener1,bdvec,bdvec1,
     * det_spin,det_spin1
      blan90=' '                                                                                                        
      bdvec='bdvec    '     
      det_spin='det_spin    ' 
      ijkl='ijkl     '   
      ijkl1=blan90            
      prefix= '          '                                                                                                                                                    
      bdvec1=blan90
      det_spin1=blan90      
      read(5,bdfil)                                                                                                                                                  
      call nomfil(59,bdvec,bdvec1,prefix,'unknown   ')    
      call nomfil(64,det_spin,det_spin1,prefix,'unknown   ')
      call nomfil(40,ijkl,ijkl1,prefix,'old       ')         
      return                                                            
      end                                                               
