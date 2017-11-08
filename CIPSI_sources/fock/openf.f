      subroutine openf
      character*80 prefix
      character*10 f25,ijkl,f50,ener
      character*90 f251,f501,ijkl1,blan90,ener1
      namelist /fokfil/ prefix,f25,ijkl,f50,ener,
     * f251,ijkl1,f501,ener1
      blan90=' '
      f25='f25       '
      ijkl='ijkl     '
      ener='ener     '
      f50='f50       '
      prefix= '          '
      f251=blan90
      f501=blan90
      ijkl1=blan90
      ener1=blan90
      read(5,fokfil)
      write(6,fokfil)
      call nomfil(25,f25,f251,prefix,'old       ')
      call nomfil(40,ijkl,ijkl1,prefix,'new       ')
      call nomfil(50,f50,f501,prefix,'old       ')
      call nomfil(33,ener,ener1,prefix,'unknown   `')
      return
      end
