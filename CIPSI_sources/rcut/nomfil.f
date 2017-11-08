      subroutine nomfil(iunit,racine,file1,prefix,status)
      character*1 prefix(80),racine(10)
      character*10 status
      character*90 file
      character*90 blan90,file1
      character*1  fil(90)
      equivalence (fil(1),file)
      blan90=' '
      file='                                           '
C      OPEN FILES
      call nomf(fil,prefix,racine)
      if(file1.ne.blan90)file=file1
      write(6,*)'ouverture de la file',iunit
      write(6,*)fil,iunit,racine,status
      OPEN (UNIT=iunit,FORM='UNFORMATTED',STATUS=status,file=file)
      return
      end
