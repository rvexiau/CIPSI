      subroutine psgrin
      implicit real*8(a-h,o-z)
      character*40 typegp,typf21,fin_de_fichier
      logical egp,ok
      common/psgrx/negp,egp,typegp(50)
      common/iofile/ ir,iw,ip,is
      common/output/nprint
      namelist/egpin/negp,typegp
      fin_de_fichier='fin de fichier'
      read(5,egpin)
      if(nprint.ne.4) 
     < write(iw,*)' number of effective group potentials',negp
      do i=1,negp
	rewind 21
    1	read(21)typf21
	if(typf21.eq.fin_de_fichier)go to 100
	if(typf21.eq.typegp(i))then
	  ok=.true.
          go to 100
        else
	  read(21)
	  read(21)
	  go to 1
        end if
  100   if(.not.ok)then
	  write(iw,*)' effective group potential',typegp(i),
     < ' is missing on file 21  program stops'
	  stop
        end if
      end do
      return   
      end
