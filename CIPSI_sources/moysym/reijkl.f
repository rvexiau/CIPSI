      subroutine reijkl(norb,nsym,ydp)                    
      implicit real*8 (a-h,o-x,z),logical*1 (y)                         
      include 'pshf.prm'
      integer*2 ij0
      real*4 bijkl,bufs                                                
      common/ij1/ij0(doa,doa),ntttt(nsymz)
      common/bij/bijkl(kget)                                            
      dimension bufs(lsize),bufd(lsize)
c     lecture de la file 40
      read(40)((ij0(i,j),i=1,norb),j=1,norb)
      read(40)nint,(ntttt(i),i=1,nsym)
      write(6,*)' nombre d''integrales lues sur f40', nint
      if(nint.gt.kget)then
	write(6,*)'le nombre d''integrales',nint,' est superieur a la'
     *  ,' taille ',kget,' du tableau bijkl progamme stoppe'
        stop
      end if
      nbuf=nint/lsize+1
      mint=0
      do k=1,nbuf
	if(ydp)then
	  read(40)bufd
	  do i=1,lsize
	    mint=mint+1
	    bijkl(mint)=bufd(i)
	  end do
	else
	  read(40)bufs
	  do i=1,lsize
	    mint=mint+1
	    bijkl(mint)=bufs(i)
	  end do
	end if
      end do
      return                                                            
      end                                                               
