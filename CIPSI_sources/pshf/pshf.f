      program pshf
      implicit real*8 (a-h,o-z)
      include 'pshf.prm'
      common/machin/isingl,nbits                                        
C changed for LINUX compatibility 2/12/97 Malte GROSS
c     character*9 day,hour
      character*26 timeday
C changed for LINUX compatibility 2/12/97 Malte GROSS
      character*8 aname,bname
      logical ymono
      common/output/nprint,itol,icut,normf,normp
      common/restar/timlim,irest,nrec,intloc,ist,jst,kst,lst
      common/iofile/ ir,iw,ip,is,iq,ih,iv,if2,if3,if4,ioda(19)
      common/times/ti,tx,tim,tom,to
      common/nshel/ex(dgp),cs(dgp),cp(dgp),cd(dgp),cf(dgp),cg(dgp),
     1kstart(ds),katom(ds),ktype(ds),kng(ds),kloc(ds),kmin(ds),kmax(ds),
     2nshell
      common/infoa/nat,ich,mul,num,nx,ne,na,nb,zan(dc),c(3,dc),en
      common/symtry/invt(48),iso(ds,12),ptr(3,144),dtr(6,288)
     1,ftr(10,480),nt
      common/scfit/aname,bname,maxit,nconvc,npunch,npi,cextrp,amix,ymono
c     common/pseudo/qn(80)
      common/pseudo/qa(doa),qb(doa),ngla(doa),nglb(doa),numgla,numglb,sc
      common/open/nc,noc(3),oc(3),alpha(3,3),beta(3,3),rland
      ir=5
      iw=6
      ip=7
      is=8
      iq=9
      ih=10
      iv=11
      if2=2
      if3=3
      if4=4
      isingl=2
c     open all files
      call openf
c
C changed 01/2016
      call fdate(timeday)
      write(6,8999) timeday
 8999 format(/,80(1h*),/,8(10h   pshf   ),/, 10x,a25,/,
     * 80(1h*))
      
!     RV 01/16 init times to zero 
      call secnd(to)
      call mole(istop)
      write(6,*) 'irest',irest
      if(irest.gt.2) call daopen
      if(tim.gt.timlim) go to 10000
      if(irest.gt.2) go to 300
c Calcul des integrales bielectroniques 
      if (.not.ymono) call jandk
      if ( istop .eq. 2 ) goto 10000
      if(tim.gt.timlim) go to 10000       
  300 if(irest.gt.3) go to 500
      call standv
      call xpsgr
      call xpsnlc
      call xpsslc  
      call xyzdip
      call chelec  
      if(tim.lt.timlim) call daclos
      if(tim.gt.timlim) go to 10000
  500 if(irest.gt.8) go to 600    
      if(bname.eq.'osrhf') then
         call scfop
      else
         if(numgla.gt.0) then
            call scfgel
         else
            call scf
         endif
      endif
  600 if(tim.gt.timlim) go to 10000
      call hfprop(bname)
10000 continue
      call daclos
      call fdate(timeday)
      write(6,8998) timeday
 8998 format(/,80(1h*),/,8(10h fin pshf ),/, 10x,a25,/,
     * 80(1h*))
C changed for LINUX compatibility 2/12/97 Malte GROSS
 9999 format(/,80(1h*),/,8(10h   pshf   ),/, 10x,a10,5x,a10,/,
     * 80(1h*))
 9998 format(/,80(1h*),/,8(10h fin pshf ),/, 10x,a10,5x,a10,/,
     * 80(1h*))
      stop
      end
