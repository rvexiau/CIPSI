      subroutine otodon(group,ntrans,nshell,c,naxis,norb,ktype,katom,
     * yf)
      implicit real*8 (a-h,o-x,z),logical*1(y)
      include 'pshf.prm'
      character*8 group,grpmin(19),grpmaj(19),dinfh,cinfv,cs,ci
      character*8 cn,dn,cnv,cnh,dnh,c1
      common nrrb,nsym,its(nsymz,nsymz),itsym(doa),itsyv(2*doa)
      common nvec,isymp(nsymz),mijkl,nrec,mblk
      common isydeg(nsymz,nsymz),ntrsy,istop,iarwr(2*doa),
     1 nocsy(nsymz),ymono,yprt,yprtv,yatom,ywvect,ynocsy
      common/om/cvc(doa*doa*2),nbp(nsymz),nbo(nsymz),ind(nsymz+1),
     y ivecsy(2*doa*nsymz)
      common /hondo/ nt,ishell(doa),kloc(doa),ktyp(doa),
     1 numc(doa),newsh(doa,nsymz),iptr(3,nsymz),idtr(6,nsymz),
     2 iftr(10,nsymz),igtr(15,nsymz),yhondo,yprth
      common /ot/ itsdh(14,14),itscv(7,7),its1(1,1),its2(2,2),its4(4,4),
     * its8(8,8),jsym1(1,1),jsym2(2,2),jsym4(4,4),jsym8(8,8),isyded(14,2
     * ),isydec(7,2 )
      data cinfv,dinfh/'CINFV','DINFH'/
      data cn,cs,ci,dn,cnv,cnh,dnh,c1
     1 /'CN','CS','CI','DN','CNV','CNH','DNH','C1'/
      data grpmaj/'C1','CS','CI','CN','S2N','CNH','CNV','DN','DNH',
     1 'DND','CINFV','DINFH','T','TH','TD','O','OH','I','IH'/
      data grpmin/'c1','cs','ci','cn','s2n','cnh','cnv','dn','dnh',
     1 'dnd','cinfv','dinfh','t','th','td','o','oh','i','ih'/
      dimension c(3,30),jsym(nsymz,nsymz)
      dimension idinf1(20),idinf2(28)
      dimension ktype(doa),katom(doa)
      character*6 ante
      common /info/ante,yex
c
      data idinf1/ 0,2,1,0,-1,0,3,0,-1,1,4,0,-1,-1,0,0,-1,-1,2,0/
     *     idinf2/ 0,2,0,0,-1,0,1,1,-1,1,2,2,-1,-1,3,3,-1,-1,4,4,
     *             -1,-1,-1,5,-1,-1,-1,6/
      nsym=nt
      do i=1,19
      if(group.eq.grpmin(i))then
      group=grpmaj(i)
      write(6,*)'groupe',group
      end if
      end do
c
      if(group.ne.cinfv.and.group.ne.dinfh) goto 4
      yp=.false.
      yd=.false.
      yf=.false.
      do 1 i=1,nshell
      if(ktype(i).eq.2) yp=.true.
      if(ktype(i).eq.3) yd=.true.
      if(ktype(i).eq.4) yf=.true.
1     continue
      if(yf) then
      yp=.true.
      yd=.true.
      end if
      if(yd) then
      yp=.true.
      end if
      if(yp) ntrans=ntrans+1
      if(yd) ntrans=ntrans+1
4     continue
      if(group.ne.dinfh) goto 5
c
c***********************************groupe dinfh
      nsym=14
      if(.not.yf) nsym=10
c     if(.not.yd) nsym=6
      yg=.false.
      ish=-norb
c
      do 10 ink=1,nsym
      yg=.not.yg
      ish=ish+norb
      jnk=(ink-1)/2
      i2=0
      do 10 i=1,nshell
      idp=idinf2 (ktype(i)+jnk*4)
      if(idp.lt.0) go to 9
      i1=ish+kloc(i)+idp
      yki=(ktype(i)-2*(ktype(i)/2)).eq.0
      ysym=yki.neqv.yg
       if(.not.ysym.and.dabs(c(3,katom(i))).lt.1.d-8) go to 9
      if(ivecsy(i1).ne.0) go to 9
      i2=i2+1
      ivecsy(i1)=i2
      if(dabs(c(3,katom(i))).gt.1.d-8) then
      i1=ish+kloc(newsh(i,2))+idp
      if(ivecsy(i1).eq.0) then
      ivecsy(i1)=i2
      if(.not.ysym) ivecsy(i1)=-i2
      end if
      end if
    9 continue
c     write(6,2345) ink,jnk,i,ktype(i),idp,kloc(i),i1,yg,yki,ysym,
c    * (ivecsy(k),k=(ink-1)*norb+1,ink*norb)
 2345 format(7i3,3l4,(/,20i3))
   10 continue
      do 681 j=1,nsym
      do 69 i=1,nsym
69    its(i,j)=itsdh(i,j)
      do 70 i=1,ntrans
70    isydeg(j,i)=isyded(j,i)
681   continue
      if(yp.or.ntrans.ne.1) goto 682
      do 683 j=1,nsym
683   isydeg(j,1)=isyded(j,2)
682   continue
      return
5     continue
 
c****************************************groupe cinfv
      if(group.ne.cinfv) goto 6
      nsym=7
      if(.not.yf) nsym=5
c     if(.not.yd) nsym=3
      ish=-norb
      do 20 ink=1,nsym
      jnk=ink-1
      i2=0
      ish=ish+norb
      do 20 i=1,nshell
      idp=idinf2 (ktype(i)+jnk*4)
      if(idp.lt.0) go to 20
      i1=ish+kloc(i)+idp
      if(ivecsy(i1).ne.0) go to 19
      i2=i2+1
      ivecsy(i1)=i2
   19 continue
c      write(6,2345) ink,jnk,i,idp,kloc(i),i1,
c     * (ivecsy(k),k=(ink-1)*norb+1,ink*norb)
  20  continue
      do 71 j=1,nsym
      do 72 i=1,nsym
72    its(i,j)=itscv(i,j)
      do 73 i=1,ntrans
73    isydeg(j,i)=isydec(j,i)
71    continue
c
      if(yp.or.ntrans.ne.1) goto 711
      do 712 j=1,nsym
712   isydeg(j,1)=isydec(j,2)
711   continue
c
      return
c
c******************************tous groupes non degeneres
6     continue
      if(group.eq.c1) go to 94
      if(naxis.ne.2.and.group.ne.cs.and.group.ne.ci) goto 90
      goto 91
90    write(6,1000)
1000                format('  programme stoppe : groupe degenere')
      stop
91    if(group.eq.cn) goto 98
      if(group.eq.cs) goto 98
      if(group.ne.ci) goto 92
98    nsym=2
      do 981 j=1,nsym
      do 981 i=1,nsym
      its(i,j)=its2(i,j)
981   jsym(i,j)=jsym2(i,j)
      goto 95
92    if(group.eq.dn) goto 97
      if(group.eq.cnv)goto 97
      if(group.ne.cnh)goto 93
97    nsym=4
c     write(6,*) '97',its4,jsym4
      do 971 j=1,nsym
      do 971 i=1,nsym
      its(i,j)=its4(i,j)
971   jsym(i,j)=jsym4(i,j)
      goto 95
93    if(group.ne.dnh) goto 94
      nsym=8
      do 931 j=1,nsym
      do 931 i=1,nsym
      its(i,j)=its8(i,j)
931   jsym(i,j)=jsym8(i,j)
      goto 95
94    if(group.ne.c1) goto 96
      nsym=1
      its(1,1)=its1(1,1)
      jsym(1,1)=jsym1(1,1)
      do 941 i=1,norb
941   ivecsy(i)=i
      return
96    write(6,1001)
1001  format(' groupe non trouve, revenir aux donnees non automatiques')
      stop
95    continue
c     do 1111 i=1,nsym
c1111 write(6,*)'jsym',(jsym(i,j),j=1,nsym)
      i1=0
      do 100 i=1,nsym
      im1n=(i-1)*norb
      i2=0
      do 110 j=1,nshell
      if(ante.eq.'  pshf')iii=2*(ktype(j)-1)+1
      if(ante.eq.'gamess')iii=(ktype(j)*(ktype(j)+1))/2
      if(yex)iii=(ktype(j)*(ktype(j)+1))/2
      do 120 k=1,iii
      i1=i1+1
      if(ivecsy(i1).ne.0) goto 120
      i2=i2+1
      do 130 l=1,nt
      k1=kloc(newsh(j,l))+k-1+im1n
      k2=jsym(i,l)
      if(ktype(j).eq.2) k2=k2*iptr(k,l)
      if(ktype(j).eq.3) k2=k2*idtr(k,l)
      if(ktype(j).eq.4) k2=k2*iftr(k,l)
      i3=i2*k2
1112  format(1x,30i4)
      if(ivecsy(k1).ne.0.and.ivecsy(k1).ne.i3) goto 140
      ivecsy(k1)=i3
      if(yprtv) write(6,1112) i,j,k,l,i1,i2,k1,k2,i3
      goto 130
140   do 150 ll=1,nt
      k4=kloc(newsh(j,ll))+k-1+im1n
150   ivecsy(k4)=0
      i2=i2-1
      goto 120
130   continue
120   continue
110   continue
100   continue
      return
      end
