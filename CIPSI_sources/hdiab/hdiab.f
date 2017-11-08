C*****************************************************************
C
C    CONSTRUCTION D UN HAMILTONIEN DIABATIQUE A POSTERIORI
C    PAR RECOUVREMENT MAXIMUM AVEC DES REFEREBCES
C
C       FILES
C       FILE 12   RECOUVREMENTS ENTRE OM ET OM DE REFERENCE(SOM)
c     vecteurs et determinants du calcul adiabatique
c        irvec=4   si cipsi   (ic=0)
c        irvec=58  bd   (ic=1)
c        irdet=60  moy
c     vecteurs et determinants de reference
c        irvec=3   si cipsi   (icr=0)
c        irvec=59  bd (icr=1)
c        irdet=61 moy
C
C
C
C********************************************************************
      program hdiab
      IMPLICIT REAL*8(A-H,O,P,R-Z),LOGICAL*1(Q)
      include 'pshf.prm'
      INTEGER*4 NE,ND,TROU,PART,NER,NDR,PR,TR
      COMMON/CIC/V(ndetz*metz),VR(ndetz*metz),h(metz,metz),
     * NE(ndetz),ND(ndetz),TROU(ntpz),PART(ntpz),
     * NER(ndetz),NDR(ndetz),TR(ntpz),PR(ntpz),
     * num(metz),numr(metz),
     * noca,NROT,NCF,METAT,NREF,NVR,QION,QANAL,QORTOM,QPRT
      dimension hr(metz,metz)
      DIMENSION e(metz),w(ndetz*metz)
      NAMELIST/HEFINP/NOM,NOMR,METAT,NREF,NVR,TR,PR,NER,NDR,
     *    NOCA,QION,QORTOM,QANAL,NROT,VR,num,numr,iprt,e,ic,icr,QPRT

c
c initialisation
c
      ic=1
      icr=1
      iprt=0
      nrot=0
      QORTOM=.FALSE.
      QMOY=.FALSE.
      QANAL=.FALSE.
      QION=.FALSE.
      QPRT=.FALSE.
      NDIM=mz
       do i=1,metz
        e(i)=0.d0
        num(i)=i
        numr(i)=i
       enddo
       write(6,*) 'openf'
      call openf
       write(6,*) ' fin openf'
       write(6,*)
      read(5,hefinp)
 
c
c   choix  energies adiabatiques externes
c
       if(e(1).gt.1.D-10) then
       do i=1,metat
        do j=1,metat
         h(i,j)=zero
        enddo
        h(i,i)=e(i)
       enddo
       endif
       write(6,*) ' energies externes ',(e(i),i=1,metat)
 
c
c   Lecture des vecteurs et determinants de l ic
c
      irvec=58
      irdet=60
c
c   lecture cipsi seulement ic=0 irvec=4 
c
      if(ic.eq.0) irvec=4
      call readic(v,h,ne,nd,trou,part,metat,ncf,irvec,irdet,ic)
      write(6,*) 'energies adiabatiques'
      write(6,*) (h(i,i),i=1,metat)
       do i=1,metat
          e(i)=h(i,i)
       enddo
c
c   Lecture des vecteurs et determinants de reference
c
      irvec=59
      irdet=61
c   lecture cipsi seulement icr=0 irvec=3 
      if(icr.eq.0) irvec=3
      call readic(vr,hr,ner,ndr,tr,pr,nvr,ncfr,irvec,irdet,icr)
      nref=ncfr
      nrefp=101
c      nrefp=ncfr
      do j=1,nvr
        do i=1,nrefp
        w(i+(j-1)*nrefp)=vr(i+(j-1)*nref)
      enddo
      enddo
      do j=1,nvr
        do i=1,nrefp
        vr(i+(j-1)*nrefp)=w(i+(j-1)*nrefp)
      enddo
      enddo
      nref=nrefp
         
c
      if(metat.eq.1) stop
      if(nrot.eq.0) nrot=nvr
c
c    selection des references
c
c
c diabatisation
c
      write(6,*)
      write(6,*) 'vecteurs soumis a rotation'
      write(6,*) ' nrot ',nrot
      write(6,*) ' etats ',(num(i),i=1,nrot)
      write(6,*) ' references ',(numr(i),i=1,nrot)
      write(6,*)
      call hrot
      ndim=metz
      call hvp(h,metat,ndim)
      write(56) ncf,metat,((v(i+(j-1)*ncf),i=1,ncf),j=1,metat),
     *              (e(i),i=1,metat),((h(i,j),i=1,metat),j=1,i)

      stop
      end
