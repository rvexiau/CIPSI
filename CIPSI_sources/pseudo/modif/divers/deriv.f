C       PROGRAM DERIVEE                                                         
C                                                                               
C                                                                               
C       ****************************************************************        
C	****************************************************************              
C                                                                               
C                                                                               
C	CE PROGRAMME GENERE PAR DERIVATION LES EXPRESSIONS FORMELLES DES              
C	INTEGRALES : 1/R**N ET X,Y,Z/R**N POUR LE CALCUL DE L'INTERAC-                
C	TION COEUR-VALENCE.                                                           
C                                                                               
C	CF. REFERENCES :                                                              
C       --------------                                                          
C	- J.C. BARTHELAT, THESE DE DOCTORAT D'ETAT, UNIVERSITE PAUL-                  
C	  SABATIER, TOULOUSE , FRANCE (1977).                                         
C                                                                               
C         ( EGALEMENT THESE DE A. SERAFINI ).                                   
C                                                                               
C                                                                               
C	................................................................              
C	................................................................              
C                                                                               
C                                                                               
C   **	NOUS POUVONS OBTENIR AINSI LES EXPRESSIONS DES INTEGRALES DE LA          
C	FORME : <S/PL/S>, <P/PL/S>, <P/PL/P>, <D/PL/S>, ... <F/PL/S>,                 
C	... (PAR COUCHE).                                                             
C		AVEC <S/, <P/ ... CENTRES EN A                                               
C		  ET /S>, /P> ... CENTRES EN B                                               
C		LE PSEUDOPOTENTIEL ATOMIQUE EST CENTRE EN C.                                 
C   **	LES FORMULES DES DERIVEES POUR LES PL (PROJECTEUR SUR LE L-IEME          
C	SOUS-ESPACE DES HARMONIQUES SPERIQUES CENTRES EN C) SONT PRO-                 
C	GRAMMES POUR P1, P2 ET P3.                                                    
C                                                                               
C                                                                               
C  ***  NOTONS QU'IL Y A UNE NUANCE ENTRE LES EXPRESSIONS DES INTEGRALES        
C       DE 1/R**N ET X,Y,Z/R**3 AU NIVEAU DE LA DEFINITION DES PL. (CF.         
C       DEMONSTRATION M.F.)                                                     
C       DE PLUS LE CAS 1/R**4 EST TRAITE POUR UN N GENERAL, C'EST-A-DIRE        
C       QUE L'ON A POSE N=0 TANDIS QUE POUR X,Y,Z/R**3, N=2 :                   
C                    <A/(Z/RC**3)/B>=<A/(COS(TETA)/RC**2)/B>                    
C                              Z : VECTEUR                                      
C                                                                               
C                                                                               
C                                                                               
C	PARAMETRES UTILISES DANS LES EXPRESSIONS ANALYTIQUES DES DERI-                
C	VEES ( CF. REFERENCES ENUMEREES CI-HAUT ) :                                   
C                                                                               
C	................................................................              
C	COEFF 2A 2B A+B DELTAS I(K,L,N) PL A AI...K AOP B BI...K BOP                  
C	................................................................              
C                                                                               
C	*** NOTONS QUE 2A ET 2B SYMBOLISENT 2*ALPHA ET 2*BETA RESPECTI-               
C	    VEMENT. A+B : ALPHA+BETA.                                                 
C	    DELTAS POUR LE DELTA DE KRONECKER.                                        
C	    AOP ET BOP SONT UTILISES LORS DES DERIVATIONS PAR RAPPORT                 
C 	    A : "C" (POUR X,Y,Z/R**3).                                               
C           ( OBTENTION DES INTEGRALES DU CHAMP ELECTRIQUE ).                   
C                                                                               
C                                                                               
C	----------------------------------------------------------------              
C	----------------------------------------------------------------              
C                                                                               
C                                                                               
C 	LES SUBROUTINES :                                                            
C	***************                                                               
C                                                                               
C       N.B. : LES SUBROUTINES DU PROGRAMME DERIVEE SONT CLASSEES PAR           
C       ****   ORDRE ALPHABETIQUE.                                              
C                                                                               
C	ADD2(IJ) : SERT A ADDITIONNER LES INTEGRALES SUPPLEMENTAIRES                  
C	           A CELLE QUE L'ON DERIVE POUR OBTENIR L'INTEGRALE                   
C	           (IJ) DESIREE.                                                      
C                                                                               
C	CASPART : TRAITE LES CAS PARTICULIERS, C'EST-A-DIRE LES CAS AVEC              
C                 AC=0, BC=0 AINSI QUE AC ET BC=0. NOUS ECRIVONS DONC :         
C                 (1) A=0, (2) B=0 ET (3) A ET B=0.                             
C                                                                               
C	CONTRAC(IJ) : CONTRACTE LES TERMES FINALS POUR UNE INTEGRALE                  
C	              (IJ) DONNEE, C'EST-A-DIRE QU'ELLE ADDITIONNE LES                
C	              TERMES QUI NE DIFFERENT QUE PAR LE COEFFICIENT.                 
C                               ( COEFFICIENT : COEFF )                         
C                                                                               
C	DERIV(ITYPE) : GENERE LES EXPRESSIONS DES DERIVEES POUR CHACUN                
C	               DES TERMES DERIVABLES DU TABLEAU TAB(NPARTAB).                 
C	               LE TYPE DE DERIVATION EST TRANSMIS PAR ARGUMENT.               
C                                                                               
C	DERIVC(IJ) : GENERE LES TERMES DES DERIVEES DE L'INTEGRALE                    
C                    (IJ=1) PAR RAPPORT A | C.                                  
C                                                                               
C       ECRIT : CETTE SUBROUTINE SERT A ECRIRE LES RESULTATS DES DERI-          
C               VATIONS SUR LA FILE 10 (FILE SEQUENTIELLE).                     
C                                                                               
C	INTEGEX : GENERE LES EXPRESSIONS DES DIFFERENTS TYPES D'INTE-                 
C	          GRALES A CALCULER.                                                  
C	                                                                              
C                                                                               
C	----------------------------------------------------------------              
C	----------------------------------------------------------------              
C                                                                               
C                                                                               
C       LES RESULTATS SONTS ECRITS SUR LA FILE 10 : " DERPL03 "                 
C       (CETTE FILE EST UNE FILE SEQUENTIELLE)                                  
C                                                                               
C                                                                               
C       ----------------------------------------------------------------        
C       ----------------------------------------------------------------        
C                                                                               
C                                                                               
C	DEFINITIONS : ( CF. CHACUNE DES SUBROUTINES)                                  
C	...........                                                                   
C                                                                               
C 	NPARTAB : LE NOMBRE TOTAL DE PARAMETRES, C'EST-A-DIRE :                      
C       *******                                                                 
C                                                                               
C	................................................................              
C       COEFF 2A 2B A+B DELTAS I(K,L,N) PL A AI...K AOP B BI...K BOP            
C	................................................................              
C                                                                               
C                                                                               
C 	NDELTA : NOMBRE DE DELTAS.                                                   
C       ******                                                                  
C                                                                               
C	IJMAX : NOMBRE MAXIMAL D'INTEGRALES QUI SERONT TRAITEES.                      
C       *****                                                                   
C                                                                               
C                                                                               
C	----------------------------------------------------------------              
C	----------------------------------------------------------------              
C                                                                               
C                                                                               
C                                                                               
C	****************************************************************              
C	****************************************************************              
C                                                                               
C                                                                               
C                                                                               
        IMPLICIT INTEGER (A-Z)                                                  
        DOUBLE PRECISION COEFF                                                  
        PARAMETER (ITYPMAX=3,INTMAX=2,NTBFIN=45000)                             
	COMMON/PARGEN/NPARTAB,NDELTA,IJMAX                                             
        COMMON/TABPAR/TABA,TABB,TABPL,TABK,TABL,TABN,TABAOP,TABBOP              
        COMMON/INT/PARIJ((ITYPMAX+1)*(ITYPMAX+2)/2,2+ITYPMAX)                   
        COMMON/DERADD2/NDEB(INTMAX,0:3,(ITYPMAX+1)*(ITYPMAX+2)/2),              
     1                 NFIN(INTMAX,0:3,(ITYPMAX+1)*(ITYPMAX+2)/2)               
        COMMON/DERI/TAB(15+ITYPMAX*(2*ITYPMAX+5)),                              
     1              TABAUX(NTBFIN,15+ITYPMAX*(2*ITYPMAX+5))                     
        COMMON/DERIVE/NTF                                                       
	COMMON/AUX/NTERMF                                                              
        COMMON/LECRIT/ENT,SOR                                                   
        COMMON/ENREG/ENR                                                        
        COMMON/INTEG/ITYPINT,ICAMAX                                             
	INTEGER*2 TAB,TABAUX,TABC(15+ITYPMAX*(2*ITYPMAX+5))                            
        CHARACTER*20 FILE9,FILE10                                               
        EQUIVALENCE (TAB(1),COEFF)                                              
	CHARACTER*5 CHAR(10)                                                           
        DATA CHAR/'<S/S>','<P/S>','<P/P>','<D/S>','<D/P>','<D/D>',              
     1            '<F/S>','<F/P>','<F/D>','<F/F>'/                              
C                                                                               
C                                                                               
C                                                                               
 4900   FORMAT(////,8X,'TABLEAU CONTENANT LES TERMES DES DERIVEES',             
     1         /,8X,41('*'),//)                                                 
 4901   FORMAT(///,1X,'INTEGRALE NO : ',I5,5X,A5,/)                             
 4902   FORMAT(1X,'NOMBRE DE TERMES FINALS BRUTS (NTERMF) : ',I5,/,1X,          
     1            'NOMBRE DE TERMES FINALS CONTRACTES (NTF) :',I5,/,1X,         
     1            'NDEB = ',I10,10X,'NFIN = ',I10,//)                           
 4904   FORMAT(1X,'NDEB = ',I10,10X,'NFIN = ',I10,//)                           
C4910	FORMAT(/,1x,'2ALPHA 2BETA ALPHA+BETA DELTAS I(K,L,N) PL A'                
C    1         ' AI...K AOP B BI...BK BOP COEFF',//)                            
C                                                                               
C                                                                               
C                                                                               
        ENT=5                                                                   
        SOR=6                                                                   
C                                                                               
C                                                                               
C                                                                               
	NPARTAB=15+ITYPMAX*(2*ITYPMAX+5)                                               
        NDELTA=ITYPMAX*(2*ITYPMAX+1)                                            
  	IJMAX=(ITYPMAX+1)*(ITYPMAX+2)/2                                              
      	NRECL=8*NPARTAB                              
CVAX
CVAX	NRECL=2*NPARTAB                                                                
CVAX
C                                                                               
C                                                                               
	TABA=12+NDELTA                                                                 
	TABB=14+NDELTA+2*ITYPMAX                                                       
	TABPL=11+NDELTA                                                                
	TABK=8+NDELTA                                                                  
	TABL=9+NDELTA                                                                  
	TABN=10+NDELTA                                                                 
	TABAOP=TABA+2*ITYPMAX+1                                                        
	TABBOP=TABB+2*ITYPMAX+1                                                        
C                                                                               
C                                                                               
C                                                                               
C       LA FILE 9 EST UNE FILE DE TRAVAIL A ACCES DIRECT.                       
C                                                                               
        OPEN(9,STATUS='SCRATCH',ACCESS='DIRECT',                   
     1                            RECL=NRECL)                                   
        OPEN(10,FILE='DERPL03',FORM='UNFORMATTED',STATUS='UNKNOWN')             
        WRITE(SOR,*) 'ON ECRIT LES RESULTATS SUR DERPL22'                       
C                                                                               
        REWIND 10                                                               
C                                                                               
C                                                                               
C                                                                               
        CALL INTEGEX                                                            
C                                                                               
C                                                                               
C                                                                               
        DO 5000 PL=0,3                                                          
        WRITE(SOR,*) 'PL=',PL                                                   
        WRITE(SOR,4900)                                                         
C                                                                               
C                                                                               
        DO 5 I=5,NPARTAB                                                        
        TAB(I)=0                                                                
 5      CONTINUE                                                                
C                                                                               
C                                                                               
        COEFF=1.0D0                                                             
        TAB(TABPL)=PL                                                           
        TAB(TABK)=PL                                                            
        TAB(TABL)=PL                                                            
C                                                                               
C                                                                               
C                                                                               
C       ITYPINT=1 POUR 1/R**N ET ITYPINT=2 POUR X,Y,Z/R**3                      
C       NTERMF=NOMBRE DE TERMES FINALS BRUTS.                                   
C       NTF=NOMBRE DE TERMES FINALS CONTRACTES.                                 
C       IJ : NUMERO DE L'INTEGRALE A CALCULER                                   
C                 EX : IJ=1 <S/PL/S>                                            
C                      IJ=2 <P/PL/S>                                            
C                      ...                                                      
C                      IJ=10 <F/PL/F>                                           
C                                                                               
C                                                                               
C                                                                               
        ITYPINT=1                                                               
        WRITE(SOR,*) 'ITYPINT=',ITYPINT                                         
C                                                                               
C                                                                               
C                                                                               
	IJ=1                                                                           
C                                                                               
C                                                                               
C	SI ITYPINT=1, IJ=1 : <S/S>.                                                   
C                                                                               
                NDEB(ITYPINT,0,IJ)=1                                            
                NFIN(ITYPINT,0,IJ)=1                                            
                WRITE(SOR,4901) IJ,CHAR(IJ)                                     
                WRITE(SOR,4904) NDEB(ITYPINT,0,IJ),NFIN(ITYPINT,0,IJ)           
C		WRITE(SOR,4910)                                                              
                ENR=1                                                           
                WRITE(9,REC=ENR) TAB                                            
C                                                                               
C                                                                               
C                                                                               
C       POUR LA DEFINITION DU TABLEAU PARIJ, CF. SUBROUTINE INTEGEX.            
C                                                                               
C                                                                               
C                                                                               
 200    DO 1000 IJ=2,IJMAX                                                      
C                                                                               
        NTERMF=0                                                                
        KL=PARIJ(IJ,1)                                                          
        AB=PARIJ(IJ,2)                                                          
        DO 500 K=NDEB(ITYPINT,0,KL),NFIN(ITYPINT,0,KL)                          
        READ(9,REC=K) TAB                                                       
        CALL DERIV(AB)                                                          
 500    CONTINUE                                                                
        CALL ADD2(IJ)                                                           
        CALL CONTRAC(IJ)                                                        
        NDEB(ITYPINT,0,IJ)=NFIN(ITYPINT,0,IJ-1)+1                               
        NFIN(ITYPINT,0,IJ)=NFIN(ITYPINT,0,IJ-1)+NTF                             
C                                                                               
        WRITE(SOR,4901) IJ,CHAR(IJ)                                             
        WRITE(SOR,4902) NTERMF,NTF,NDEB(ITYPINT,0,IJ),NFIN(ITYPINT,0,IJ)        
C	WRITE(SOR,4910)                                                               
C                                                                               
 1000   CONTINUE                                                                
C                                                                               
C                                                                               
C	TRAITEMENT DES CAS PARTICULIERS :                                             
C                                                                               
        ICAMAX=3                                                                
        CALL CASPART                                                            
C                                                                               
C                                                                               
C                                                                               
        IF(ITYPINT.EQ.2) GO TO 9999                                             
        ITYPINT=2                                                               
        WRITE(SOR,*) 'ITYPINT=',ITYPINT                                         
        CALL DERIVC(PL)                                                         
        GO TO 200                                                               
C                                                                               
C                                                                               
 9999   CALL ECRIT                                                              
C                                                                               
C                                                                               
 5000   CONTINUE                                                                
C                                                                               
C                                                                               
        STOP                                                                    
        END                                                                     
C                                                                               
C                                                                               
C                                                                               
C       ****************************************************************        
C       ****************************************************************        
C                                                                               
C                                                                               
C                                                                               
        SUBROUTINE ADD2(IJ)                                                     
C                                                                               
C       ................................................................        
C                                                                               
C	CETTE SUBROUTINE SERT A ADDITIONNER LES INTEGRALES SUPPLEMEN-                 
C	TAIRES A CELLE QUE L'ON DERIVE POUR OBTENIR  L' INTEGRALE 'IJ'                
C	DESIREE. LE NUMERO DE L'INTEGRALE A ADDITIONNER EST CONTENU                   
C	DANS PARIJ(IJ,3). DE PLUS, ON INVERSE LE ROLE DE CERTAINS IN-                 
C	DICES SI L'INTEGRALE A ADDITIONNER N'EST PAS CONNUE.                          
C	<IJK/LMN>= 1/2*B [ DEL/DEL BN <IJK/LM> + DEL LN <IJK/M> +                     
C		                                DEL MN <IJK/L> ]                             
C		INVERSER LE ROLE DE L ET DE M.                                               
C                                                                               
C       ................................................................        
C                                                                               
C                                                                               
C                                                                               
        IMPLICIT INTEGER (A-Z)                                                  
	DOUBLE PRECISION COEFF                                                         
        PARAMETER (ITYPMAX=3,INTMAX=2,NTBFIN=45000)                             
	COMMON/PARGEN/NPARTAB,NDELTA,IJMAX                                             
        COMMON/DERI/TAB(15+ITYPMAX*(2*ITYPMAX+5)),                              
     1              TABAUX(NTBFIN,15+ITYPMAX*(2*ITYPMAX+5))                     
        COMMON/DERIVE/NTF                                                       
        COMMON/DERADD2/NDEB(INTMAX,0:3,(ITYPMAX+1)*(ITYPMAX+2)/2),              
     1                 NFIN(INTMAX,0:3,(ITYPMAX+1)*(ITYPMAX+2)/2)               
        COMMON/INT/PARIJ((ITYPMAX+1)*(ITYPMAX+2)/2,2+ITYPMAX)                   
	COMMON/AUX/NTERMF                                                              
        COMMON/INTEG/ITYPINT                                                    
        INTEGER*2 TAB,TABAUX                                                    
	EQUIVALENCE (TAB(1),COEFF)                                                     
C                                                                               
C                                                                               
C                                                                               
C       POUR LA DEFINITION DU TABLEAU PARIJ, CF. SUBROUTINE INTEGEX,            
C                                                                               
C                                                                               
        IF(PARIJ(IJ,2).LT.4)THEN                                                
                L=0                                                             
                            ELSE                                                
                L=1                                                             
        ENDIF                                                                   
C                                                                               
C                                                                               
	IF(PARIJ(IJ,3).EQ.0) GO TO 3000                                                
C                                                                               
C                                                                               
	ITYPE=PARIJ(IJ,2)                                                              
	M=PARIJ(IJ,3)                                                                  
	TABAIP=12+NDELTA                                                               
	TABBIP=14+NDELTA+2*ITYPMAX                                                     
C                                                                               
C	N=NOMBRE DE FOIS QUE LA BOUCLE DO 2500 SERA EXECUTEE, C'EST-A-                
C	DIRE LE NOMBRE D'INTEGRALE(S) A ADDITIONNER.                                  
C       SI L'ON DERIVE A GAUCHE N=ITYPE-1 ET SI L'ON DERIVE A DROITE,           
C	N=ITYPE-ITYPMAX-1.                                                            
C                                                                               
	IF(ITYPE.LE.ITYPMAX)THEN                                                       
		N=ITYPE-1                                                                     
		            ELSE                                                              
		N=ITYPE-ITYPMAX-1                                                             
	ENDIF                                                                          
	NADD=0                                                                         
C                                                                               
C                                                                               
	DO 2500 I=1,N                                                                  
        NADD=NADD+1                                                             
	IF(ITYPE.LE.ITYPMAX)THEN                                                       
		II=I                                                                          
                            ELSE                                                
		II=ITYPMAX+I                                                                  
	ENDIF                                                                          
C                                                                               
	DO 2000 J=NDEB(ITYPINT,0,M),NFIN(ITYPINT,0,M)                                  
        READ(9,REC=J) TAB                                                       
C                                                                               
C	LES DELTAS SONT CONTENUS DANS PARIJ(IJ,3+NADD) :                              
C                                                                               
	DELTA=PARIJ(IJ,3+NADD)                                                         
C                                                                               
C	IL N'Y A PAS D'INVERSION POUR LA DERNIERE INTEGRALE :                         
C                                                                               
	IF(NADD.EQ.N) GO TO 1000                                                       
C                                                                               
C                                                                               
C	                                                                              
C	INVERSION DES ROLES DE AI ET AK-1 :                                           
C                                                                               
C       - INVERSION DES INDICES DES A ET DES B :                                
C                                                                               
	IAUX=TAB(TABAIP+ITYPE-1)                                                       
	TAB(TABAIP+ITYPE-1)=TAB(TABAIP+II)                                             
	TAB(TABAIP+II)=IAUX                                                            
	IAUX=TAB(TABBIP+ITYPE-1)                                                       
	TAB(TABBIP+ITYPE-1)=TAB(TABBIP+II)                                             
        TAB(TABBIP+II)=IAUX                                                     
C                                                                               
C	-INVERSION DES INDICES DES DELTAS :                                           
C                                                                               
	DO 300 JJ=1,2*ITYPMAX+1                                                        
C                                                                               
C       ON SAUTE LE CAS JJ=II CAR PAS DE DELTA DIAGONAL.                        
C       ON SAUTE JJ=ITYPE-1 CAR ON VEUT INVERSER LE ROLE DE II ET DE            
C       ITYPE-1                                                                 
C                                                                               
	IF(JJ.EQ.II .OR. JJ.EQ.ITYPE-1) GO TO 300                                      
	IF(JJ.GT.II) GO TO 250                                                         
	JI=(II-1)*(II-2)/2+JJ                                                          
	KI=(ITYPE-2)*(ITYPE-3)/2+JJ                                                    
	IAUX=TAB(7+JI)                                                                 
        TAB(7+JI)=TAB(7+KI)                                                     
        TAB(7+KI)=IAUX                                                          
	GO TO 300                                                                      
 250	IF(JJ.GT.ITYPE-1) GO TO 275                                                
	JI=(JJ-1)*(JJ-2)/2+II                                                          
	KI=(ITYPE-2)*(ITYPE-3)/2+JJ                                                    
	IAUX=TAB(7+JI)                                                                 
	TAB(7+JI)=TAB(7+KI)                                                            
	TAB(7+KI)=IAUX                                                                 
	GO TO 300                                                                      
 275	JI=(JJ-1)*(JJ-2)/2+II                                                      
	KI=(JJ-1)*(JJ-2)/2+ITYPE-1                                                     
	IAUX=TAB(7+JI)                                                                 
	TAB(7+JI)=TAB(7+KI)                                                            
	TAB(7+KI)=IAUX                                                                 
 300	CONTINUE                                                                   
C                                                                               
C                                                                               
C	DIVISION PAR 2*ALPHA OU PAR 2*BETA SELON LE TYPE DE DERIVATION.               
C  	ADDITION DES TERMES SUPPLEMENTAIRES. CES TERMES SONT MIS DANS               
C	LE TABLEAU TABAUX(NTERMF,NPARTAB).                                            
C                                                                               
 1000   TAB(DELTA)=1                                                            
	TAB(5+L)=TAB(5+L)-1                                                            
	NTERMF=NTERMF+1                                                                
	IF(NTERMF.GT.NTBFIN) GO TO 3100                                                
	DO 1500 K=1,NPARTAB                                                            
	TABAUX(NTERMF,K)=TAB(K)                                                        
 1500	CONTINUE                                                                  
C                                                                               
 2000   CONTINUE                                                                
C                                                                               
C                                                                               
 2500	CONTINUE                                                                  
C                                                                               
C                                                                               
 3000   RETURN                                                                  
 3100   PRINT*,'NTBFIN INSUFFISANT, NTERMF = ',NTERMF                           
C                                                                               
C                                                                               
        STOP                                                                    
        END                                                                     
C                                                                               
C                                                                               
C                                                                               
C       ****************************************************************        
C       ****************************************************************        
C                                                                               
C                                                                               
C                                                                               
        SUBROUTINE CASPART                                                      
C                                                                               
C       ----------------------------------------------------------------        
C                                                                               
C       CETTE SUBROUTINE TRAITE LES CAS PARTICULIERS, C'EST-A-DIRE :            
C               ICAS( ) : ( 1 ) A=0                                             
C               .......   ( 2 ) B=0                                             
C                         ( 3 ) A ET B=0                                        
C                                                                               
C       *** WARNING : IL FAUDRAIT PEUT-ETRE APPELE LA SUBROUTINE :              
C           *******        " CONTRAC "...                                       
C                                                                               
C       ----------------------------------------------------------------        
C                                                                               
C	DEFINITIONS :                                                                 
C	***********                                                                   
C                                                                               
C	NTFCP : NOMBRE DE TERMES FINALS POUR UN CAS PARTICULIER.                      
C                                                                               
C	TABCP(NTFCP,NPARTAB) : TABLEAU CONTENANT LES TERMES FINALS DE                 
C	                       CHAQUE CAS PARTICULIER.                                
C                                                                               
C       ---------------------------------------------------------------         
C                                                                               
        IMPLICIT INTEGER (A-Z)                                                  
        PARAMETER (ITYPMAX=3,INTMAX=2,NTBFIN=45000,DENU=15)                     
        DOUBLE PRECISION COEFF,TABDENU(0:DENU)                                  
	COMMON/PARGEN/NPARTAB,NDELTA,IJMAX                                             
        COMMON/TABPAR/TABA,TABB,TABPL,TABK,TABL,TABN,TABAOP,TABBOP              
        COMMON/DERI/TAB(15+ITYPMAX*(2*ITYPMAX+5)),                              
     1              TABAUX(NTBFIN,15+ITYPMAX*(2*ITYPMAX+5))                     
        COMMON/DERIVE/NTF                                                       
        COMMON/DERADD2/NDEB(INTMAX,0:3,(ITYPMAX+1)*(ITYPMAX+2)/2),              
     1                 NFIN(INTMAX,0:3,(ITYPMAX+1)*(ITYPMAX+2)/2)               
        COMMON/LECRIT/ENT,SOR                                                   
        COMMON/ENREG/ENR                                                        
        COMMON/INTEG/ITYPINT                                                    
        INTEGER*2 TAB,TABAUX                                                    
	DIMENSION PKL(2),PAB(2)                                                        
        EQUIVALENCE (TAB(1),COEFF)                                              
	CHARACTER*5 CHAR(10)                                                           
        DATA CHAR/'<S/S>','<P/S>','<P/P>','<D/S>','<D/P>','<D/D>',              
     1            '<F/S>','<F/P>','<F/D>','<F/F>'/                              
C                                                                               
C                                                                               
C                                                                               
  700   FORMAT(/////,1X,'CAS PARTICULIERS , 1 : (A=0), 2 : (B=0) ET'            
     1                ' 3 : (A ET B=0)',/,1X,57('*'))                           
  800   FORMAT(////,1X,'CAS NO : ',I5,/,1X,8('*'),//)                           
  850   FORMAT(1X,'INTEGRALE NO : ',I5,'CAS : ',I5,'P<0 : FAUX')                
  851   FORMAT(1X,'NTFCP > NTBFIN, NTFCP =',I5)                                 
  900   FORMAT(8X,'TABLEAU CONTENANT LES TERMES DES DERIVEES',                  
     1         ' : TABCP(NTFIN,24)',/,8X,59('*'),//)                            
  901   FORMAT(///,1X,'INTEGRALE NO : ',I5,5X,A5,/)                             
  902   FORMAT(1X,'NOMBRE DE TERMES FINALS (NTFCP) : ',I5,/,1X,                 
     1            'NDEB1 = ',I10,10X,'NFIN1 = ',I10,//)                         
C 903	FORMAT(/,1x,'2ALPHA 2BETA ALPHA+BETA DELTAS I(K,L,N) PL A'                
C    1         ' AI...K AOP B BI...BK BOP COEFF',//)                            
C                                                                               
C                                                                               
        WRITE(SOR,700)                                                          
C                                                                               
C                                                                               
C	ON GENERE LES VALEURS DES 1*3...(2K+1)                                        
C                                                                               
        TABDENU(0)=1                                                            
        DO 5 I=1,15                                                             
        L=2*I+1                                                                 
        TABDENU(I)=L*TABDENU(I-1)                                               
 5      CONTINUE                                                                
C                                                                               
C                                                                               
C	PKL : ADRESSES DES K ET L DE I(KLN).                                          
C	PAB : ADRESSES DES A ET DES B.                                                
C                                                                               
	PKL(1)=8+NDELTA                                                                
	PKL(2)=9+NDELTA                                                                
	PAB(1)=12+NDELTA                                                               
	PAB(2)=14+NDELTA+2*ITYPMAX                                                     
C                                                                               
C                                                                               
C  	ON TRAITE LES TROIS CAS ( (1): A=0, (2): B=0 ET (3) A ET B=0 )              
C	L'UN A LA SUITE DE L'AUTRE. (CF. BOUCLE DO 500).                              
C                                                                               
C                                                                               
        DO 500 ICAS=1,3                                                         
        WRITE(SOR,800) ICAS                                                     
C                                                                               
C                                                                               
        DO 400 IJ=1,IJMAX                                                       
        ENR1=ENR                                                                
        NTFCP=0                                                                 
C                                                                               
C                                                                               
        DO 300 I=NDEB(ITYPINT,0,IJ),NFIN(ITYPINT,0,IJ)                          
C                                                                               
        READ(9,REC=I) TAB                                                       
C                                                                               
C                                                                               
        IF(ICAS.EQ.3) GO TO 100                                                 
C                                                                               
C                                                                               
C	LES 2 PREMIERS CAS : (1) A=0, (2) B=0 :                                       
C       .....................................                                   
C                                                                               
	NOKL=PKL(ICAS)                                                                 
	VALKL=TAB(NOKL)                                                                
	SOM=VALKL                                                                      
	DO 15 II=0,2*ITYPMAX+1                                                         
	SOM=SOM+TAB(PAB(ICAS)+II)                                                      
	TAB(PAB(ICAS)+II)=0                                                            
 15	CONTINUE                                                                    
C                                                                               
        IF(SOM) 30,35,300                                                       
 30             WRITE(SOR,850)IJ,ICAS                                           
                STOP                                                            
 35             NTFCP=NTFCP+1                                                   
		COEFF=COEFF/TABDENU(VALKL)                                                    
		TAB(10+NDELTA)=TAB(10+NDELTA)-TAB(NOKL)                                       
                IF(TAB(TABPL).NE.0)THEN                                         
			WRITE(SOR,*) 'ERREUR DANS CASPART'                                           
                        WRITE(SOR,*) 'ICAS=',ICAS,'IJ=',IJ                      
			WRITE(SOR,*) 'TAB(TABPL)=',TAB(TABPL)                                        
                        STOP                                                    
                ENDIF                                                           
		TAB(NOKL)=-100                                                                
                ENR=ENR+1                                                       
                WRITE(9,REC=ENR) TAB                                            
		GO TO 300                                                                     
C                                                                               
C                                                                               
C	CAS NO. 3 : A ET B=0 :                                                        
C       ....................                                                    
C                                                                               
 100	NOK=PKL(1)                                                                 
	VALK=TAB(NOK)                                                                  
	NOL=PKL(2)                                                                     
	VALL=TAB(NOL)                                                                  
	SOM=VALK+VALL                                                                  
	DO 110 II=0,2*ITYPMAX+1                                                        
	SOM=SOM+TAB(PAB(1)+II)+TAB(PAB(2)+II)                                          
	TAB(PAB(1)+II)=0                                                               
	TAB(PAB(2)+II)=0                                                               
 110	CONTINUE                                                                   
        IF(SOM) 150,175,300                                                     
 150            WRITE(SOR,850) IJ,ICAS                                          
                STOP                                                            
 175            NTFCP=NTFCP+1                                                   
		COEFF=COEFF/(TABDENU(VALK)*TABDENU(VALL)*2)                                   
		TAB(10+NDELTA)=3-TAB(10+NDELTA)+TAB(NOK)+TAB(NOL)                             
		TAB(NOK)=-100                                                                 
		TAB(NOL)=-100                                                                 
C                                                                               
C                                                                               
        ENR=ENR+1                                                               
        WRITE(9,REC=ENR) TAB                                                    
C                                                                               
 300    CONTINUE                                                                
C                                                                               
C                                                                               
	NDEB(ITYPINT,ICAS,IJ)=ENR1+1                                                   
	NFIN(ITYPINT,ICAS,IJ)=ENR1+NTFCP                                               
C                                                                               
        WRITE(SOR,901) IJ,CHAR(IJ)                                              
        WRITE(SOR,902) NTFCP,NDEB(ITYPINT,ICAS,IJ),NFIN(ITYPINT,ICAS,IJ)        
 400    CONTINUE                                                                
C                                                                               
C                                                                               
 500    CONTINUE                                                                
C                                                                               
C                                                                               
        RETURN                                                                  
        END                                                                     
C                                                                               
C                                                                               
C                                                                               
C       ****************************************************************        
C       ****************************************************************        
C                                                                               
C                                                                               
        SUBROUTINE CONTRAC(IJ)                                                  
C                                                                               
C       ................................................................        
C                                                                               
C       CETTE SUBROUTINE CONTRACTE LES 'NTERMF' TERMES FINALS POUR UNE          
C	INTEGRALE (IJ) DONNEE, C'EST-A-DIRE QU'ELLE ADDITIONNE LES TER-               
C	MES QUI NE DIFFERENT QUE PAR LE COEFFICIENT.                                  
C                                                                               
C       ................................................................        
C                                                                               
C                                                                               
C	DEFINITIONS :                                                                 
C       ***********                                                             
C                                                                               
C       NTF : NOMBRE DE TERMES FINALS CONTRACTES POUR UNE INTEGRALE             
C	      (IJ) DONNEE.                                                            
C                                                                               
C       INDEX(NTFIN) : TABLEAU CONTENANT LE NUMERO DE L'INDEX (NOINDEX)         
C                      DE CHAQUE TERME FINAL DE TABAUX(NTERMF,NPARTAB).         
C                                                                               
C                                                                               
C	...............................................................               
C                                                                               
C                                                                               
C                                                                               
        IMPLICIT INTEGER (A-Z)                                                  
        DOUBLE PRECISION COEFF,COEFF1,EPS,COEFFI                                
        PARAMETER (ITYPMAX=3,INTMAX=2,NTBFIN=45000)                             
	COMMON/PARGEN/NPARTAB,NDELTA,IJMAX                                             
        COMMON/DERI/TAB(15+ITYPMAX*(2*ITYPMAX+5)),                              
     1              TABAUX(NTBFIN,15+ITYPMAX*(2*ITYPMAX+5))                     
        COMMON/DERIVE/NTF                                                       
        COMMON/DERADD2/NDEB(INTMAX,0:3,(ITYPMAX+1)*(ITYPMAX+2)/2),              
     1                 NFIN(INTMAX,0:3,(ITYPMAX+1)*(ITYPMAX+2)/2)               
	COMMON/AUX/NTERMF                                                              
        COMMON/INTEG/ITYPINT                                                    
        COMMON/ENREG/ENR                                                        
        COMMON/LECRIT/ENT,SOR                                                   
        EQUIVALENCE (TAB(1),COEFF),(TAB1(1),COEFF1)                             
        DIMENSION INDEX(NTBFIN)                                                 
        INTEGER*2 TAB,TABAUX,TAB1(15+ITYPMAX*(2*ITYPMAX+5))                     
C                                                                               
C                                                                               
C                                                                               
        EPS=1.D-6                                                               
C                                                                               
C                                                                               
        DO 5 I=1,NTERMF                                                         
        INDEX(I)=0                                                              
 5      CONTINUE                                                                
C                                                                               
C                                                                               
C                                                                               
        NOINDEX=0                                                               
C                                                                               
C                                                                               
C       ATTRIBUTION D'UN NUMERO D'INDEX A CHACUN DES TERMES FINALS. LES         
C       TERMES QUI NE DIFFERENT QUE PAR LE COEFF PORTERONT LE MEME NUME-        
C       RO D'INDEX : NOINDEX.                                                   
C                                                                               
C                                                                               
        DO 50 N=1,NTERMF-1                                                      
        NOINDEX= NOINDEX+1                                                      
        IF(INDEX(N).NE.0) GO TO 50                                              
        DO 15 J=1,NPARTAB                                                       
        TAB(J)=TABAUX(N,J)                                                      
 15     CONTINUE                                                                
        DO 40 I=N+1,NTERMF                                                      
        IF(INDEX(I).EQ.0)THEN                                                   
                DO 25 K=5,NPARTAB                                               
                IF(TAB(K).NE.TABAUX(I,K)) GO TO 40                              
 25             CONTINUE                                                        
                INDEX(I)=NOINDEX                                                
                INDEX(N)=NOINDEX                                                
        ENDIF                                                                   
 40     CONTINUE                                                                
 50     CONTINUE                                                                
C                                                                               
C                                                                               
C       ADDITION DES TERMES FINALS QUI NE DIFFERENT QUE PAR LE COEFFI-          
C	CIENT ET CONTRACTION DU NOMBRE DE TERMES FINALS (NTERMF DEVIENT               
C	NTF).                                                                         
C	IL EST A NOTER QUE LORSQU'UN TERME A ETE ADDITIONNE, SON INDEX                
C	EST MIS A -1.                                                                 
C                                                                               
C                                                                               
        NTF=0                                                                   
        DO 200 I=1,NTERMF                                                       
        DO 60 K=1,NPARTAB                                                       
 60     TAB(K)=TABAUX(I,K)                                                      
        IF(INDEX(I).EQ.-1) GO TO 200                                            
        IF(INDEX(I).EQ.0) GO TO 100                                             
        DO 80 J=I+1,NTERMF                                                      
        IF(INDEX(I).EQ.INDEX(J))THEN                                            
                DO 65 K=1,4                                                     
 65             TAB1(K)=TABAUX(J,K)                                             
                COEFF=COEFF+COEFF1                                              
                INDEX(J)=-1                                                     
        ENDIF                                                                   
 80     CONTINUE                                                                
        IF(DABS(COEFF).LT.EPS) GO TO 200                                        
 100    NTF=NTF+1                                                               
        ENR=ENR+1                                                               
        WRITE(9,REC=ENR) TAB                                                    
 200    CONTINUE                                                                
C                                                                               
C                                                                               
C                                                                               
        RETURN                                                                  
        END                                                                     
C                                                                               
C                                                                               
C                                                                               
C       ****************************************************************        
C       ****************************************************************        
C                                                                               
C                                                                               
C                                                                               
        SUBROUTINE DERIV(ITYPE)                                                 
C                                                                               
C       ................................................................        
C                                                                               
C	CETTE SUBROUTINE GENERE LES EXPRESSIONS DES DERIVEES POUR CHA-                
C	CUN DES TERMES DERIVABLES DU TABLEAU : TAB(NPARTAB)                           
C                                                                               
C	LE TYPE DE DERIVATION DOIT ETRE TRANSMIS EN ARGUMENT.                         
C                                                                               
C	................................................................              
C                                                                               
C	DEFINITIONS :                                                                 
C	***********	                                                                  
C                                                                               
C	TAB(NPARTAB) :                                                                
C                                                                               
C       [1-4]  5  6   7   NDELTA  +1+1+1  +1+1   ...   +1 +1   ...    +1        
C	COEFF 2A 2B (A+B) DELTAS I(K,L,N) PL A AI...AK AOP B BI...BK BOP              
C                                                                               
C                                                                               
C       NOMBRE DE CASES TOTALES = NPARTAB                                       
C	NDELTA=ITYPMAX*(2*ITYPMAX+1)                                                  
C	AI...AK = 2*ITYPMAX  (IDEM POUR BI...BK)                                      
C                                                                               
C			NOTONS QUE 2A ET 2B SYMBOLISENT 2*ALPHA ET                                  
C                       2*BETA RESPECTIVEMENT.                                  
C                                                                               
C	................................................................              
C                                                                               
C                                                                               
C                                                                               
        IMPLICIT INTEGER (A-Z)                                                  
        DOUBLE PRECISION COEFF1,COEFF                                           
        PARAMETER (ITYPMAX=3,INTMAX=2,NTBFIN=45000)                             
	COMMON/PARGEN/NPARTAB,NDELTA,IJMAX                                             
        COMMON/DERI/TAB(15+ITYPMAX*(2*ITYPMAX+5)),                              
     1              TABAUX(NTBFIN,15+ITYPMAX*(2*ITYPMAX+5))                     
        COMMON/DERADD2/NDEB(INTMAX,0:3,(ITYPMAX+1)*(ITYPMAX+2)/2),              
     1                 NFIN(INTMAX,0:3,(ITYPMAX+1)*(ITYPMAX+2)/2)               
	COMMON/AUX/NTERMF                                                              
        COMMON/INTEG/ITYPINT                                                    
        COMMON/LECRIT/ENT,SOR                                                   
        INTEGER*2 TAB,TABAUX,TAB1(15+ITYPMAX*(2*ITYPMAX+5))                     
        EQUIVALENCE (TAB1(1),COEFF1),(TAB(1),COEFF)                             
C                                                                               
C                                                                               
C	DIVISION PAR 2*ALPHA SI L'ON DERIVE A GAUCHE OU PAR 2*BETA SI                 
C	L'ON DERIVE A DROITE.                                                         
C                                                                               
C                                                                               
        IF(ITYPE.LE.ITYPMAX)THEN                                                
                TAB(5)=TAB(5)-1                                                 
                      ELSE                                                      
                TAB(6)=TAB(6)-1                                                 
        ENDIF                                                                   
C                                                                               
C                                                                               
C	ON DEFINIT LE NOMBRE DE PARAMETRES DERIVABLES PAR : INPAR ET                  
C 	CES PARAMETRES SERONT CONSIDERES DANS L'ORDRE SUIVANT :                      
C                                                                               
C	PL, A, AI...AK, AOP, B, BI...BK, BOP, I(K,L,N) ET EXP.                        
C                                                                               
C	INPAR=1+1+(2*ITYPMAX)+1+1+(2*ITYPMAX)+1+1+1                                   
C	     =7+4*ITYPMAX                                                             
C                                                                               
C	CELA NOUS PERMET DONC DE TRAITER LES PARAMETRES DERIVABLES DIF-               
C	FERENTS DE ZERO DANS UN PREMIER TEMPS. LES 'INPAR' ET 'INPAR-1'               
C	(L'INTEGRALE IKLN ET L'EXPONENTIELLE RESPECTIVEMENT) SONT ALORS               
C	CONSIDERES POUR CHAQUE INTEGRALE (IJ), CES DERNIERS DEVANT ETRE               
C	DERIVES DANS TOUS LES CAS.                                                    
C	ON ASSIGNE A CHACUN DES PARAMETRES DERIVABLES (NPAR=1,INPAR), UN              
C	NOUVEAU NUMERO (NPARP) SERVANT A LES REGROUPER ET CE DERNIER SE-              
C	RA UTILISE DANS LE " GO TO CALCULE ".                                         
C                                                                               
C                                                                               
	INPAR=7+4*ITYPMAX                                                              
C                                                                               
C                                                                               
        DO 2000 NPAR=1,INPAR                                                    
	NPARP=NPAR                                                                     
        IF(NPAR.GT.INPAR-2)THEN                                                 
		NPARP=7+NPAR-INPAR                                                            
		GO TO 15                                                                      
	ENDIF                                                                          
        IF(TAB(10+NDELTA+NPAR).EQ.0) GO TO 2000                                 
C                                                                               
C                                                                               
	IF(NPAR.EQ.1) GO TO 15                                                         
C                                                                               
	IF(ITYPE.LE.ITYPMAX)THEN                                                       
		IF(NPAR.GT.(2*ITYPMAX+3)) GO TO 2000                                          
		IF(NPAR.GT.2 .AND. NPAR.LE.(2*ITYPMAX+3)) NPARP=3                             
			    ELSE                                                                     
		IF(NPAR.LT.(2*ITYPMAX+4)) GO TO 2000                                          
		IF(NPAR.EQ.(2*ITYPMAX+4)) NPARP=4                                             
		IF(NPAR.GT.(2*ITYPMAX+4)) NPARP=5                                             
	ENDIF                                                                          
C                                                                               
C                                                                               
C	DEFINITIONS DES NUMEROS DES CASES DES PARAMETRES : A, B, AI,                  
C	BI, PL ET I(KLN) DANS LE TABLEAU TAB(NPARTAB). LES (AI) SONT                  
C	DEFINIS PAR RAPPORT A : "A" ET AU "TYPE DE DERIVATION". IDEM                  
C	POUR LES (BI).                                                                
C                                                                               
C                                                                               
 15	TABA=12+NDELTA                                                              
	TABB=14+NDELTA+2*ITYPMAX                                                       
	TABAI=TABA+ITYPE                                                               
	TABBI=TABB+ITYPE                                                               
	TABPL=11+NDELTA                                                                
	TABK=8+NDELTA                                                                  
	TABL=9+NDELTA                                                                  
	TABN=10+NDELTA                                                                 
C                                                                               
C                                                                               
        DO 20 I=1,NPARTAB                                                       
        TAB1(I)=TAB(I)                                                          
 20     CONTINUE                                                                
C                                                                               
C                                                                               
	GO TO (500,100,300,200,400,900,1000) NPARP                                     
C                                                                               
C	      (PL  A   AI  B   BI  IKLN EXP)                                          
C       ................................................................        
C                                                                               
C                                                                               
C       DERIVEE DE A                                                            
C       ------------                                                            
C                                                                               
 100    NTERMP=0                                                                
        COEFF1=COEFF1*TAB1(TABA)                                                
        TAB1(5)=TAB1(5)+2                                                       
        TAB1(7)=TAB1(7)-1                                                       
        TAB1(TABAI)=TAB1(TABAI)+1                                               
        TAB1(TABA)=TAB1(TABA)-2                                                 
        GO TO 1500                                                              
C                                                                               
C                                                                               
C       DERIVEE DE B                                                            
C       ------------                                                            
C                                                                               
 200    NTERMP=0                                                                
        COEFF1=COEFF1*TAB1(TABB)                                                
        TAB1(6)=TAB1(6)+2                                                       
        TAB1(7)=TAB1(7)-1                                                       
        TAB1(TABBI)=TAB1(TABBI)+1                                               
        TAB1(TABB)=TAB1(TABB)-2                                                 
        GO TO 1500                                                              
C                                                                               
C                                                                               
C       DERIVEES DES A(I,J,K,L,M,N) ET B(I,J,K,L,M,N)                           
C       ---------------------------------------------                           
C                                                                               
 300                    K=NPAR-2                                                
			GO TO 450                                                                    
C                                                                               
C                                                                               
 400                    K=NPAR-(4+2*ITYPMAX)                                    
C                                                                               
C                                                                               
 450	NTERMP=0                                                                   
    	COEFF1=TAB1(10+NDELTA+NPAR)*COEFF1                                         
	TAB1(10+NDELTA+NPAR)=TAB1(10+NDELTA+NPAR)-1                                    
        IF(K.EQ.ITYPE) GO TO 1500                                               
	IF(K.GT.ITYPE)THEN                                                             
		TDELTA=(K-1)*(K-2)/2+ITYPE+7                                                  
		TAB1(TDELTA)=1                                                                
                      ELSE                                                      
		TDELTA=(ITYPE-1)*(ITYPE-2)/2+K+7                                              
		TAB1(TDELTA)=1                                                                
	ENDIF                                                                          
	GO TO 1500                                                                     
C                                                                               
C                                                                               
C       DERIVEES DE P1, P2 ET P3                                                
C       ************************                                                
C                                                                               
 500    PL=TAB(TABPL)                                                           
    	GO TO (600,700,800) PL                                                     
C                                                                               
C       DERIVEE DE P1                                                           
C       -------------                                                           
C                                                                               
C       1ER TERME                                                               
C                                                                               
 600    NTERMP=1                                                                
        TAB1(5)=TAB1(5)+1                                                       
        TAB1(6)=TAB1(6)+1                                                       
        IF(ITYPE.GT.ITYPMAX)THEN                                                
                TAB1(TABAI)=TAB1(TABAI)+1                                       
                      ELSE                                                      
                TAB1(TABBI)=TAB1(TABBI)+1                                       
        ENDIF                                                                   
        TAB1(7)=TAB1(7)-1                                                       
        TAB1(TABA)=TAB1(TABA)-1                                                 
        TAB1(TABB)=TAB1(TABB)-1                                                 
        TAB1(TABPL)=0                                                           
        GO TO 1500                                                              
C                                                                               
C       2EME TERME                                                              
C                                                                               
 610    COEFF1=-COEFF1                                                          
        TAB1(7)=TAB1(7)-1                                                       
        IF(ITYPE.LE.ITYPMAX)THEN                                                
                TAB1(5)=TAB1(5)+2                                               
                TAB1(TABA)=TAB1(TABA)-2                                         
                TAB1(TABAI)=TAB1(TABAI)+1                                       
                      ELSE                                                      
                TAB1(6)=TAB1(6)+2                                               
                TAB1(TABB)=TAB1(TABB)-2                                         
                TAB1(TABBI)=TAB1(TABBI)+1                                       
        ENDIF                                                                   
C	TAB1(TABPL)=1                                                                 
        GO TO 1500                                                              
C                                                                               
C                                                                               
C       DERIVEE DE P2                                                           
C       -------------                                                           
C                                                                               
C       1ER TERME                                                               
C                                                                               
 700    NTERMP=3                                                                
        COEFF1=COEFF1*3                                                         
        TAB1(5)=TAB1(5)+1                                                       
        TAB1(6)=TAB1(6)+1                                                       
        TAB1(7)=TAB1(7)-1                                                       
        TAB1(TABA)=TAB1(TABA)-1                                                 
        TAB1(TABB)=TAB1(TABB)-1                                                 
        TAB1(TABPL)=1                                                           
        IF(ITYPE.LE.ITYPMAX)THEN                                                
                TAB1(TABBI)=TAB1(TABBI)+1                                       
                      ELSE                                                      
                TAB1(TABAI)=TAB1(TABAI)+1                                       
        ENDIF                                                                   
        GO TO 1500                                                              
C                                                                               
C       2EME TERME                                                              
C                                                                               
 710    COEFF1=COEFF1*-2                                                        
        TAB1(7)=TAB1(7)-1                                                       
C       TAB1(TABPL)=2                                                           
        IF(ITYPE.LE.ITYPMAX)THEN                                                
                TAB1(TABAI)=TAB1(TABAI)+1                                       
                TAB1(5)=TAB1(5)+2                                               
                TAB1(TABA)=TAB1(TABA)-2                                         
                      ELSE                                                      
                TAB1(TABBI)=TAB1(TABBI)+1                                       
                TAB1(6)=TAB1(6)+2                                               
                TAB1(TABB)=TAB1(TABB)-2                                         
        ENDIF                                                                   
        GO TO 1500                                                              
C                                                                               
C       3EME TERME                                                              
C                                                                               
 720    COEFF1=-COEFF1                                                          
        TAB1(TABPL)=0                                                           
        TAB1(7)=TAB1(7)-1                                                       
        IF(ITYPE.LE.ITYPMAX)THEN                                                
                TAB1(TABAI)=TAB1(TABAI)+1                                       
                TAB1(5)=TAB1(5)+2                                               
                TAB1(TABA)=TAB1(TABA)-2                                         
                      ELSE                                                      
                TAB1(TABBI)=TAB1(TABBI)+1                                       
                TAB1(6)=TAB1(6)+2                                               
                TAB1(TABB)=TAB1(TABB)-2                                         
        ENDIF                                                                   
        GO TO 1500                                                              
C                                                                               
C                                                                               
C       DERIVEE DE P3                                                           
C       -------------                                                           
C                                                                               
C       1ER TERME                                                               
C                                                                               
 800    NTERMP=6                                                                
        COEFF1=COEFF1*-3                                                        
        TAB1(7)=TAB1(7)-1                                                       
C       TAB1(TABPL)=3                                                           
        IF(ITYPE.LE.ITYPMAX)THEN                                                
		TAB1(TABAI)=TAB1(TABAI)+1                                                     
                TAB1(5)=TAB1(5)+2                                               
                TAB1(TABA)=TAB1(TABA)-2                                         
                      ELSE                                                      
                TAB1(TABBI)=TAB1(TABBI)+1                                       
                TAB1(6)=TAB1(6)+2                                               
                TAB1(TABB)=TAB1(TABB)-2                                         
        ENDIF                                                                   
        GO TO 1500                                                              
C                                                                               
C       2EME TERME                                                              
C                                                                               
 810    COEFF1=COEFF1*5                                                         
        TAB1(TABPL)=2                                                           
        TAB1(5)=TAB1(5)+1                                                       
        TAB1(6)=TAB1(6)+1                                                       
        TAB1(7)=TAB1(7)-1                                                       
        TAB1(TABA)=TAB1(TABA)-1                                                 
        TAB1(TABB)=TAB1(TABB)-1                                                 
        IF(ITYPE.GT.ITYPMAX)THEN                                                
                TAB1(TABAI)=TAB1(TABAI)+1                                       
                      ELSE                                                      
                TAB1(TABBI)=TAB1(TABBI)+1                                       
        ENDIF                                                                   
        GO TO 1500                                                              
C                                                                               
C       3EME TERME                                                              
C                                                                               
 820    COEFF1=COEFF1*-3                                                        
        TAB1(TABPL)=1                                                           
        TAB1(7)=TAB1(7)-1                                                       
        IF(ITYPE.LE.ITYPMAX)THEN                                                
                TAB1(TABAI)=TAB1(TABAI)+1                                       
                TAB1(5)=TAB1(5)+2                                               
                TAB1(TABA)=TAB1(TABA)-2                                         
                      ELSE                                                      
                TAB1(TABBI)=TAB1(TABBI)+1                                       
                TAB1(6)=TAB1(6)+2                                               
                TAB1(TABB)=TAB1(TABB)-2                                         
        ENDIF                                                                   
        GO TO 1500                                                              
C                                                                               
C       4EME TERME                                                              
C                                                                               
 830    TAB1(5)=TAB1(5)+1                                                       
        TAB1(6)=TAB1(6)+1                                                       
        TAB1(7)=TAB1(7)-1                                                       
        TAB1(TABA)=TAB1(TABA)-1                                                 
        TAB1(TABB)=TAB1(TABB)-1                                                 
        TAB1(TABPL)=0                                                           
        IF(ITYPE.GT.ITYPMAX)THEN                                                
                TAB1(TABAI)=TAB1(TABAI)+1                                       
                      ELSE                                                      
                TAB1(TABBI)=TAB1(TABBI)+1                                       
        ENDIF                                                                   
        GO TO 1500                                                              
C                                                                               
C                                                                               
C       DERIVEE DE I KLN                                                        
C       ****************                                                        
C                                                                               
C       1ER TERME                                                               
C                                                                               
 900    NTERMP=10                                                               
        TAB1(7)=TAB1(7)-1                                                       
        TAB1(TABN)=TAB1(TABN)-1                                                 
        IF(ITYPE.LE.ITYPMAX)THEN                                                
                TAB1(TABAI)=TAB1(TABAI)+1                                       
                TAB1(5)=TAB1(5)+2                                               
                TAB1(TABA)=TAB1(TABA)-1                                         
                TAB1(TABK)=TAB1(TABK)+1                                         
                      ELSE                                                      
                TAB1(TABBI)=TAB1(TABBI)+1                                       
                TAB1(6)=TAB1(6)+2                                               
                TAB1(TABB)=TAB1(TABB)-1                                         
                TAB1(TABL)=TAB1(TABL)+1                                         
        ENDIF                                                                   
        GO TO 1500                                                              
C                                                                               
C       2EME TERME                                                              
C                                                                               
 910    IF((ITYPE.LE.ITYPMAX .AND. TAB1(TABK).EQ.0).OR.                         
     1     (ITYPE.GT.ITYPMAX .AND. TAB1(TABL).EQ.0)) GO TO 2000                 
        TAB1(7)=TAB1(7)-1                                                       
        IF(ITYPE.LE.ITYPMAX)THEN                                                
                TAB1(TABAI)=TAB1(TABAI)+1                                       
                COEFF1=COEFF1*TAB1(TABK)                                        
                TAB1(5)=TAB1(5)+2                                               
                TAB1(TABA)=TAB1(TABA)-2                                         
                      ELSE                                                      
                TAB1(TABBI)=TAB1(TABBI)+1                                       
                COEFF1=COEFF1*TAB1(TABL)                                        
                TAB1(6)=TAB1(6)+2                                               
                TAB1(TABB)=TAB1(TABB)-2                                         
        ENDIF                                                                   
        GO TO 1500                                                              
C                                                                               
C                                                                               
C       DERIVEES DES EXP(-ALPHA*(AC**2) ET                                      
C                   EXP(-BETA*(BC**2)                                           
C       **********************************                                      
C                                                                               
C                                                                               
 1000   NTERMP=0                                                                
        COEFF1=-COEFF1                                                          
        IF(ITYPE.LE.ITYPMAX)THEN                                                
                TAB1(TABAI)=TAB1(TABAI)+1                                       
                TAB1(5)=TAB1(5)+1                                               
                      ELSE                                                      
                TAB1(TABBI)=TAB1(TABBI)+1                                       
                TAB1(6)=TAB1(6)+1                                               
        ENDIF                                                                   
        GO TO 1500                                                              
C                                                                               
C                                                                               
C	NTERMF=NOMBRE DE TERMES FINALS POUR UNE INTEGRALE 'IJ' DONNEE.                
C                                                                               
C	NTERMP= VARIABLE UTILISEE DANS LE "GO TO CALCULE" PERMETTANT                  
C	D'ADDITIONNER TOUS LES TERMES D'UN MEME PARAMETRE DERIVABLE.                  
C                                                                               
C			EX : I(KLN) : 2 TERMES                                                      
C	LES RESULTATS DE LA SUBROUTINE 'DERIV' SONT MIS DANS LE TA-                   
C	BLEAU : TABAUX(NTERMF,NPARTAB).                                               
C                                                                               
C                                                                               
 1500   NTERMP=NTERMP+1                                                         
        NTERMF=NTERMF+1                                                         
        IF(NTERMF.GT.NTBFIN) GO TO 2100                                         
        DO 1800 I=1,NPARTAB                                                     
        TABAUX(NTERMF,I)=TAB1(I)                                                
        TAB1(I)=TAB(I)                                                          
 1800   CONTINUE                                                                
        GO TO (2000,610,2000,710,720,2000,810,820,830,2000,                     
     1         910,2000) NTERMP                                                 
C                                                                               
C                                                                               
 2000   CONTINUE                                                                
C                                                                               
C                                                                               
C                                                                               
        RETURN                                                                  
 2100   PRINT*,'NTBFIN INSUFFISANT, NTERMF= ',NTERMF                            
        STOP                                                                    
        END                                                                     
                                                                                
C                                                                               
C                                                                               
C       ****************************************************************        
C       ****************************************************************        
C                                                                               
C                                                                               
        SUBROUTINE DERIVC(L)                                                    
C                                                                               
C       ................................................................        
C                                                                               
C       CETTE SUBROUTINE GENERE LES TERMES DE L'INTEGRALE :                     
C                     < S / X,Y,Z/RC**3 PLC / S >     (1)                       
C                 CF. DEMONSTRATION M.F.                                        
C                                                                               
C       ................................................................        
C                                                                               
        IMPLICIT INTEGER (A-Z)                                                  
        DOUBLE PRECISION COEFF1,COEFF2,SIGN                                     
        PARAMETER (ITYPMAX=3,INTMAX=2,NTBFIN=45000)                             
	COMMON/PARGEN/NPARTAB,NDELTA,IJMAX                                             
        COMMON/TABPAR/TABA,TABB,TABPL,TABK,TABL,TABN,TABAOP,TABBOP              
        COMMON/DERIVE/NTF                                                       
        COMMON/DERADD2/NDEB(INTMAX,0:3,(ITYPMAX+1)*(ITYPMAX+2)/2),              
     1                 NFIN(INTMAX,0:3,(ITYPMAX+1)*(ITYPMAX+2)/2)               
        COMMON/DERI/TAB(15+ITYPMAX*(2*ITYPMAX+5)),                              
     1              TABAUX(NTBFIN,15+ITYPMAX*(2*ITYPMAX+5))                     
        COMMON/LECRIT/ENT,SOR                                                   
        COMMON/ENREG/ENR                                                        
        COMMON/INTEG/ITYPINT                                                    
        INTEGER*2 TAB1(15+ITYPMAX*(2*ITYPMAX+5)),TAB2(15+ITYPMAX*               
     1            (2*ITYPMAX+5)),TAB,TABAUX                                     
        EQUIVALENCE (TAB1(1),COEFF1),(TAB2(1),COEFF2)                           
C                                                                               
C                                                                               
 4900   FORMAT(//,28X,'DERIVEE PAR RAPPORT A | C',/,28X,25('-'))                
 4901   FORMAT(///,1X,'INTEGRALE NO :     1      <S/S>',/)                      
 4902   FORMAT(1X,'NOMBRE DE TERMES FINALS (2*L+2) : ',I5,/,1X,                 
     1            'NDEB = ',I10,10X,'NFIN = ',I10,//)                           
C4903	FORMAT(/,1x,'2ALPHA 2BETA ALPHA+BETA DELTAS I(K,L,N) PL A'                
C    1         ' AI...K AOP B BI...BK BOP COEFF',//)                            
C                                                                               
C                                                                               
C                                                                               
        COEFF1=DFLOAT(2*L+1)                                                    
        SIGN=-1.0                                                               
        LL=L+1                                                                  
        NO2=ENR                                                                 
C                                                                               
C                                                                               
C                                                                               
        DO 400 NTERM=1,L+1                                                      
C                                                                               
C                                                                               
        IF(NO2.GT.NTBFIN) GO TO 9999                                            
C                                                                               
C                                                                               
        NO2=NO2+2                                                               
        NO1=NO2-1                                                               
C                                                                               
C                                                                               
        DO 20 I=5,NPARTAB                                                       
        TAB1(I)=0                                                               
 20     CONTINUE                                                                
        LL=LL-1                                                                 
        TAB1(TABN)=2                                                            
        TAB1(TABPL)=LL                                                          
        SIGN=-SIGN                                                              
        COEFF1=COEFF1*SIGN                                                      
        DO 100 K=1,NPARTAB                                                      
        TAB2(K)=TAB1(K)                                                         
 100    CONTINUE                                                                
C                                                                               
C                                                                               
        IF(COEFF1.GT.0.0D0)THEN                                                 
		MI=L+1                                                                        
		MP=L                                                                          
                          ELSE                                                  
		MI=L                                                                          
		MP=L+1                                                                        
        ENDIF                                                                   
C                                                                               
C                                                                               
        TAB1(5)=1                                                               
	TAB1(TABAOP)=1                                                                 
	TAB1(TABA)=-1                                                                  
	TAB1(TABK)=MI                                                                  
 	TAB1(TABL)=MP                                                                 
C                                                                               
C                                                                               
        TAB2(6)=1                                                               
	TAB2(TABBOP)=1                                                                 
	TAB2(TABB)=-1                                                                  
	TAB2(TABK)=MP                                                                  
	TAB2(TABL)=MI                                                                  
C                                                                               
C                                                                               
        DO 250 J=1,NPARTAB                                                      
        TABAUX(NO1,J)=TAB1(J)                                                   
        TABAUX(NO2,J)=TAB2(J)                                                   
 250    CONTINUE                                                                
C                                                                               
C                                                                               
        COEFF1=DABS(COEFF1)-2.0D0                                               
C                                                                               
C                                                                               
 400	CONTINUE                                                                   
C                                                                               
C                                                                               
C                                                                               
        NDEB(ITYPINT,0,1)=ENR+1                                                 
        NFIN(ITYPINT,0,1)=ENR+(2*L+2)                                           
C                                                                               
        WRITE(SOR,4900)                                                         
        WRITE(SOR,4901)                                                         
        WRITE(SOR,4902) (2*L+2),NDEB(ITYPINT,0,1),NFIN(ITYPINT,0,1)             
C       WRITE(SOR,4903)                                                         
        DO 7000 I=NDEB(ITYPINT,0,1),NFIN(ITYPINT,0,1)                           
        DO 6000 J=1,NPARTAB                                                     
 6000   TAB1(J)=TABAUX(I,J)                                                     
        ENR=ENR+1                                                               
        WRITE(9,REC=ENR)TAB1                                                    
 7000   CONTINUE                                                                
C                                                                               
C                                                                               
C                                                                               
        RETURN                                                                  
 9999   PRINT*,'L=',L,'NTBFIN INSUFFISANT, NO2=',NO2                            
        STOP                                                                    
        END                                                                     
C                                                                               
C                                                                               
C       ****************************************************************        
C       ****************************************************************        
C                                                                               
C                                                                               
C                                                                               
        SUBROUTINE ECRIT                                                        
C                                                                               
C       ................................................................        
C                                                                               
C       CETTE SUBROUTINE A POUR FONCTION D'ECRIRE LES RESULTATS DES DERI-       
C       VEES POUR CHACUN DES PL SUR LA FILE 10 DE LA FACON SUIVANTE :           
C                                                                               
C       CAS 1/R**N : NDEB, NFIN                                                 
C       **********   CASP=0,1,2,3                                               
C                                                                               
C       CAS X,Y,Z/R**3 : NDEB, NFIN                                             
C       **************   CASP=0,1,2,3                                           
C                                                                               
C                          TOTAL : 8 ENREGISTREMENTS                            
C                                                                               
C       ................................................................        
C                                                                               
C                                                                               
C                                                                               
        IMPLICIT INTEGER(A-Z)                                                   
        PARAMETER (ITYPMAX=3,INTMAX=2,NTBFIN=45000,lbloc=16384)
        COMMON/PARGEN/NPARTAB,NDELTA,IJMAX                                      
        COMMON/INTEG/ITYPINT,ICAMAX                                             
        COMMON/DERADD2/NDEB(INTMAX,0:3,(ITYPMAX+1)*(ITYPMAX+2)/2),              
     1                 NFIN(INTMAX,0:3,(ITYPMAX+1)*(ITYPMAX+2)/2)               
        COMMON/LECRIT/ENT,SOR                                                   
        INTEGER*2 TAB(NTBFIN*(15+ITYPMAX*(2*ITYPMAX+5)))                        
        COMMON/DERI/TAB                                                         
C                                                                               
C                                                                               
C                                                                               
        NPAR=NTBFIN*(15+ITYPMAX*(2*ITYPMAX+5))                                  
        WRITE(10) (((NDEB(I,J,K),NFIN(I,J,K),K=1,IJMAX),J=0,ICAMAX),            
     1               I=1,ITYPINT)                                               
C                                                                               
C                                                                               
        K=1                                                                     
C                                                                               
C                                                                               
        DO 10 I=1,NFIN(ITYPINT,ICAMAX,IJMAX)                                    
        READ(9,REC=I) (TAB(J),J=K,(K+NPARTAB-1))                                
        K=K+NPARTAB                                                             
        IF(K.GT.NPAR)THEN                                                       
		WRITE(SOR,*) 'TABLEAU DANS ECRIT TROP PETIT'                                  
                STOP                                                            
        ENDIF                                                                   
 10     CONTINUE                                                                
C                                                                               
C                                                                               
        nbloc=(k+npartab-1)/lbloc
        if((k+NPARTAB-1).ne.lbloc*nbloc) nbloc=nbloc+1
        do 20 i=1,nbloc
        mdeb=(i-1)*lbloc+1
        mfin=mdeb+lbloc
        if(i.eq.nbloc) mfin=k+npartab-1
        WRITE(10) (TAB(L),L=mdeb,mfin)
   20   continue                                           
C       WRITE(SOR,*) (TAB(L),L=1,K+NPARTAB-1)                                   
C                                                                               
C                                                                               
        RETURN                                                                  
        END                                                                     
C                                                                               
C                                                                               
C                                                                               
C	****************************************************************              
C	****************************************************************              
C                                                                               
C                                                                               
C                                                                               
        SUBROUTINE INTEGEX                                                      
C                                                                               
C       ...............................................................         
C                                                                               
C	CETTE SUBROUTINE GENERE LES EXPRESSIONS DES DIFFERENTS TYPES D'               
C       INTEGRALES A CALCULER ( PAR COUCHE ). C'EST-A-DIRE QUE POUR UNE         
C	INTEGRALE " IJ ", ON DEFINIT LES PARAMETRES SUIVANTS : LE NUME-               
C       RO DE L'INTEGRALE A DERIVER POUR OBTENIR CELLE DESIREE, LE(S)           
C	NUMERO(S) DE(S) INTEGRALE(S) A ADDITIONNER ET LEUR(S) DELTA(S)                
C	CORRESPONDANT(S).                                                             
C	CES PARAMETRES SONT CONTENUS DANS LE TABLEAU PARIJ(IJMAX,NPMAX).              
C                                                                               
C	...............................................................               
C                                                                               
C                                                                               
C	CHAQUE INTEGRALE " IJ " EST REPRESENTEE PAR UN COUPLE D'INDICE                
C	(I,J) [ CF. LES BOUCLES DO 1000 ET DO 950 ], PAR EXEMPLE :                    
C                                                                               
C			(1,1) : <S/S>                                                               
C			(2,1) : <P/S>                                                               
C			(2,2) : <P/P>                                                               
C			(3,1) : <D/S>                                                               
C			(3,2) : <D/P>                                                               
C			 ...     ...                                                                
C                                                                               
C                                                                               
C     	DEFINITIONS :                                                            
C       ***********                                                             
C                                                                               
C	ITYPMAX :  LE TYPE MAXIMAL DES ORBITALES ETUDIEES, PAR EX. :                  
C                                                                               
C                       " S " , ITYPMAX=0                                       
C                       " P " , ITYPMAX=1                                       
C                       " D " , ITYPMAX=2                                       
C                       " F " , ITYPMAX=3                                       
C                       " G " , ITYPMAX=4                                       
C                        ...       ...                                          
C                                                                               
C       IJMAX : NOMBRE TOTAL D'INTEGRALES .                                     
C  			IJMAX=(ITYPMAX+1)*(ITYPMAX+2)/2                                           
C                                                                               
C	NPMAX : LE NOMBRE MAXIMAL DE PARAMETRES QUE CONTIENT LE TABLEAU               
C               PARIJ(IJMAX,NPMAX).                                             
C 			NPMAX=3+(ITYPMAX-1)                                                        
C			     =2+ITYPMAX                                                             
C                                                                               
C	PARIJ(IJMAX,NPMAX) : POUR CHACUNE DES INTEGRALES " IJ ", CE TA-               
C                            BLEAU CONTIENT LES PARAMETRES SUIVANTS :           
C                                                                               
C  		PARIJ(IJ,1) : LE NUMERO DE L'INTEGRALE A DERIVER                           
C                             POUR OBTENIR L'INTEGRALE " IJ " DE-               
C                             SIREE.                                            
C               PARIJ(IJ,2) : LE TYPE DE DERIVATION ( A1,2,3...K                
C                             OU B1,2,3...K).                                   
C 		PARIJ(IJ,3) : LE NUMERO DE L'INTEGRALE A ADDITIONNER.                       
C		PARIJ(IJ,3+M...NPMAX) : LE(S) NUMERO(S) DE(S) DELTA(S).                      
C		                                                                             
C                                                                               
C	...............................................................               
C                                                                               
C                                                                               
C                                                                               
    	IMPLICIT INTEGER(A-Z)                                                      
   	PARAMETER (ITYPMAX=3)                                                       
	COMMON/PARGEN/NPARTAB,NDELTA,IJMAX                                             
        COMMON/INT/PARIJ((ITYPMAX+1)*(ITYPMAX+2)/2,2+ITYPMAX)                   
	COMMON/LECRIT/ENT,SOR                                                          
	CHARACTER*5 CHAR(10)                                                           
        DATA CHAR/'<S/S>','<P/S>','<P/P>','<D/S>','<D/P>','<D/D>',              
     1            '<F/S>','<F/P>','<F/D>','<F/F>'/                              
C                                                                               
C                                                                               
 4900	FORMAT(5X,'PARAMETRES (PARIJ(IJMAX,NPMAX)) POUR CHAQUE'                   
     1         ' INTEGRALE  IJ',/,5X,57('*'),///)                               
 5000	FORMAT(1X,//,'IJ= ',I5,5X,A5,/,5X,'ON DERIVE : ',I5,5X,                   
     1        'PAR RAPPORT A :',I5,5X,/,5X,'ON ADDITIONNE LE(S)'                
     1        'TERME(S) :',I5,/,5X,'LE(S) DELTA(S) :')                          
 5001   FORMAT(5X,15I5)                                                         
C                                                                               
C                                                                               
	WRITE(SOR,4900)                                                                
C                                                                               
C                                                                               
        NPMAX=2+ITYPMAX                                                         
C                                                                               
	DO 200 I=1,IJMAX                                                               
	DO 100 J=1,NPMAX                                                               
	PARIJ(I,J)=0                                                                   
 100	CONTINUE                                                                   
 200	CONTINUE                                                                   
C                                                                               
C                                                                               
	IJ=1                                                                           
	DO 1000 I=2,ITYPMAX+1                                                          
 	DO 950 J=1,I                                                                  
	IJ=IJ+1                                                                        
C                                                                               
C	ON DERIVE A GAUCHE QUAND J=1 ( 1S ), SINON ON DERIVE A DROITE                 
C	DANS TOUS LES AUTRES CAS.                                                     
C                                                                               
	IF(J.EQ.1)THEN                                                                 
		PARIJ(IJ,1)=(I-1)*(I-2)/2+J                                                   
		PARIJ(IJ,2)=I-1                                                               
C                                                                               
C	ON N'A DES TERMES SUPPLEMENTAIRES QU'A PARTIR DES " D ", C'EST-               
C	A-DIRE QUAND (I-2)>0 LORSQUE L'ON DERIVE A GAUCHE OU QUAND                    
C 	(J-2)>0 LORSQUE L'ON DERIVE A DROITE.                                        
C                                                                               
		IF((I-2).GT.0)THEN                                                            
			PARIJ(IJ,3)=(I-2)*(I-3)/2+J                                                  
			AK=I-1                                                                       
			DO 400 M=1,AK-1                                                              
			PARIJ(IJ,3+M)=7+((AK-1)*(AK-2)/2+M)                                          
 400			CONTINUE                                                                 
            	ENDIF                                                              
                  ELSE                                                          
		PARIJ(IJ,1)=(I*(I-1)/2+J)-1                                                   
		PARIJ(IJ,2)=ITYPMAX+J-1                                                       
		IF((J-2).GT.0)THEN                                                            
 			PARIJ(IJ,3)=IJ-2                                                            
			AK=ITYPMAX+J-1                                                               
			M=0                                                                          
			DO 600 L=ITYPMAX+1,AK-1                                                      
			M=M+1                                                                        
			PARIJ(IJ,3+M)=7+((AK-1)*(AK-2)/2+L)                                          
 600			CONTINUE                                                                 
		ENDIF                                                                         
	ENDIF                                                                          
C                                                                               
C                                                                               
	WRITE(SOR,5000) IJ,CHAR(IJ),PARIJ(IJ,1),PARIJ(IJ,2),PARIJ(IJ,3)                
        WRITE(SOR,5001) (PARIJ(IJ,K),K=4,2+ITYPMAX)                             
C                                                                               
C                                                                               
 950    CONTINUE                                                                
 1000   CONTINUE                                                                
        RETURN                                                                  
        END                                                                     
C                                                                               
C                                                                               
C       ****************************************************************        
C       ****************************************************************        

