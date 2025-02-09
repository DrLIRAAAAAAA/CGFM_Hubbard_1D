!*
	INCLUDE 'MODALL-F'
!!!
!*
!* ############################################################
!*               	BOGOTA  01/11/2019
!* ############################################################
!*
!*  H= -t\sum<> [c^{+}c +h.c. ] + U\sum_i X_i
!*
!*       * HUBBARD MODEL- ATOMIC METHOD
!*#############################################################
	PROGRAM AMHUBBARD
	USE PRINCIPAL
	IMPLICIT NONE
        EXTERNAL FD2

	INTEGER ::ITES,ICHARXI,IU1XI,J,LR,LCR,ICB,ISB,NA1,NA4
	
	REAL(KIND=KIND(1.0D0)) :: SFUP,SFDP,SVACP,SD2P,STOTP,SFDN,SFUN,SVACN,SD2N,STOTN
	REAL(KIND=KIND(1.0D0)) :: SDS,SF,SFS,COMP,SVAC,EI,COMPP,COMPN
        REAL(KIND=KIND(1.0D0)) :: UI,UF,DMT,AMI,AMF,DNT,AGAP,SD2,FD2
        REAL(KIND=8):: EPS,DEK,DMI,DMF,FMIN,EG,GAPBA,EGBA,ACOPtR,UINTERVALO
        REAL(KIND=8):: C1,C2,C3,C4,D2,AU
	INTEGER :: NT,I,K,LFINAL,MT,INORM
	INTEGER :: IS,NS,IJ,IK
        character(len=25) :: correlation,chainsize,date
        character(len=100) :: filename3,filename4	
	COMMON/BA/GAPBA,EGBA
	COMMON/HOPR/ACOPtR
	COMMON/INTERV/UINTERVALO
	COMMON/GAP/AGAP
	EXTERNAL FMU
	COMMON/N/LFINAL	
	COMMON/NORMALIZACAO/INORM
        COMMON/arquivo/correlation,chainsize,date


!*
!*################################################################
!			PARAMETROS
	LFIN=5
!	LFINAL=2
	IZI=1
	IDOWN=1  ! IDOWN=0 (OFF); 1 (ON)
	ACOPtI=1.0D0
	ACOPtR=1.0D0
        D=0.5D0
        ACOPH=0.0D0
	ACOPVI=0.0D0   !first-neighbor electron-electron interaction
	ACPI=0.0D0      !delta ionic Hubbard model
	
        !hopping
	ACt(1)=1.D0
	ACt(2)=1.D0
	ACt(3)=1.D0
	ACt(4)=1.D0
	ACt(5)=1.D0
	ACt(6)=1.D0
	ACt(7)=1.D0
	ACt(8)=1.D0
	ACt(9)=1.D0
	ACt(10)=1.D0
	ACt(LFIN)=0.D0

        EPS=1.D-8
	AMUI=0.0D0
	BETAI=1.D4
	
	C1=2.445D0
        C2=2.581D0
        C3=0.090D0
        C4=0.22D0
!*
!*#################################################################
!	OPEN (44,FILE='saida.dat')
!	OPEN (45,FILE='Arbore.dat')

        WRITE(chainsize,'(I5.1)') LFIN
        call date_and_time(DATE=date)

        filename3=""//trim(adjustl(date))//"_aHartree_LFIN="//&
        trim(adjustl(chainsize))//".dat"
        OPEN (14, FILE=filename3)
        filename4=""//trim(adjustl(date))//"_aocup_LFIN="//&
        trim(adjustl(chainsize))//".dat"
        OPEN (15, FILE=filename4)

!        OPEN (14,FILE='aHartre.dat')
!        OPEN (15,FILE='aocupacao.dat')


	PI=4.0*DATAN(1.D0)


	AMI=0.0D0
	AMF=5.0D0
	MT=1
	DNT=(AMF-AMI)/(MT-1)
	IF(MT.EQ.1)DNT=0.D0
	DO IJ=1,MT
        AMUI=AMI+(IJ-1)*DNT

        UI=4.0D0
        UF=16.0D0
	NT=1
	DMT=(UF-UI)/(NT-1)
	IF(NT.EQ.1)DMT=0.D0
	DO I=1,NT
        ACOPUI=UI+(I-1)*DMT
        AU=ACOPUI
        UINTERVALO=ACOPUI
        ACOPnI=-0.5D0*ACOPUI-AMUI
        CALL EGGAPBA(GAPBA,EGBA)
        
!        WRITE(chainsize,'(I5.1)') LFINAL
        WRITE(correlation,'(10F15.2)') UINTERVALO
!        call date_and_time(DATE=date)
        
!        filename4=""//trim(adjustl(date))//"_aground_agap_LFINAL="//&
!        trim(adjustl(chainsize))//"_ACOPU="//trim(adjustl(correlation))//".dat"
!        OPEN (68, FILE=filename4)
!        filename3=""//trim(adjustl(date))//"_aocup_LFINAL="//&
!        trim(adjustl(chainsize))//"_ACOPU="//trim(adjustl(correlation))//".dat"
!        OPEN (59,FILE=filename3)
        
        DMI=-1.0D0
        DMF=0.0D0
!        DEK=Fmin(DMI,DMF,FD2,EPS)
!        ACOPtI=DEK
        WRITE(6,*)'RESULTADO FINAL RECALCULADO POR GROUND'
!        CALL GROUND(EG,SFUP,SFDP,SVACP,SD2P,SFDN,SFUN,SVACN,SD2N)

!        CALL DENSI2(SFUP,SFDP,SVACP,SD2P,SFDN,SFUN,SVACN,SD2N)


        CALL OCUPT(SFUP,SFDP,SVACP,SD2P,STOTP,SFDN,SFUN,SVACN,SD2N,STOTN)
          COMPP=SVACP+SFUP+SFDP+SD2P
          COMPN=SVACN+SFUN+SFDN+SD2N
          D2=0.25D0*(1.D0+C1*AU)/(1.D0+C2*AU+C3*AU*AU+C4*AU*AU*AU)


!          EG=4.D0*ACOPtI/PI+0.5D0*ACOPUI*(SD2P+SD2N)
!	WRITE(68,121)ACOPUI,ACOPt,ACOPtR,EG,EGBA,AGAP,GAPBA
	WRITE(14,121)ACOPUI,AGAP,GAPBA,EG,EGBA,SVACP,SFUP,SFDP,SD2P,D2,COMPP!,SVACN,SFUN,SFDN,SD2N,COMPN

!        WRITE(15,121)ACOPUI,1.D0/BETA,ACOPtR,AMUI,SVACP+SVACN,SFUP+SFDP,&
!        SFUN+SFDN,COMPP+COMPN,SFUP-SFDP,SFUN-SFDN,SVACP,SFUP,SFDP,SD2P,&
!        D2,COMPP,SVACN,SFUN,SFDN,SD2N,COMPN
!        WRITE(15,121)AMUI,SVACP,SFUP,SFDP,SD2P,STOTP,COMPP,SVACN,SFUN,SFDN,SD2N,STOTN,COMPN, &
        WRITE(15,121)AMUI,0.5D0*(SVACP+SVACN),0.5D0*(SFUP+SFUN),0.5D0*(SFDP+SFDN),0.5D0*(SD2P+SD2N), &
        0.5D0*(COMPP+COMPN),0.5D0*(SFUP+SFUN+SFDP+SFDN+2.D0*(SD2P+SD2N))
121   FORMAT(30F9.5)
        
        ENDDO
	ENDDO

        
	PRINT*,'JEJE'
500	END PROGRAM AMHUBBARD
!*
!**********##############################################**********
!**********#############   R O T I N A S  ###############**********
!**********##############################################**********
	INCLUDE 'EQRT2SA.f'
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          REAL(KIND=8) FUNCTION FD2(ATH)
          USE PRINCIPAL
          IMPLICIT NONE
          REAL(KIND=8):: C1,C2,C3,C4,D2
          REAL(KIND=8):: ATH,U
          REAL(KIND=8):: SFUP,SFDP,SVACP,SD2P,STOTP,SFUN,SFDN,SVACN,SD2N,STOTN
          ACOPtI=ATH
          U=ACOPUI
          C1=2.445D0
          C2=2.581D0
          C3=0.090D0
          C4=0.22D0
          D2=0.25D0*(1.D0+C1*U)/(1.D0+C2*U+C3*U*U+C4*U*U*U)
          CALL OCUPT(SFUP,SFDP,SVACP,SD2P,STOTP,SFDN,SFUN,SVACN,SD2N,STOTN)
          FD2=DABS(D2-SD2P)
!          WRITE(6,*)'DOUBLE=',TH,SD2P,D2
          RETURN
          END

        SUBROUTINE EGGAPBA(GAPBA,EGBA)
        USE PRINCIPAL
        implicit none
        REAL(KIND=8) FF1(3001),FF2(3001)
        real(kind=8):: U,TH,AY,BY,R,H1,GAPBA,EGBA,AX,F,G,RES1,RES2
        INTEGER:: N1,IX,IZ
        external f,g
!        open(22,file='gap_BA.dat')
!        open(13,file='e0_BA.dat')
!        open(123,file='integrando.dat')
        U=ACOPUI
        TH=ACOPtI
        N1=3001
       AY=0.001D0
       BY=30.D0
       R=BY-AY
      H1=R/(N1-1)
        DO 1 IX=1,N1
          AX=AY+(IX-1)*H1
          FF1(IX)=F(AX)
          IF(AX.LT.1.D0) FF1(IX)=0.d0
          FF2(IX)=G(AX)
!        write(123,*)AX,FF1(IX),FF2(IX)
10      format(20F15.8)
1         CONTINUE

          CALL SIMP(N1,H1,FF1,RES1)
          CALL SIMP(N1,H1,FF2,RES2)
          GAPBA=RES1
          EGBA=RES2
          write(22,10)ACOPUI,GAPBA
          write(13,10)ACOPUI,EGBA
        RETURN
        end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Function f(x)
        USE PRINCIPAL
        implicit none
        REAL(KIND=8) f, x
        f = (16.D0*(ACOPtI)*(ACOPtI)/ACOPUI)*(SQRT((x*x)-1.D0))/(DSINH(2.D0*pi*(ACOPtI)*x/ACOPUI))
        return

        end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Function g(x)
        USE PRINCIPAL
        implicit none
        REAL(KIND=8) g, x
        g =-4.D0*(bessel_j0(x)*bessel_j1(x))/(x*(1.d0+dexp(x*ACOPUI/2.D0)))
        return

        end




          REAL FUNCTION FMU(DEQ)
          USE PRINCIPAL
          IMPLICIT NONE
          REAL(KIND=8):: EGBA,GAPBA,EG,DEQ,ACOPtR
          REAL(KIND=8):: SFUP,SFDP,SVACP,SD2P,STOTP,SFDN,SFUN,SVACN,SD2N,STOTN
          COMMON/BA/GAPBA,EGBA 
          COMMON/HOPR/ACOPtR
          ACOPtR=DEQ
          CALL GROUND(EG,SFUP,SFDP,SVACP,SD2P,SFDN,SFUN,SVACN,SD2N)
          FMU=DABS(EGBA-EG)
          WRITE(6,*)'CCCC=',ACOPtR,EGBA,EG
          RETURN
          END


          SUBROUTINE GROUND(EG,SFUP,SFDP,SVACP,SD2P,SFDN,SFUN,SVACN,SD2N)
          USE PRINCIPAL
	  USE AMODO
	  USE AMODHSIS
	  USE FGREEN
          IMPLICIT NONE
          EXTERNAL FERM
	  REAL(KIND=8):: FERM
          COMPLEX(kind=8):: ZW
          REAL(KIND=8)::FF1(3001),RES10,RES1,SDS,SVAC,SF,SD2,EG
          REAL(KIND=8)::ETA,AY,BY,R,H1,VA,VB,ACOPtR
          REAL(KIND=8):: RES2,RES3,FF2(3001),FF3(3001),RES20,RES30,UINTERVALO
          REAL(KIND=8):: SFUP,SFDP,SVACP,SD2P,STOTP,SFDN,SFUN,SVACN,SD2N,STOTN,COMPP,COMPN
          REAL(KIND=8):: GDIFCT,GDIFCFCT,SA,SD,SV,SD1
          INTEGER:: N1,I,IY,LFINAL
          character(len=25) :: correlation,chainsize,date
          character(len=100) :: filename6
          COMMON/HOPR/ACOPtR
          COMMON/INTERV/UINTERVALO
          COMMON/N/LFINAL
          COMMON/arquivo/correlation,chainsize,date
        
          filename6=""//trim(adjustl(date))//"_aintegrando_LFIN="//&
          trim(adjustl(chainsize))//"_ACOPU="//trim(adjustl(correlation))//".dat"
          OPEN (67,FILE=filename6)


          N1=101
          VA=0.001D0
          VB=ACOPUI
!          VA(2)=-0.5D0
!          VB(2)=-0.10D0
!          VA(3)=-0.10D0
!          VB(3)=0.1D0
!          VA(4)=0.1D0
!          VB(4)=0.5D0
!          VA(5)=0.5D0
!          VB(5)=PI
          RES1=0.D0
          DO 10 I=1,1
          AY=VA
          BY=VB
          R=BY-AY
          H1=R/(N1-1)
          IF(N1 .EQ. 1)H1=0.D0
          DO 1 IY=1,N1
          ACOPUI=AY+(IY-1)*H1
          ACOPnI=-0.5D0*ACOPUI-AMUI
!          CALL DENSI2(SFUP,SFDP,SVACP,SD2P,SFDN,SFUN,SVACN,SD2N)
          CALL OCUPT(SFUP,SFDP,SVACP,SD2P,STOTP,SFDN,SFUN,SVACN,SD2N,STOTN)
          COMPP=SVACP+SFUP+SFDP+SD2P
          COMPN=SVACN+SFUN+SFDN+SD2N
!        WRITE(59,111)ACOPU,1.D0/BETA,ACOPtR,AMUI,SFUP,SFDP,SVACP,SD2P,COMPP,SFUN,SFDN,SVACN,SD2N,COMPN
          FF1(IY)=0.5D0*(SD2P+SD2N)
          WRITE(67,111)ACOPU,FF1(IY),0.5D0*(SVACP+SVACN),0.5D0*(SFUP+SFUN),0.5D0*(SFDP+SFDN), &
          0.5D0*(SD2P+SD2N),0.5D0*(COMPP+COMPN),0.5D0*(SFUP+SFUN+SFDP+SFDN+2.D0*(SD2P+SD2N))
1         CONTINUE
          CALL SIMP(N1,H1,FF1,RES10)
          RES1=RES1+RES10
10        CONTINUE
          EG=(-4.D0*ACOPtR/PI)+RES1
!          WRITE(6,111)ACOPU,ACOPt,EG,4.D0*ACOPtR/PI,RES1
          !WRITE(68,111)ACOPU,ACOPt,EG
111       FORMAT(30F13.8)
!        close(67)
        RETURN
        END

	
        real(KIND=8) function funcDensi(x, y)
        USE PRINCIPAL
        IMPLICIT NONE
        COMPLEX (kind=8)::ZW,ZMI,ZG11,ZG13,ZG31,ZG33,ZGF,ZA1,ZA12,ZAINT,ZGF0
        COMPLEX (kind=8)::ZGFP,ZG22,ZG24,ZG42,ZG44,ZGFN
	real(kind=8):: x,y
        real(kind=8):: sq,EFT,ek,wp,wm,arg1,arg2,anFp,anFm
!        real(kind=8):: AMU,TH,E0,U,BETA
 	real(kind=8)::R2,AL
        real(kind=8)::DELK2,t1,t2,A,B,A2
        real(kind=8):: ALFA
        real(kind=8):: WLP,WLM,AAF,BBF,AAC,BBC,AAW,RSI,RSIG,AISIG,ABSIG
        real(kind=8):: HNIJ,ETTA
        REAL(KIND=8):: tx, ty, txy, ax, ay
        REAL(KIND=8)::RFS,RFS13,RFS3,US
!        COMMON/DADOS/AMU,TH,E0,U,BETA
        COMMON/BOSE/R2,AL
        COMMON/ALF/ALFA
        COMMON/OCUPA/HNIJ
        COMMON/DDDD/AAW
        COMMON/TB/tx, ty, txy, ax, ay
        COMMON/RESI/RFS(12),RFS13(12),RFS3(12),US(12)
        OPEN (39,FILE='akxkyc.dat')
        ETTA=1.D-3
        ZW=DCMPLX(AAW,ETTA)
        TX=1.D0
        TY=1.D0
        t1 = -2.D0 * tx * dcos(x)
        t2 = -2.D0 * ty * dcos(y)
        ek = t1 + t2 - amu
        CALL EXATA
        CALL GKONDO(ZW,ZGF0,ZG11,ZG13,ZG31,ZG33,ZGFP,ZG22,ZG24,ZG42,ZG44,ZGFN)
        ZA1=ZMI+2.D0*ACOPtI*(DCOS(X)+DCOS(Y))
        ZA12=ZA1*ZA1
        A2=-4.D0*ACOPtI*ACOPtI
        ZAINT=1.D0/CDSQRT(ZA12+A2)
        IF(DIMAG(ZAINT).LT.0.D0)ZAINT=DCONJG(ZAINT)
        funcDensi=DIMAG(ZAINT)
        end function funcDensi


      SUBROUTINE OCUPT(SFUP,SFDP,SVACP,SD2P,STOTP,SFDN,SFUN,SVACN,SD2N,STOTN)
      USE PRINCIPAL
      IMPLICIT NONE
      real(kind=8):: AI(4),A(3),B(3),SFUP,SFDP,SVACP,SD2P,STOTP,COMPP,SFUN,SFDN,SVACN,SD2N,STOTN,COMPN
      REAL(KIND=8):: RESUP,RESDP,RESVP,RES2P,RESTP,RESDN,RESUN,RESVN,RES2N,RESTN,U
      integer :: I,INORM
      COMMON /INTERVALO/AI
      COMMON/NORMALIZACAO/INORM
      COMMON/INTERV/U

      CALL DENSI2(SFUP,SFDP,SVACP,SD2P,SFDN,SFUN,SVACN,SD2N)
        SFUP=0.D0
        SFDP=0.D0
        SVACP=0.D0
        SD2P=0.D0
        STOTP=0.D0
        SFUN=0.D0
        SFDN=0.D0
        SVACN=0.D0
        SD2N=0.D0
        STOTN=0.D0

!        A(1)=AI(1)-4.D0
!        B(1)=AI(2)
!        A(2)=AI(3)
!        B(2)=AI(4)+4.D0

        A(1)=AI(1)-15.D0
!        B(1)=AI(2)
!        A(2)=AI(2)
!        B(2)=AI(3)
!        A(3)=AI(3)
        B(1)=AI(4)+15.D0

        PRINT*,"OOOOOOOO=",AI(1),AI(2),AI(3),AI(4)

      DO I=1,1
      CALL OCUP(A(I),B(I),RESUP,RESDP,RESVP,RES2P,RESTP,RESDN,RESUN,RESVN,RES2N,RESTN)
          SFUP=SFUP+RESUP
          SFDP=SFDP+RESDP
          SVACP=SVACP+RESVP
          SD2P=SD2P+RES2P
          STOTP=STOTP+RESTP
          SFDN=SFDN+RESDN
          SFUN=SFUN+RESUN
          SVACN=SVACN+RESVN
          SD2N=SD2N+RES2N
          STOTN=STOTN+RESTN
      ENDDO
      RETURN
      END



      SUBROUTINE OCUP(A,B,SFUP,SFDP,SVACP,SD2P,STOTP,SFDN,SFUN,SVACN,SD2N,STOTN)
      USE PRINCIPAL
      IMPLICIT NONE
      EXTERNAL  funcDensi
      EXTERNAL ZFER
      COMPLEX (kind=8):: ZW,ZGA,ZG11,ZG13,ZG31,ZG33,ZG22,ZG24,ZG42,ZG44,ZGFN,ZGFP,ZGF0,ZMI,ZA1,ZA12,ZAINT,ZGFD,ZF(10001),ZGF2D
      COMPLEX (kind=8):: VZGF3(10001),VZGF(10001),VZGFP(10001),VZGFN(10001),VZGF0(10001),VZGF1(10001),VZGF13(10001),VZGF31(10001)
      COMPLEX (kind=8):: VZGF30(10001),ZG110,ZG130,ZG310,ZG330,VZGF10(10001),VZGF130(10001),VZGF310(10001)
      COMPLEX (kind=8):: VZW(10001),ZFER,VZGF2D(10001),VZGF2(10001),VZGF24(10001),VZGF42(10001),VZGF4(10001),ZGAM(10001)
      real(kind=8):: funcDensi
      real(kind=8):: FF2(10001),RES1,RES2,RES20,RF2D(10001)
      real(kind=8):: CC,OI,OF,ETTA,DMT,AW,RFP(10001),RFN(10001),RFNUP(10001),RFNVP(10001),RFN2P(10001),RFNDP(10001),RF0(10001)
      COMPLEX (kind=8):: R11,R13,R33,R2,R3
      COMPLEX (KIND=8):: ZARGUI,ZGA11,ZGA13,ZGA31,ZGA33,ZGA22,ZGA24,ZGA42,ZGA44,ZGAMA,ZGD,ZINT
      COMPLEX (KIND=8):: ZM11,ZM13,ZM31,ZM33,ZM22,ZM24,ZM42,ZM44,ZMI2,ZTETA,ZM110,ZM130,ZM310,ZM330
      real(kind=8):: RFNDN(10001),RFNVN(10001),RFN2N(10001),RFNUN(10001),RFNTN(10001),RFNTP(10001)
      INTEGER::NT,I,J,N1,IY,IJ,IK,LFINAL,K,IZIFINAL,INORM
      REAL(KIND=8):: FF1(10001),VA(2),VB(2),RESV,REST,AY,BY,R,H1,AKY,A2,RES10,VAW(10001),RESA1,RESA2
      REAL(KIND=8):: RESD,SA1,SA2,SFUP,SFDP,SVACP,SD2P,SFDN,SFUN,SVACN,SD2N,STOTP,STOTN
      REAL(KIND=8):: RESUP,RESVP,RES2P,RESDP,RESDN,RESVN,RES2N,RESUN,COMPP,COMPN,A,B,RESTP,RESTN
      real(kind=8):: TH,E0,U,RF2,ROF3,CPI2,RF3,AGAP,AGAP1,AGAP2,ACOPtR,UINTERVALO
      real(kind=8):: axa, axb, aya, ayb,CPI
      character(len=25) :: chainsize,site,hopping,correlation,localenergy,chemicalpotential,temperature
      character(len=25) :: hour,minutes,seconds,milliseconds
      character(len=100) :: filename,filename2
      integer :: values(8)
!      COMMON/DADOS/AMU,TH,E0,U,BETA
      COMMON/DDDD/AW
      COMMON/HOPR/ACOPtR
      COMMON/INTERV/UINTERVALO
      COMMON/N/LFINAL
      COMMON/INDI/IZI
      COMMON/NORMALIZACAO/INORM


        CC=1.D0/PI
        CPI=1.D0/(2.D0*PI)
        CPI2=CPI*CPI
        axa=-PI
        axb=PI
        aya=-PI
        ayb=PI
        NT=10001
        ETTA=5.0D-3

        DO IJ=1,NT
        VZGF0(IJ)=DCMPLX(0.D0,0.D0)
        VZGF10(IJ)=DCMPLX(0.D0,0.D0)
        VZGF130(IJ)=DCMPLX(0.D0,0.D0)
        VZGF310(IJ)=DCMPLX(0.D0,0.D0)
        VZGF30(IJ)=DCMPLX(0.D0,0.D0)
        VZGFP(IJ)=DCMPLX(0.D0,0.D0)
        VZGFN(IJ)=DCMPLX(0.D0,0.D0)
        VZGF1(IJ)=DCMPLX(0.D0,0.D0)
        VZGF13(IJ)=DCMPLX(0.D0,0.D0)
        VZGF31(IJ)=DCMPLX(0.D0,0.D0)
        VZGF3(IJ)=DCMPLX(0.D0,0.D0)
        VZGF2(IJ)=DCMPLX(0.D0,0.D0)
        VZGF24(IJ)=DCMPLX(0.D0,0.D0)
        VZGF42(IJ)=DCMPLX(0.D0,0.D0)
        VZGF4(IJ)=DCMPLX(0.D0,0.D0)
        ZF(IJ)=DCMPLX(0.D0,0.D0)

        ENDDO

!       DO 546 LFIN=2,LFINAL
            IF(MOD(LFIN,2).EQ.0) THEN
                IZIFINAL = LFIN/2
            ELSE
                IZIFINAL = (LFIN+1)/2
            ENDIF
!            IF(ACt(LFIN).NE.0.D0)IZIFINAL=1

        WRITE(6,*)'NNNNNNNN=',INORM

     DO 786 IZI=1,IZIFINAL
     PRINT*,'IZI--OCUP-->',IZI





          AY=A
          BY=B
          R=BY-AY
          DMT=0.D0
          IF(NT.EQ.1)GOTO 20
          DMT=R/(NT-1)
20        CONTINUE

          DO 241 J=1,NT
          VAW(J)=AY+(J-1)*DMT
          VZW(J)=DCMPLX(VAW(J),ETTA)
          ZF(J)=ZFER(DCMPLX(AY+(J-1)*DMT,ETTA))

!          CALL GKONDO(VZW(J),ZGF0,ZG11,ZG13,ZG31,ZG33,ZGFP,ZG110,ZG130,ZG310,ZG330,ZGFN)
!          VZGF0(J)=VZGF0(J)+ZGF0
!          VZGF10(J)=VZGF10(J)+ZG110
!          VZGF30(J)=VZGF30(J)+ZG330
!          VZGF130(J)=VZGF130(J)+ZG130
!          VZGF310(J)=VZGF310(J)+ZG310

          CALL GKONDO(VZW(J),ZGF0,ZG11,ZG13,ZG31,ZG33,ZGFP,ZG22,ZG24,ZG42,ZG44,ZGFN)
          VZGF1(J)=VZGF1(J)+ZG11
          VZGF3(J)=VZGF3(J)+ZG33
          VZGF13(J)=VZGF13(J)+ZG13
          VZGF31(J)=VZGF31(J)+ZG31
          VZGF2(J)=VZGF2(J)+ZG22
          VZGF24(J)=VZGF24(J)+ZG24
          VZGF42(J)=VZGF42(J)+ZG42
          VZGF4(J)=VZGF4(J)+ZG44

241     CONTINUE
786     CONTINUE
!546     CONTINUE


          DO 956 J=1,NT
          VAW(J)=AY+(J-1)*DMT
          ZF(J)=ZFER(DCMPLX(AY+(J-1)*DMT,ETTA))

          ZM11=VZGF1(J)/INORM
          ZM13=VZGF13(J)/INORM
          ZM31=VZGF31(J)/INORM
          ZM33=VZGF3(J)/INORM

          ZGAMA=ZM11+ZM13+ZM31+ZM33
          ZTETA=ZM11*ZM33-ZM13*ZM31

!         ZMI=1.D0/ZGAMA
!         ZMI2=ZMI*ZMI
!         ZARGUI=CDSQRT(ZMI2-4.D0*ACOPtR*ACOPtR)
!         ZINT=1.D0/ZARGUI
!         IF(DIMAG(ZINT).GT.0.D0)ZINT=DCONJG(ZINT)
!
!         VZGF1(J)=(ZM11*ZINT/ZGAMA)
!         VZGF13(J)=(ZM13*ZINT/ZGAMA)
!         VZGF31(J)=(ZM31*ZINT/ZGAMA)
!         VZGF3(J)=(ZM33*ZINT/ZGAMA)
!         VZGFP(J)=VZGF1(J)+VZGF13(J)+VZGF31(J)+VZGF3(J)

          ZARGUI=(1.D0+(D)*ZGAMA)/(1.D0-(D)*ZGAMA)
          ZINT=(0.5D0/D)*CDLOG(ZARGUI)
          IF(DIMAG(ZINT).GT.0.D0)ZINT=DCONJG(ZINT)

          VZGF1(J)=ZTETA/ZGAMA+(ZM11-ZTETA/ZGAMA)*ZINT/ZGAMA
          VZGF13(J)=-ZTETA/ZGAMA+(ZM13+ZTETA/ZGAMA)*ZINT/ZGAMA
          VZGF31(J)=-ZTETA/ZGAMA+(ZM31+ZTETA/ZGAMA)*ZINT/ZGAMA
          VZGF3(J)=ZTETA/ZGAMA+(ZM33-ZTETA/ZGAMA)*ZINT/ZGAMA
          VZGFP(J)=VZGF1(J)+VZGF13(J)+VZGF31(J)+VZGF3(J)

 !         RFNUP(J)=-CC*DIMAG((VZGFP(J)-VZGF3(J))*ZF(J))
 !         RFNVP(J)=-CC*DIMAG((VZGFP(J)-VZGF3(J))*(1.D0-ZF(J)))
          RFN2P(J)=-CC*DIMAG(VZGF3(J)*ZF(J))
          RFNDP(J)=-CC*DIMAG(VZGF3(J)*(1.D0-ZF(J)))

          RFNUP(J)=-CC*DIMAG(VZGF1(J)*ZF(J))
          RFNVP(J)=-CC*DIMAG(VZGF1(J)*(1.D0-ZF(J)))
!          RFN2P(J)=-CC*DIMAG((VZGFP(J)-VZGF1(J))*ZF(J))
!          RFNDP(J)=-CC*DIMAG((VZGFP(J)-VZGF1(J))*(1.D0-ZF(J)))
          RFNTP(J)=-CC*DIMAG(VZGFP(J)*ZF(J))


          ZM22=VZGF2(J)/INORM
          ZM24=VZGF24(J)/INORM
          ZM42=VZGF42(J)/INORM
          ZM44=VZGF4(J)/INORM

          ZGAMA=ZM22+ZM24+ZM42+ZM44
          ZTETA=ZM22*ZM44-ZM24*ZM42

!         ZMI=1.D0/ZGAMA
!         ZMI2=ZMI*ZMI
!         ZARGUI=CDSQRT(ZMI2-4.D0*ACOPtR*ACOPtR)
!         ZINT=1.D0/ZARGUI
!         IF(DIMAG(ZINT).GT.0.D0)ZINT=DCONJG(ZINT)

!         VZGF2(J)=(ZM22*ZINT/ZGAMA)
!         VZGF24(J)=(ZM24*ZINT/ZGAMA)
!         VZGF42(J)=(ZM42*ZINT/ZGAMA)
!         VZGF4(J)=(ZM44*ZINT/ZGAMA)
!         VZGFN(J)=VZGF2(J)+VZGF24(J)+VZGF42(J)+VZGF4(J)

          ZARGUI=(1.D0+(D)*ZGAMA)/(1.D0-(D)*ZGAMA)
          ZINT=(0.5D0/D)*CDLOG(ZARGUI)
          IF(DIMAG(ZINT).GT.0.D0)ZINT=DCONJG(ZINT)

          VZGF2(J)=ZTETA/ZGAMA+(ZM22-ZTETA/ZGAMA)*ZINT/ZGAMA
          VZGF24(J)=-ZTETA/ZGAMA+(ZM24+ZTETA/ZGAMA)*ZINT/ZGAMA
          VZGF42(J)=-ZTETA/ZGAMA+(ZM42+ZTETA/ZGAMA)*ZINT/ZGAMA
          VZGF4(J)=ZTETA/ZGAMA+(ZM44-ZTETA/ZGAMA)*ZINT/ZGAMA
          VZGFN(J)=VZGF2(J)+VZGF24(J)+VZGF42(J)+VZGF4(J)

 !         RFNUN(J)=-CC*DIMAG((VZGFN(J)-VZGF2(J))*(1.D0-ZF(J)))
 !         RFNVN(J)=-CC*DIMAG((VZGFN(J)-VZGF4(J))*(1.D0-ZF(J)))
          RFN2N(J)=-CC*DIMAG(VZGF4(J)*ZF(J))
          RFNUN(J)=-CC*DIMAG(VZGF4(J)*(1.D0-ZF(J)))

          RFNDN(J)=-CC*DIMAG(VZGF2(J)*ZF(J))
          RFNVN(J)=-CC*DIMAG(VZGF2(J)*(1.D0-ZF(J)))
!          RFN2N(J)=-CC*DIMAG((VZGFN(J)-VZGF2(J))*ZF(J))
!          RFNDN(J)=-CC*DIMAG((VZGFN(J)-VZGF4(J))*ZF(J))
          RFNTN(J)=-CC*DIMAG(VZGFN(J)*ZF(J))

956     CONTINUE

          CALL SIMP(NT,DMT,RFNUP,RESUP)
          CALL SIMP(NT,DMT,RFNVP,RESVP)
          CALL SIMP(NT,DMT,RFN2P,RES2P)
          CALL SIMP(NT,DMT,RFNDP,RESDP)
          CALL SIMP(NT,DMT,RFNTP,RESTP)

          CALL SIMP(NT,DMT,RFNDN,RESDN)
          CALL SIMP(NT,DMT,RFNVN,RESVN)
          CALL SIMP(NT,DMT,RFN2N,RES2N)
          CALL SIMP(NT,DMT,RFNUN,RESUN)
          CALL SIMP(NT,DMT,RFNTN,RESTN)

          SFUP=RESUP
          SFDP=RESDP
          SVACP=RESVP
          SD2P=RES2P
          STOTP=RESTP

          SFDN=RESDN
          SFUN=RESUN
          SVACN=RESVN
          SD2N=RES2N
          STOTN=RESTN


      RETURN
      END



      SUBROUTINE DENSI2(SFUP,SFDP,SVACP,SD2P,SFDN,SFUN,SVACN,SD2N)
      USE PRINCIPAL
      IMPLICIT NONE
      EXTERNAL  funcDensi
      EXTERNAL ZFER
      COMPLEX (kind=8):: ZW,ZGA,ZG11,ZG13,ZG31,ZG33,ZG22,ZG24,ZG42,ZG44,ZGFN,ZGFP,ZGF0,ZMI,ZA1,ZA12,ZAINT,ZGFD,ZF(20001),ZGF2D
      COMPLEX (kind=8):: VZGF3(20001),VZGFP(20001),VZGFN(20001),VZGF0(20001),VZGF1(20001),VZGF13(20001),VZGF31(20001)
      COMPLEX (kind=8):: ZG330,ZG110,ZG130,ZG310,ZGAMA0
      COMPLEX (kind=8):: VZGF10(20001),VZGF130(20001),VZGF310(20001),VZGF30(20001)            
      COMPLEX (kind=8):: VZW(20001),ZFER,VZGF2D(20001),VZGF2(20001),VZGF24(20001),VZGF42(20001),VZGF4(20001)
      real(kind=8):: funcDensi
      real(kind=8):: FF2(20001),RES1,RES2,RES20,RF2D(20001),AI(4)
      real(kind=8):: CC,OI,OF,ETTA,DMT,AW,RFP(20001),RFN(20001),RFNUP(20001),RFNVP(20001),RFN2P(20001),RFNDP(20001),RF0(20001)
      COMPLEX (kind=8):: R11,R13,R33,R2,R3
      COMPLEX (KIND=8):: ZARGUI,ZGA11,ZGA13,ZGA31,ZGA33,ZGA22,ZGA24,ZGA42,ZGA44,ZGAMA,ZGD,ZINT
      COMPLEX (KIND=8):: ZM11,ZM13,ZM31,ZM33,ZM22,ZM24,ZM42,ZM44,ZMI2,ZTETA
      real(kind=8):: RFNDN(20001),RFNVN(20001),RFN2N(20001),RFNUN(20001)
      INTEGER::NT,I,J,N1,IY,IJ,IK,LFINAL,K,IZIFINAL,INORM,JJ
      REAL(KIND=8):: FF1(20001),VA(5),VB(5),RESV,REST,AY,BY,R,H1,AKY,A2,RES10,VAW(20001),RESA1,RESA2
      REAL(KIND=8):: RESD,SA1,SA2,SFUP,SFDP,SVACP,SD2P,SFDN,SFUN,SVACN,SD2N
      REAL(KIND=8):: RESUP,RESVP,RES2P,RESDP,RESDN,RESVN,RES2N,RESUN,COMPP,COMPN
      real(kind=8):: TH,E0,U,RF2,ROF3,CPI2,RF3,AGAP,AGAP1,AGAP2,AGAP3,AGAP4,ACOPtR,GAPBA,EGBA,UINTERVALO
      real(kind=8):: axa, axb, aya, ayb,CPI
      character(len=25) :: chainsize,site,hopping,correlation,localenergy,chemicalpotential,temperature
      character(len=25) :: hour,minutes,seconds,milliseconds
      character(len=100) :: filename,filename2
      integer :: values(8),IL
!      COMMON/DADOS/AMU,TH,E0,U,BETA
      COMMON/DDDD/AW
      COMMON/BA/GAPBA,EGBA
      COMMON/HOPR/ACOPtR
      COMMON/INTERV/UINTERVALO
      COMMON /INTERVALO/AI
      COMMON/GAP/AGAP
      COMMON/N/LFINAL
      COMMON/NORMALIZACAO/INORM
      character(len=25) :: date
      COMMON/arquivo/correlation,chainsize,date

        INORM=0
        CC=1.D0/PI
        CPI=1.D0/(2.D0*PI)
        CPI2=CPI*CPI
        axa=-PI
        axb=PI
        aya=-PI
        ayb=PI
        NT=20001
        ETTA=1.D-8

        DO IJ=1,NT
        VZGF0(IJ)=DCMPLX(0.D0,0.D0)
        VZGF10(IJ)=DCMPLX(0.D0,0.D0)
        VZGF130(IJ)=DCMPLX(0.D0,0.D0)
        VZGF310(IJ)=DCMPLX(0.D0,0.D0)
        VZGF30(IJ)=DCMPLX(0.D0,0.D0)
        VZGFP(IJ)=DCMPLX(0.D0,0.D0)
        VZGFN(IJ)=DCMPLX(0.D0,0.D0)
        VZGF1(IJ)=DCMPLX(0.D0,0.D0)
        VZGF13(IJ)=DCMPLX(0.D0,0.D0)
        VZGF31(IJ)=DCMPLX(0.D0,0.D0)
        VZGF3(IJ)=DCMPLX(0.D0,0.D0)
        VZGF2(IJ)=DCMPLX(0.D0,0.D0)
        VZGF24(IJ)=DCMPLX(0.D0,0.D0)
        VZGF42(IJ)=DCMPLX(0.D0,0.D0)
        VZGF4(IJ)=DCMPLX(0.D0,0.D0)
        ZF(IJ)=DCMPLX(0.D0,0.D0)
        VZGF2D(IJ)=DCMPLX(0.D0,0.D0)
        ENDDO

!       DO 546 LFIN=2,LFINAL
            IF(MOD(LFIN,2).EQ.0) THEN
                IZIFINAL = LFIN/2
            ELSE
                IZIFINAL = (LFIN+1)/2
            ENDIF
!            IF(ACt(LFIN).NE.0.D0)IZIFINAL=1
        INORM=INORM+IZIFINAL
        WRITE(6,*)'NNNNNNNN=',INORM
        CALL RESIDUOS
     DO 786 IZI=1,IZIFINAL
     PRINT*,'IZI--DENSI2-->',IZI
!        CALL RESIDUOS

        CALL MEOCUP
        IF(IDOWN.EQ.1)CALL MEOCDOWN

        OI=-15.0D0-UINTERVALO
        OF=15.0D0+UINTERVALO
        R=OF-OI
        DMT=0.D0
        IF(NT.EQ.1)GOTO 20
        DMT=R/(NT-1)
20      CONTINUE

        DO 241 J=1,NT
          VAW(J)=OI+(J-1)*DMT
          VZW(J)=DCMPLX(VAW(J),ETTA)
          ZF(J)=ZFER(DCMPLX(OI+(J-1)*DMT,ETTA))

!          CALL GKONDO(VZW(J),ZGF0,ZG11,ZG13,ZG31,ZG33,ZGFP,ZG110,ZG130,ZG310,ZG330,ZGFN)
!          VZGF0(J)=VZGF0(J)+ZGF0
!          VZGF10(J)=VZGF10(J)+ZG110
!          VZGF30(J)=VZGF30(J)+ZG330
!          VZGF130(J)=VZGF130(J)+ZG130
!          VZGF310(J)=VZGF310(J)+ZG310

          CALL GKONDO(VZW(J),ZGF0,ZG11,ZG13,ZG31,ZG33,ZGFP,ZG22,ZG24,ZG42,ZG44,ZGFN)
          VZGF1(J)=VZGF1(J)+ZG11
          VZGF3(J)=VZGF3(J)+ZG33
          VZGF13(J)=VZGF13(J)+ZG13
          VZGF31(J)=VZGF31(J)+ZG31
          VZGF2(J)=VZGF2(J)+ZG22
          VZGF24(J)=VZGF24(J)+ZG24
          VZGF42(J)=VZGF42(J)+ZG42
          VZGF4(J)=VZGF4(J)+ZG44

!         i=10
!         CALL trapzd2d(axa, axb, aya, ayb, ROF3, funcDensi,  i)
!         RF3=CC*CPI2*ROF3
!         WRITE(6,*)'PASSEI POR AQUI',J,RF3
!         VA(1)=-PI
!         VB(1)=-0.50D0
!         VA(2)=-0.5D0
!         VB(2)=-0.10D0
!         VA(3)=-0.10D0
!         VB(3)=0.1D0
!         VA(4)=0.1D0
!         VB(4)=0.5D0
!         VA(5)=0.5D0
!         VB(1)=PI
!          DO 10 I=1,1
!          AY=VA(I)
!          BY=VB(I)
!          R=BY-AY
!          H1=R/(N1-1)
!          DO 1 IY=1,N1
!          AKY=AY+(IY-1)*H1
!          ZA1=ZMI+2.D0*ACOPtI*DCOS(AKY)
!          ZA12=ZA1*ZA1
!          A2=-4.D0*ACOPtI*ACOPtI
!          ZAINT=1.D0/CDSQRT(ZA12+A2)
!          IF(DIMAG(ZAINT).LT.0.D0)ZAINT=DCONJG(ZAINT)
!          FF1(IY)=DREAL(ZAINT)
!          FF2(IY)=DIMAG(ZAINT)
!          WRITE(12,121)AKY,FF1(IY)
!1         CONTINUE
!          CALL SIMP(N1,H1,FF1,RES10)
!          CALL SIMP(N1,H1,FF2,RES20)
!          RES1=RES1+RES10
!          RES2=RES2+RES20
!10        CONTINUE
!          RES1=CPI*RES1
!            RES2=CPI*RES2
!            ZGF2D=DCMPLX(RES1,RES2)
!            VZGF2D(J)=VZGF2D(J)+ZGF2D
!          WRITE(6,121)AW,RF1,RF2,RF3

241     CONTINUE
786     CONTINUE
!546     CONTINUE


          DO 956 J=1,NT
          VAW(J)=OI+(J-1)*DMT
          ZF(J)=ZFER(DCMPLX(OI+(J-1)*DMT,ETTA))

          ZM11=VZGF1(J)/INORM
          ZM13=VZGF13(J)/INORM
          ZM31=VZGF31(J)/INORM
          ZM33=VZGF3(J)/INORM

          ZGAMA=ZM11+ZM13+ZM31+ZM33
          ZTETA=ZM11*ZM33-ZM13*ZM31

!         ZMI=1.D0/ZGAMA
!         ZMI2=ZMI*ZMI
!         ZARGUI=CDSQRT(ZMI2-4.D0*ACOPtR*ACOPtR)
!         ZINT=1.D0/ZARGUI
!         IF(DIMAG(ZINT).GT.0.D0)ZINT=DCONJG(ZINT)
!
!         VZGF1(J)=(ZM11*ZINT/ZGAMA)
!         VZGF13(J)=(ZM13*ZINT/ZGAMA)
!         VZGF31(J)=(ZM31*ZINT/ZGAMA)
!         VZGF3(J)=(ZM33*ZINT/ZGAMA)
!         VZGFP(J)=VZGF1(J)+VZGF13(J)+VZGF31(J)+VZGF3(J)

          ZARGUI=(1.D0+(D)*ZGAMA)/(1.D0-(D)*ZGAMA)
          ZINT=(0.5D0/D)*CDLOG(ZARGUI)
          IF(DIMAG(ZINT).GT.0.D0)ZINT=DCONJG(ZINT)

          VZGF1(J)=ZTETA/ZGAMA+(ZM11-ZTETA/ZGAMA)*ZINT/ZGAMA
          VZGF13(J)=-ZTETA/ZGAMA+(ZM13+ZTETA/ZGAMA)*ZINT/ZGAMA
          VZGF31(J)=-ZTETA/ZGAMA+(ZM31+ZTETA/ZGAMA)*ZINT/ZGAMA
          VZGF3(J)=ZTETA/ZGAMA+(ZM33-ZTETA/ZGAMA)*ZINT/ZGAMA
          VZGFP(J)=VZGF1(J)+VZGF13(J)+VZGF31(J)+VZGF3(J)

          RFP(J)=-CC*DIMAG(VZGFP(J))


          ZM22=VZGF2(J)/INORM
          ZM24=VZGF24(J)/INORM
          ZM42=VZGF42(J)/INORM
          ZM44=VZGF4(J)/INORM

          ZGAMA=ZM22+ZM24+ZM42+ZM44
          ZTETA=ZM22*ZM44-ZM24*ZM42

!         ZMI=1.D0/ZGAMA
!         ZMI2=ZMI*ZMI
!         ZARGUI=CDSQRT(ZMI2-4.D0*ACOPtR*ACOPtR)
!         ZINT=1.D0/ZARGUI
!         IF(DIMAG(ZINT).GT.0.D0)ZINT=DCONJG(ZINT)

!         VZGF2(J)=(ZM22*ZINT/ZGAMA)
!         VZGF24(J)=(ZM24*ZINT/ZGAMA)
!         VZGF42(J)=(ZM42*ZINT/ZGAMA)
!         VZGF4(J)=(ZM44*ZINT/ZGAMA)
!         VZGFN(J)=VZGF2(J)+VZGF24(J)+VZGF42(J)+VZGF4(J)

          ZARGUI=(1.D0+(D)*ZGAMA)/(1.D0-(D)*ZGAMA)
          ZINT=(0.5D0/D)*CDLOG(ZARGUI)
          IF(DIMAG(ZINT).GT.0.D0)ZINT=DCONJG(ZINT)

          VZGF2(J)=ZTETA/ZGAMA+(ZM22-ZTETA/ZGAMA)*ZINT/ZGAMA
          VZGF24(J)=-ZTETA/ZGAMA+(ZM24+ZTETA/ZGAMA)*ZINT/ZGAMA
          VZGF42(J)=-ZTETA/ZGAMA+(ZM42+ZTETA/ZGAMA)*ZINT/ZGAMA
          VZGF4(J)=ZTETA/ZGAMA+(ZM44-ZTETA/ZGAMA)*ZINT/ZGAMA
          VZGFN(J)=VZGF2(J)+VZGF24(J)+VZGF42(J)+VZGF4(J)

          RFN(J)=-CC*DIMAG(VZGFN(J))


!          IMPLEMENTAR CÁLCULO DAS FG DE ORDEM ZERO AQUI
!          RF0(J)=-CC*DIMAG(VZGF0(J))


!          IF(ABS(ACOPU-UINTERVALO).LT.10.D-8.AND.UINTERVALO.NE.0.D0) THEN
!          filename=""//trim(adjustl(date))//"_adensity_LFIN="//&
!          trim(adjustl(chainsize))//"_ACOPU="//trim(adjustl(correlation))//".dat"
!          OPEN (12,FILE=filename)
!          WRITE(12,121)VAW(J),RFP(J),RFN(J)!,RF0(J)
!          END IF

121   FORMAT(30F12.6)

956     CONTINUE

          JJ=0
          DO IK=1,(NT+1)/2
          IF(RFP(IK).GT.10.D3*ETTA)THEN
          AGAP1=VAW(IK)
          JJ=JJ+1
          AI(JJ)=AGAP1
          GO TO 778
          END IF
          ENDDO
778       CONTINUE

          DO IK=(NT+1)/2,1,-1
          IF(RFP(IK).GT.10.D3*ETTA)THEN
          AGAP2=VAW(IK)
          JJ=JJ+1
          AI(JJ)=AGAP2
          GO TO 790
          END IF
          ENDDO
790       CONTINUE

          DO IK=(NT+1)/2,NT
          IF(RFP(IK).GT.10.D3*ETTA)THEN
          AGAP3=VAW(IK)
          JJ=JJ+1
          AI(JJ)=AGAP3
          GO TO 779
          END IF
          ENDDO
779       CONTINUE

          AGAP=AGAP3-AGAP2

          DO IK=NT,1,-1
          IF(RFP(IK).GT.10.D3*ETTA)THEN
          AGAP4=VAW(IK)
          JJ=JJ+1
          AI(JJ)=AGAP4
          GO TO 791
          END IF
          ENDDO
791       CONTINUE
          DO IL=1,4
          PRINT*,'AAAAAAAAAA=',AGAP,AI(IL)
          ENDDO


      RETURN
      END




        subroutine trapzd2d(a, b, ya, yb, stxy, func, i)

        INTEGER i
       real(kind=8):: a,b

        real(kind=8):: ya, yb

        real(kind=8):: func
        external func

        integer j,k

       real(kind=8):: delx, dely

        integer it
        real(kind=8):: tnm

        real(kind=8):: x,y
        real(kind=8):: stxy

        it = 2**(i-2)

        tnm = 1.* it

        if (i.eq.1) then
           stxy =0.5*(b-a)*( yb - ya )*(func( b, yb )+func(a,ya ))

         else

            delx = ( b - a ) / tnm
            x = a + ( 0.5 * delx )

            stxy=0.
            do 11 j = 1 , it

               dely = ( yb - ya ) / tnm
               y = ya + ( 0.5 * dely )

               do 21 k = 1 , it
                  stxy = stxy +  delx * dely * func( x, y )
                  y = y + dely
    21         enddo

               x = x + delx
    11      enddo

        endif
        end subroutine trapzd2d



	
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE GKONDO(ZWD,ZGF0,ZM11,ZM13,ZM31,ZM33,ZGFP,ZM22,ZM24,ZM42,ZM44,ZGFN)
	USE PRINCIPAL
	USE AMODO
	USE AMODHSIS
	USE FGREEN
	IMPLICIT NONE
        
	COMPLEX(KIND=KIND(1.0D0)) :: ZWD,ZGFAT11,ZGFAT13,ZGFAT33,ZM011,ZM033,ZARG33,ZGFATD11,ZGFATD13,ZGFATD33,ZDENO
	
	COMPLEX(KIND=KIND(1.0D0)) ::ZGFAT44,ZGFAT22,ZGFAT24

	COMPLEX(KIND=KIND(1.0D0)) :: ZGA11,ZGA13,ZGA31,ZGA33,ZARG13,ZARG31,ZM33,ZM11,ZM13,ZM31,ZM0,ACOPn0
	
	COMPLEX(KIND=KIND(1.0D0)) :: ZGAD11,ZGAD13,ZGAD31,ZGAD33,ZGAD,Z1AD,Z1BD,Z1CD,Z1DD,ZDETAD
	
	COMPLEX(KIND=KIND(1.0D0)) :: ZMD11,ZMD13,ZMD31,ZMD33,ZGAMAD,ZGAMATOTAL,ZM11TOTAL,ZM13TOTAL,ZM31TOTAL,ZM33TOTAL

	COMPLEX(KIND=KIND(1.0D0)) :: ZDETA,ZA,ZB,ZC,ZD,Z1A,Z1B,Z1C,Z1D,ZARG,ZM2,ZM330,ZM110,ZM130,ZM310,ZGA110,ZGA130,ZGA310,ZGA330

	COMPLEX(KIND=KIND(1.0D0)) :: ZDETB,ZGF,ZG11,ZG33,ZG13,ZG31,ZGA,ZARG1,ZARG2,ZARGUI,ZGD11,ZGD13,ZGD31,ZGD33,ZGFD
	
	COMPLEX(KIND=KIND(1.0D0)) :: ZG22,ZG24,ZG42,ZG44,ZGFN,ZGFP,ZM22,ZM24,ZM42,ZM44,ZGA22,ZGA24,ZGA42,ZGA44
	
	COMPLEX(KIND=KIND(1.0D0)) :: ZARG0,ZGF0,ZARG033,ZG330,ZGAMA,ZGD,ZINT,ZTETA,ZGAMA0,ZINT0,ZARGUI0,ZARGUID,ZGDD,ZINTD,ZMID,ZMID2
	
	COMPLEX(KIND=KIND(1.0D0)) :: ZARG011,ZG110,ZARG11,ZMI,ZMI2,ZMI0,ZMI20,ZG130,ZG310,ZG220,ZG240,ZG420,ZG440
	
        REAL(KIND=KIND(1.0D0)) :: OI,OF,DMT,AKI,AKF,DK,EK,ACOPtR
        
	REAL(KIND=KIND(1.0D0)) :: TTIL,ETTA,AW,CC,DEF,RENOR,RENOR33
	
	REAL(KIND=KIND(1.0D0)) :: RFS(12),RFS13(12),RFS3(12),US(12)
	
!        REAL(KIND=KIND(1.0D0)) :: AMU,TH,TH0,E0,U,BETA
	
	REAL(KIND=KIND(1.0D0)) :: E1(4),DXP(4),AZMON,AM11,AM33,EMINMON
	
	INTEGER :: I,MT,NT,J,I1,J1,NBX,NBY,JBTX,JBTY,G,LFINAL
	
	COMMON/HOPR/ACOPtR
	COMMON/N/LFINAL
	COMMON/INDI/IZI 
 
        
!       VARIÁVEIS REMOVIDAS DAS DECLARAÇÕES POIS ESTÃO NO MÓDULO:
!       D,EMIN

        ACOPn0=ACOPn

	AZMON=0.0D0
	E1(1)=0.D0
	E1(2)=ACOPn0
	E1(3)=ACOPn0
	E1(4)=2.D0*ACOPn0+ACOPU
	EMINMON=E1(1)
	DO 50 I=2,4
	IF(EMINMON.GT.E1(I)) EMINMON=E1(I)
50	CONTINUE
	DO I=1,4
	E1(I)=E1(I)-EMINMON
	ENDDO
	DO 61 I=1,4
	DXP(I) = DEXP(-BETA*E1(I))
	AZMON=AZMON+DXP(I)
61	CONTINUE

	AM11=(DXP(1)+DXP(2))/AZMON
	AM33=(DXP(3)+DXP(4))/AZMON
	ZM011=AM11/(ZWD-ACOPn0)
	ZM033=AM33/(ZWD-ACOPn0-ACOPU)
!	ZGAMA0=ZM011+ZM033
!	ZMI0=1.D0/ZGAMA0
!	ZMI20=ZMI0*ZMI0
!	ZARG0=CDSQRT(ZMI20-4.D0*ACOPt*ACOPt)	
!	IF(DIMAG(ZARG0).LT.0.D0)ZARG0=DCONJG(ZARG0)
!	ZINT0=1.D0/ZARG0
!	ZG110=ZM011*ZINT0/ZGAMA0	
!	ZG330=ZM033*ZINT0/ZGAMA0
!!	IF(DIMAG(ZG110).GT.0.D0)ZG110=DCONJG(ZG110)
!!	IF(DIMAG(ZG330).GT.0.D0)ZG330=DCONJG(ZG330)	
!	ZGF0= ZG110 + ZG330

	ZGA110=ZM011
	ZGA330=ZM033
	ZGA130=DCMPLX(0.D0,0.D0)
	ZGA310=DCMPLX(0.D0,0.D0)
      ZM110=ZGA110
      ZM130=ZGA130
      ZM310=ZGA310
      ZM330=ZGA330
      
		 ZGAMA0=ZM110+ZM130+ZM310+ZM330
         ZMI0=1.D0/ZGAMA0
         ZMI20=ZMI0*ZMI0
!         ZARGUI0=CDSQRT(ZMI20-4.D0*ACOPtR*ACOPtR)
!         ZINT0=1.D0/ZARGUI0
		 ZARGUI0=(1.D0+(D-AMU)*ZGAMA0)/(1.D0-(D+AMU)*ZGAMA0) 
		 ZINT0=(0.5D0/D)*CDLOG(ZARGUI0)        
         IF(DIMAG(ZARGUI0).LT.0.D0)ZARGUI0=DCONJG(ZARGUI0)
         ZG110=ZM110*ZINT0/ZGAMA0
         ZG130=ZM130*ZINT0/ZGAMA0
         ZG310=ZM310*ZINT0/ZGAMA0
         ZG330=ZM330*ZINT0/ZGAMA0
         ZGF0=ZG110+ZG130+ZG310+ZG330
      
!       ZTETA=ZM11*ZM33-ZM13*ZM31
!       ZARGUI=(1.D0+(D-AMU)*ZGAMA)/(1.D0-(D+AMU)*ZGAMA)
!          ZINT=(0.5D0/D)*CDLOG(ZARGUI)
! !         IF(DIMAG(ZINT).GT.0.D0)ZINT=DCONJG(ZINT)
!          ZG110=ZTETA/ZGAMA+(ZM11-ZTETA/ZGAMA)*ZINT/ZGAMA
!          ZG130=-ZTETA/ZGAMA+(ZM13+ZTETA/ZGAMA)*ZINT/ZGAMA
!          ZG310=-ZTETA/ZGAMA+(ZM31+ZTETA/ZGAMA)*ZINT/ZGAMA
!          ZG330=ZTETA/ZGAMA+(ZM33-ZTETA/ZGAMA)*ZINT/ZGAMA
         ZG220=ZG110
         ZG240=ZG130
         ZG420=ZG310
         ZG440=ZG330

!        IF(LFIN.GT.2)THEN
!        ZM011=DCMPLX(0.D0,0.D0)
!        ZM033=DCMPLX(0.D0,0.D0)

!	ZG110=DCMPLX(0.D0,0.D0)
!	ZG130=DCMPLX(0.D0,0.D0)
!	ZG310=DCMPLX(0.D0,0.D0)
!	ZG330=DCMPLX(0.D0,0.D0)
!	ZG220=DCMPLX(0.D0,0.D0)
!	ZG240=DCMPLX(0.D0,0.D0)
!	ZG420=DCMPLX(0.D0,0.D0)
!	ZG440=DCMPLX(0.D0,0.D0)
!	END IF

        ZGF0=ZG110+ZG130+ZG310+ZG330
	
!	CALL EXATA
!        
!	RFS=DM11
!	RFS13=DM13
!	RFS3=DM33
!	US=DUT
!        
!	ZGFATD11=RFS(1)/(ZWD-US(1))+RFS(2)/(ZWD-US(2))+RFS(3)/(ZWD-US(3))+ &
!	RFS(4)/(ZWD-US(4))+RFS(5)/(ZWD-US(5))+RFS(6)/(ZWD-US(6))+ &
!	RFS(7)/(ZWD-US(7))+RFS(8)/(ZWD-US(8))+RFS(9)/(ZWD-US(9))+ &
!	RFS(10)/(ZWD-US(10))+RFS(11)/(ZWD-US(11))+RFS(12)/(ZWD-US(12))
!	ZGFATD13=RFS13(1)/(ZWD-US(1))+RFS13(2)/(ZWD-US(2))+ &
!	RFS13(3)/(ZWD-US(3))+RFS13(4)/(ZWD-US(4))+RFS13(5)/(ZWD-US(5))+ &
!	RFS13(6)/(ZWD-US(6))+RFS13(7)/(ZWD-US(7))+RFS13(8)/(ZWD-US(8))+ &
!	RFS13(9)/(ZWD-US(9))+RFS13(10)/(ZWD-US(10))+RFS13(11)/(ZWD-US(11))+ &
!	RFS13(12)/(ZWD-US(12))
!	ZGFATD33=RFS3(1)/(ZWD-US(1))+RFS3(2)/(ZWD-US(2))+ &
!	RFS3(3)/(ZWD-US(3))+RFS3(4)/(ZWD-US(4))+RFS3(5)/(ZWD-US(5))+ &
!	RFS3(6)/(ZWD-US(6))+RFS3(7)/(ZWD-US(7))+RFS3(8)/(ZWD-US(8))+ &
!	RFS3(9)/(ZWD-US(9))+RFS3(10)/(ZWD-US(10))+ &
!	RFS3(11)/(ZWD-US(11))+RFS3(12)/(ZWD-US(12))
         
	ZGFAT11=DCMPLX(0.D0,0.D0)
	ZGFAT33=DCMPLX(0.D0,0.D0)
	ZGFAT13=DCMPLX(0.D0,0.D0)
	ZGFAT22=DCMPLX(0.D0,0.D0)
	ZGFAT24=DCMPLX(0.D0,0.D0)
	ZGFAT44=DCMPLX(0.D0,0.D0)
	ZDENO=DCMPLX(0.D0,0.D0)
	RENOR=0.0D0

	DO 100 G=1,P
        NBX=IB1CUPSISX(V(G,1))
        NBY=IB1CUPSISY(V(G,1))
       ZDENO=ZWD + WPHOS(NBY)%op3(V(G,3)) - WPHOS(NBX)%op3(V(G,2))
!        ZDENO=ZWD - WPHOS(NBY)%op3(V(G,3)) + WPHOS(NBX)%op3(V(G,2))

!        PRINT*,V(G,1),V(G,2),V(G,3),CG33(V(G,1))%op(V(G,2),V(G,3))

        ZGFAT33=ZGFAT33+CG33(V(G,1))%op(V(G,2),V(G,3))/ZDENO

100	END DO

        DO 200 G=1,PP
        NBX=IB1CUPSISX(VV(G,1))
        NBY=IB1CUPSISY(VV(G,1))
       ZDENO=ZWD + WPHOS(NBY)%op3(VV(G,3)) - WPHOS(NBX)%op3(VV(G,2))
!        ZDENO=ZWD - WPHOS(NBY)%op3(VV(G,3)) + WPHOS(NBX)%op3(VV(G,2))

!        PRINT*,VV(G,1),VV(G,2),VV(G,3),CG11(VV(G,1))%op(VV(G,2),VV(G,3))

        ZGFAT11=ZGFAT11+CG11(VV(G,1))%op(VV(G,2),VV(G,3))/ZDENO

200     END DO

        DO 300 G=1,PPP
        NBX=IB1CUPSISX(VVV(G,1))
        NBY=IB1CUPSISY(VVV(G,1))
       ZDENO=ZWD + WPHOS(NBY)%op3(VVV(G,3)) - WPHOS(NBX)%op3(VVV(G,2))
!        ZDENO=ZWD - WPHOS(NBY)%op3(VVV(G,3)) + WPHOS(NBX)%op3(VVV(G,2))

!        PRINT*,VVV(G,1),VVV(G,2),VVV(G,3),CG13(VVV(G,1))%op(VVV(G,2),VVV(G,3))

        ZGFAT13=ZGFAT13+CG13(VVV(G,1))%op(VVV(G,2),VVV(G,3))/ZDENO

300     END DO          !#BLOQUES-CUPSIS


    IF(IDOWN.EQ.1)THEN
	DO 400 G=1,L
	NBX=IB1CDOWNSISX(S(G,1))
	NBY=IB1CDOWNSISY(S(G,1))
	ZDENO=ZWD + WPHOS(NBY)%op3(S(G,3)) - WPHOS(NBX)%op3(S(G,2))
!	ZDENO=ZWD - WPHOS(NBY)%op3(S(G,3)) + WPHOS(NBX)%op3(S(G,2))
	
!        PRINT*,S(G,1),S(G,2),S(G,3),CG44(S(G,1))%op(S(G,2),S(G,3))

	ZGFAT44=ZGFAT44+CG44(S(G,1))%op(S(G,2),S(G,3))/ZDENO
	
400	END DO

        DO 500 G=1,LL
        NBX=IB1CDOWNSISX(SS(G,1))
        NBY=IB1CDOWNSISY(SS(G,1))
       ZDENO=ZWD + WPHOS(NBY)%op3(SS(G,3)) - WPHOS(NBX)%op3(SS(G,2))
!        ZDENO=ZWD - WPHOS(NBY)%op3(SS(G,3)) + WPHOS(NBX)%op3(SS(G,2))

!        PRINT*,SS(G,1),SS(G,2),SS(G,3),CG22(SS(G,1))%op(SS(G,2),SS(G,3))

        ZGFAT22=ZGFAT22+CG22(SS(G,1))%op(SS(G,2),SS(G,3))/ZDENO

500     END DO

        DO 600 G=1,LLL
        NBX=IB1CDOWNSISX(SSS(G,1))
        NBY=IB1CDOWNSISY(SSS(G,1))
       ZDENO=ZWD + WPHOS(NBY)%op3(SSS(G,3)) - WPHOS(NBX)%op3(SSS(G,2))
!        ZDENO=ZWD - WPHOS(NBY)%op3(SSS(G,3)) + WPHOS(NBX)%op3(SSS(G,2))

!        PRINT*,SSS(G,1),SSS(G,2),SSS(G,3),CG24(SSS(G,1))%op(SSS(G,2),SSS(G,3))

        ZGFAT24=ZGFAT24+CG24(SSS(G,1))%op(SSS(G,2),SSS(G,3))/ZDENO

600     END DO		!#BLOQUES-CDOWNSIS
    ENDIF
	
!	ZGAD11=ZGFATD11
!	ZGAD13=ZGFATD13
!	ZGAD31=ZGFATD13
!	ZGAD33=ZGFATD33
!	
!	ZGAD=ZGAD11+ZGAD13+ZGAD31+ZGAD33
!!	Z1AD=1.D0+ACOPt*ZGAD11
!!	Z1BD=1.D0+ACOPt*ZGAD33
!!	Z1CD=ACOPt*ZGAD13
!!	Z1DD=ACOPt*ZGAD31
!	ZDETAD=1.D0-ACOPt*ZGAD
!	ZGDD=ZGAD11*ZGAD33-ZGAD13*ZGAD31
!	ZMD11=(ZGAD11+ACOPt*ZGDD)/ZDETAD
!	ZMD13=(ZGAD13-ACOPt*ZGDD)/ZDETAD
!	ZMD31=(ZGAD31-ACOPt*ZGDD)/ZDETAD
!	ZMD33=(ZGAD33+ACOPt*ZGDD)/ZDETAD
!	ZGAMAD=ZMD11+ZMD13+ZMD31+ZMD33
!         
!!	ZM11TOTAL=ZMD11
!!	ZM13TOTAL=ZMD13
!!	ZM31TOTAL=ZMD31
!!	ZM33TOTAL=ZMD33
!!	ZGAMATOTAL=ZM11TOTAL+ZM13TOTAL+ZM33TOTAL
!         
!	ZMID=1.D0/ZGAMAD
!	ZMID2=ZMID*ZMID
!	ZARGUID=CDSQRT(ZMID2-4.D0*ACOPtR*ACOPtR)
!	IF(DIMAG(ZARGUID).LT.0.D0)ZARGUID=DCONJG(ZARGUID)
!	ZINTD=1.D0/ZARGUID
!         
!	ZGD11=(ZG110+ZMD11*ZINTD/ZGAMAD)
!	ZGD13=(ZG130+ZMD13*ZINTD/ZGAMAD)
!	ZGD31=(ZG310+ZMD31*ZINTD/ZGAMAD)
!	ZGD33=(ZG330+ZMD33*ZINTD/ZGAMAD)
!	ZGFD=ZGD11+ZGD13+ZGD31+ZGD33
         
	ZGA11=ZGFAT11
	ZGA13=ZGFAT13
	ZGA31=ZGFAT13
	ZGA33=ZGFAT33

      ZGA=ZGA11+ZGA13+ZGA31+ZGA33
      ZTETA=1.D0-ACOPt*ZGA
      ZGD=ZGA11*ZGA33-ZGA13*ZGA31
!      IF(IZI.EQ.1)THEN
!      ZM11=ZM110+(ZGA11+ACOPt*ZGD)/ZTETA
!      ZM13=ZM130+(ZGA13-ACOPt*ZGD)/ZTETA
!      ZM31=ZM310+(ZGA31-ACOPt*ZGD)/ZTETA
!      ZM33=ZM330+(ZGA33+ACOPt*ZGD)/ZTETA
!      ELSE       	     

!      ZM11=(ZGA11+ACOPt*ZGD)/ZTETA
!      ZM13=(ZGA13-ACOPt*ZGD)/ZTETA
!      ZM31=(ZGA31-ACOPt*ZGD)/ZTETA
!      ZM33=(ZGA33+ACOPt*ZGD)/ZTETA  
      
      ZM11=ZGA11
      ZM13=ZGA13
      ZM31=ZGA31
      ZM33=ZGA33 
      
         
      ZGAMA=ZM11+ZM13+ZM31+ZM33
 !     ENDIF
      
      
 
         ZMI=1.D0/ZGAMA
         ZMI2=ZMI*ZMI
         ZARGUI=CDSQRT(ZMI2-4.D0*ACOPtR*ACOPtR)
         IF(DIMAG(ZARGUI).LT.0.D0)ZARGUI=DCONJG(ZARGUI)
         ZINT=1.D0/ZARGUI
         ZG11=ZM11*ZINT/ZGAMA
         ZG13=ZM13*ZINT/ZGAMA
         ZG31=ZM31*ZINT/ZGAMA
         ZG33=ZM33*ZINT/ZGAMA
      
!       ZTETA=ZM11*ZM33-ZM13*ZM31
!       ZARGUI=(1.D0+(D-AMU)*ZGAMA)/(1.D0-(D+AMU)*ZGAMA)
!          ZINT=(0.5D0/D)*CDLOG(ZARGUI)
! !         IF(DIMAG(ZINT).GT.0.D0)ZINT=DCONJG(ZINT)
!          ZG11=ZG110+ZTETA/ZGAMA+(ZM11-ZTETA/ZGAMA)*ZINT/ZGAMA
!          ZG13=ZG130-ZTETA/ZGAMA+(ZM13+ZTETA/ZGAMA)*ZINT/ZGAMA
!          ZG31=ZG310-ZTETA/ZGAMA+(ZM31+ZTETA/ZGAMA)*ZINT/ZGAMA
!          ZG33=ZG330+ZTETA/ZGAMA+(ZM33-ZTETA/ZGAMA)*ZINT/ZGAMA
         ZGFP=ZG11+ZG13+ZG31+ZG33
         
    IF(IDOWN.EQ.1)THEN
      ZGA22=ZGFAT22
      ZGA24=ZGFAT24
      ZGA42=ZGFAT24
      ZGA44=ZGFAT44

      ZM22=ZGA22
      ZM24=ZGA24
      ZM42=ZGA42
      ZM44=ZGA44

!      ZGA=ZGA22+ZGA24+ZGA42+ZGA44
!      ZTETA=1.D0-ACOPt*ZGA
!      ZGD=ZGA22*ZGA44-ZGA24*ZGA42
!      ZM22=(ZGA22+ACOPt*ZGD)/ZTETA
!      ZM24=(ZGA24-ACOPt*ZGD)/ZTETA
!      ZM42=(ZGA42-ACOPt*ZGD)/ZTETA
!      ZM44=(ZGA44+ACOPt*ZGD)/ZTETA
!      ZGAMA=ZM22+ZM24+ZM42+ZM44
         
         ZMI=1.D0/ZGAMA
         ZMI2=ZMI*ZMI
         ZARGUI=CDSQRT(ZMI2-4.D0*ACOPtR*ACOPtR)
         IF(DIMAG(ZARGUI).LT.0.D0)ZARGUI=DCONJG(ZARGUI)
         ZINT=1.D0/ZARGUI
         ZG22=ZM22*ZINT/ZGAMA
         ZG24=ZM24*ZINT/ZGAMA
         ZG42=ZM42*ZINT/ZGAMA
         ZG44=ZM44*ZINT/ZGAMA
         
!          ZTETA=ZM22*ZM44-ZM24*ZM42
!          ZARGUI=(1.D0+(D-AMU)*ZGAMA)/(1.D0-(D+AMU)*ZGAMA)
!          ZINT=(0.5D0/D)*CDLOG(ZARGUI)
! !         IF(DIMAG(ZINT).GT.0.D0)ZINT=DCONJG(ZINT)
!          ZG22=ZG220+ZTETA/ZGAMA+(ZM22-ZTETA/ZGAMA)*ZINT/ZGAMA
!          ZG24=ZG240-ZTETA/ZGAMA+(ZM24+ZTETA/ZGAMA)*ZINT/ZGAMA
!          ZG42=ZG420-ZTETA/ZGAMA+(ZM42+ZTETA/ZGAMA)*ZINT/ZGAMA
!          ZG44=ZG440+ZTETA/ZGAMA+(ZM44-ZTETA/ZGAMA)*ZINT/ZGAMA
         ZGFN=ZG22+ZG24+ZG42+ZG44
    END IF
         
	END SUBROUTINE GKONDO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE RESIDUOS
	USE PRINCIPAL
	USE AMODO
	USE AMODHSIS
	USE FGREEN
	IMPLICIT NONE
!	COMPLEX(KIND=KIND(1.0D0)) :: ZWDI,ZGF0,ZG11,ZG33,ZGFD,ZGF

	ACOPt=ACOPtI
	ACOPU=ACOPUI
	ACOPn=ACOPnI
	ACOPV=ACOPVI
	ACPIHM=ACPI
	AMU=AMUI
	BETA=BETAI

	!IONIC HUBBARD MODEL
	ACS(1)=-ACPIHM
	ACS(2)=ACPIHM
	ACS(3)=-ACPIHM
	ACS(4)=ACPIHM
	ACS(5)=-ACPIHM
	ACS(6)=ACPIHM
	ACS(7)=-ACPIHM
	ACS(8)=ACPIHM
	ACS(9)=-ACPIHM
	ACS(10)=ACPIHM
	
	
	CALL VALIN 
	
90	LAT=1
	LBL=0
!	ITES=0
!*************           INFINITO                  ***************
100	LAT=LAT+1
	LBL=LBL+1
	     
	PRINT*,'#############','RED=',LAT,'######'
!	WRITE(44,*)'#############','RED=',LAT,'######'

	IF(LAT>2)CALL LEE
	
	CALL ESTRUCTURA
	CALL CONECTIONS   
	CALL MULTI
	CALL THSIS

         
	IF(LAT.EQ.LFIN)THEN
	
!	 CALL MEOCUP
!	 IF(IDOWN.EQ.1)CALL MEOCDOWN
	 
	 GO TO 900
	END IF

	CALL TCUPES1
	IF(IDOWN.EQ.1)CALL TCDOWNES1
	CALL TNSIS
	CALL TCUPSIS
	CALL TCUPSISQ
	CALL TCDOWNSIS
	CALL TCDOWNSISQ
	CALL DEALLOC_IGSIS
	CALL DEALLOC_ESP1
	CALL TRANSFNEW
  
	GOTO 100	
	
900	PRINT*,'OK-RESIDUOS'
	
	END SUBROUTINE RESIDUOS
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 	SUBROUTINE OCUP(AMU,SDS,SVAC,SF,SD2)
! 	IMPLICIT NONE
! 	EXTERNAL FERM
! 	EXTERNAL ZFER
! 	COMPLEX (kind=8):: Z,ZGA,ZG11,ZG33,ZGF,ZGF0,ZGFD,ZFER,ZF
! 	COMPLEX (kind=8):: ZG13,ZG31,ZG33,ZGFP,ZG22,ZG24,ZG42,ZG44,ZGFN
! 	real(kind=8):: SD,SDS,A,B,R,H1,EPSIL,NPTS,ICHECK,RELERR
! 	real(kind=8):: ETTA,CC,FF1,FF2,FF3,FF4,FERM,SF0,SF,SVAC,SD2,AMU
! 	real(kind=8):: RES1,RES2,RES3,RES4,PI,VA(5),VB(5)
! 	INTEGER:: I,J,N1
! 	DIMENSION FF1(2001),FF2(2001),FF3(2001),FF4(2001)
! 	DATA PI /3.14159265358979323D0/
! 	CC=-1.D0/PI
! 	EPSIL=1.D-4
! 	ETTA=1.D-4
! 	SD=0.D0
! 	SF0=0.D0
! 	SVAC=0.D0
! 	SD2=0.D0
! 	N1=2001
! 	VA(1)=-5.0D0
! 	VB(1)=-1.0D0
! 	VA(2)=-1.0D0
! 	VB(2)=-0.20D0
! 	VA(3)=-0.20D0
! 	VB(3)=0.20D0
! 	VA(4)=0.20D0
! 	VB(4)=1.0D0
! 	VA(5)=1.0D0
! 	VB(5)=5.0D0
! 	DO 10 I=1,5
! 	 A=VA(I)
! 	 B=VB(I)
! 	 R=B-A
! 	 H1=R/(N1-1)
! 	 DO 1 J=2,N1
! 	 Z=DCMPLX(A+(J-1)*H1,ETTA)
! 	 CALL GKONDO(Z,ZGF0,ZG11,ZG13,ZG31,ZG33,ZGFP,ZG22,ZG24,ZG42,ZG44,ZGFN)
! !*******************************************
! 	 ZF=ZFER(Z)
! !	 ZGF=ZGF0
! 	 FF1(J)=DIMAG(ZGF)
! 	 FF2(J)=DIMAG(ZGF*ZF)
! 	 FF3(J)=DIMAG(ZGF*(1.D0-ZF))
! 1	 FF4(J)=DIMAG(ZG33*ZF)
! 	 CALL SIMP(N1,H1,FF1,RES1)
! 	 CALL SIMP(N1,H1,FF2,RES2)
! 	 CALL SIMP(N1,H1,FF3,RES3)
! 	 CALL SIMP(N1,H1,FF4,RES4)
! !***********************************
! 	 SD=SD+RES1
! 	 SF0=SF0+RES2
! 	 SVAC=SVAC+RES3
! 	 SD2=SD2+RES4
! 10	 CONTINUE
! 	 SDS=CC*SD
! 	 SF=CC*SF0
! 	 SVAC=CC*SVAC
! 	 SD2=CC*SD2
! 	END SUBROUTINE OCUP
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        SUBROUTINE EXATA

!       ESTA ROTINA CALCULA OS RESÍDUOS DO DÍMERO DE FORMA ANALÍTICA

	USE PRINCIPAL
	USE AMODO
	USE FGREEN
	IMPLICIT NONE
	EXTERNAL FERM
	REAL(KIND=8):: FERM
	REAL(KIND=8):: EMINDIM,AZ,E(16),DXP(16),TT
!        REAL(KIND=8):: DM11(12),DM13(12),DM33(12),DUT(12)
!        REAL(KIND=8):: AMU,TH0,TH,E0,U,BETA
	REAL(KIND=8):: A,B,C,R2
	REAL(KIND=8):: C01,C02,C03,C04
	INTEGER:: I,KK

	TT=1.D0/BETA
	R2=DSQRT(2.D0)
	C=DSQRT((16.D0*ACOPt*ACOPt)+(ACOPU*ACOPU))
	A=DSQRT(2.D0*(((16.D0*ACOPt*ACOPt)/((C-ACOPU)*(C-ACOPU)))+1.D0))
	B=DSQRT(2.D0*(((16.D0*ACOPt*ACOPt)/((C+ACOPU)*(C+ACOPU)))+1.D0))
	C01=(4.0D0*ACOPt)/(B*(C+ACOPU)*R2)
	C02=(1.0D0)/(B*R2)
	C03=(4.0D0*ACOPt)/(A*(C-ACOPU)*R2)
	C04=(1.0D0)/(A*R2)

	E(1)=0.D0
	E(2)=ACOPn+ACOPt
	E(3)=ACOPn+ACOPt
	E(4)=ACOPn-ACOPt
	E(5)=ACOPn-ACOPt
	E(6)=2.D0*ACOPn
	E(7)=2.D0*ACOPn
	E(8)=2.D0*ACOPn
	E(9)=2.D0*ACOPn+ACOPU
	E(10)=2.D0*ACOPn+((ACOPU+C)/2.D0)
	E(11)=2.D0*ACOPn+((ACOPU-C)/2.D0)
	E(12)=3.D0*ACOPn+ACOPU+ACOPt
	E(13)=3.D0*ACOPn+ACOPU+ACOPt
	E(14)=3.D0*ACOPn+ACOPU-ACOPt
	E(15)=3.D0*ACOPn+ACOPU-ACOPt
	E(16)=4.D0*ACOPn+2.D0*ACOPU
	EMINDIM=E(1)
	DO 50 I=2,16
	IF(EMINDIM.GT.E(I)) EMINDIM=E(I)
50	CONTINUE
	KK=0
	DO 51 I=1,16
	KK=KK+1
	IF(EMINDIM.EQ.E(I))GO TO 52
51	CONTINUE
52	CONTINUE

!---------------------------PARTITION FUNCTION---------------------
	AZ=0.0D0
	DO I=1,16
	E(I)=E(I)-EMINDIM
	AZ=AZ+DEXP(-BETA*E(I))
	ENDDO
	DUT(1)=ACOPn+ACOPt
	DUT(2)=ACOPn-ACOPt
	DUT(3)=ACOPn+ACOPt+ACOPU
	DUT(4)=ACOPn-ACOPt+ACOPU
	DUT(5)=ACOPn+ACOPt+((ACOPU+C)/2.D0)
	DUT(6)=ACOPn+ACOPt+((ACOPU-C)/2.D0)
	DUT(7)=ACOPn-ACOPt+((ACOPU+C)/2.D0)
	DUT(8)=ACOPn-ACOPt+((ACOPU-C)/2.D0)
	DUT(9)=ACOPn+ACOPt+ACOPU-((ACOPU+C)/2.D0)
	DUT(10)=ACOPn+ACOPt+ACOPU-((ACOPU-C)/2.D0)
	DUT(11)=ACOPn-ACOPt+ACOPU-((ACOPU+C)/2.D0)
	DUT(12)=ACOPn-ACOPt+ACOPU-((ACOPU-C)/2.D0)
	DO 61 I=1,16
	DXP(I) = DEXP(-BETA*E(I))
61	CONTINUE
	DM11(1)=(0.5D0*(DXP(1)+DXP(2))+0.5D0*(DXP(4)+DXP(6))+ &
	0.25D0*(DXP(5)+DXP(8))+0.25D0*(DXP(9)+DXP(12)))/AZ
	DM11(2)=(0.5D0*(DXP(1)+DXP(4))+0.5D0*(DXP(2)+DXP(6))+ &
	0.25D0*(DXP(3)+DXP(8))+0.25D0*(DXP(9)+DXP(14)))/AZ
	DM11(3)=0.0D0
	DM11(4)=0.0D0
	DM11(5)=(C01*C01*(DXP(5)+DXP(10)))/AZ
	DM11(6)=(C03*C03*(DXP(5)+DXP(11)))/AZ
	DM11(7)=(C01*C01*(DXP(3)+DXP(10)))/AZ
	DM11(8)=(C03*C03*(DXP(3)+DXP(11)))/AZ
	DM11(9)=(C02*C02*(DXP(10)+DXP(12)))/AZ
	DM11(10)=(C04*C04*(DXP(11)+DXP(12)))/AZ
	DM11(11)=(C02*C02*(DXP(10)+DXP(14)))/AZ
	DM11(12)=(C04*C04*(DXP(11)+DXP(14)))/AZ
        
	DM13(1)=0.0D0
	DM13(2)=0.0D0
	DM13(3)=0.0D0
	DM13(4)=0.0D0
	DM13(5)=((C02*(-C01))*(DXP(5)+DXP(10)))/AZ
	DM13(6)=((C04*C03)*(DXP(5)+DXP(11)))/AZ
	DM13(7)=((C02*C01)*(DXP(3)+DXP(10)))/AZ
	DM13(8)=((C04*(-C03))*(DXP(3)+DXP(11)))/AZ
	DM13(9)=((C02*(-C01))*(DXP(10)+DXP(12)))/AZ
	DM13(10)=((C04*C03)*(DXP(11)+DXP(12)))/AZ
	DM13(11)=((C02*C01)*(DXP(10)+DXP(14)))/AZ
	DM13(12)=((C04*(-C03))*(DXP(11)+DXP(14)))/AZ

	DM33(1)=0.0D0
	DM33(2)=0.0D0
	DM33(3)=(0.25D0*(DXP(5)+DXP(9))+0.5D0*(DXP(7)+DXP(13))+&
	0.25D0*(DXP(8)+DXP(12))+0.5D0*(DXP(15)+DXP(16)))/AZ
	DM33(4)=(0.25D0*(DXP(3)+DXP(9))+0.5D0*(DXP(7)+DXP(15))+&
	0.25D0*(DXP(8)+DXP(14))+0.5D0*(DXP(13)+DXP(16)))/AZ
	DM33(5)=(C02*C02*(DXP(5)+DXP(10)))/AZ
	DM33(6)=(C04*C04*(DXP(5)+DXP(11)))/AZ
	DM33(7)=(C02*C02*(DXP(3)+DXP(10)))/AZ
	DM33(8)=(C04*C04*(DXP(3)+DXP(11)))/AZ
	DM33(9)=(C01*C01*(DXP(10)+DXP(12)))/AZ
	DM33(10)=(C03*C03*(DXP(11)+DXP(12)))/AZ
	DM33(11)=(C01*C01*(DXP(10)+DXP(14)))/AZ
	DM33(12)=(C03*C03*(DXP(11)+DXP(14)))/AZ
        
	END SUBROUTINE EXATA
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE SIMP(N,H,F,RES)
	USE PRINCIPAL
	IMPLICIT NONE
	REAL(KIND=8):: F
	REAL(KIND=8):: H,RES,SP,SI
	INTEGER:: N,I,M
	DIMENSION F(N)
	M=N-3
	SP=F(2)
	SI=F(3)
	DO 1 I=4,M,2
	SP=SP+F(I)
1	SI=SI+F(I+1)
	RES=H*(F(1)+2.D0*(SI+2.D0*(SP+F(N-1)))+F(N))/3.D0
	END SUBROUTINE SIMP
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	REAL(KIND=8) function  FERM(AX)
	USE PRINCIPAL
	IMPLICIT NONE
	REAL(KIND=8):: AX,BEXP,ARRG,T,ABSX
!	REAL(KIND=8):: AMU,TH,E0,U,BETA
!	COMMON/DADOS/AMU,TH,E0,U,BETA
	BEXP=35.0D0
	T=1.D0/BETA
	ABSX=DABS(AX)
	IF(ABSX - T*BEXP)250,200,200
200	IF(AX)210,220,230
210	FERM=1.0D0
	RETURN
220	FERM=0.5D0
	RETURN
230	FERM=0.0D0
	RETURN
250	ARRG=AX*BETA
	FERM=1.0D0/(1.0D0+DEXP(ARRG))
300	RETURN
	END FUNCTION FERM
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	COMPLEX (KIND=8) function  ZFER(ZW)
	USE PRINCIPAL
	IMPLICIT NONE
	EXTERNAL FERM
	COMPLEX (kind=8):: ZW,Z,ZEX,ZUM
	REAL(KIND=8):: FERM
	REAL(KIND=8):: BEXP,T,RX,RY,ABSX
!	REAL(KIND=8):: AMU,TH,E0,U,BETA
!	COMMON/DADOS/AMU,TH,E0,U,BETA
	DATA BEXP/35.D0/
	DATA ZUM/(1.D0,0.D0)/
	Z=ZW
	ZEX=Z*BETA
	T=1.D0/BETA
	RX=DREAL(Z)
	RY=DIMAG(Z)
	IF(RY.NE.0.D0)GO TO 100
	ZFER=DCMPLX(FERM(RX),0.D0)
	RETURN
100	ABSX=DABS(RX)
	IF(ABSX - T*BEXP) 250,200,200
200	IF(RX)210,220,230
210	ZFER= (1.0D0,0.D0)
	RETURN
220	ZFER=(0.5D0,0.D0)
	RETURN
230	ZFER=(0.D0,0.D0)
	RETURN
250	CONTINUE
	ZFER=ZUM/(ZUM+CDEXP(ZEX))
	END FUNCTION ZFER
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC6CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE MEOCUP
	USE PRINCIPAL
	USE AMODESP1
	USE AMODESP2
	USE AMODIGSIS
	USE AMODO
	USE AMODHSIS
	USE AMODOLO
	
	IMPLICIT NONE
	INTEGER :: J,JJX,JJY,IU1XX,IU1YY,ICHARXX,ICHARYY,IFIND,IU1T,IEA,IEB
	INTEGER :: ICHART,JRIGHT,JLEFT,IB1,ISIG,IB2,M,NBX,NBY
	INTEGER :: IIB2,JIB2,I1,I2,IX11,J1,J2,IY11,JJ,JJYNEW,J11
	INTEGER :: K,JJXNEW,I11,IIB1,JIB1,NS1,JBTX,JBTY
	REAL(KIND=KIND(1.0D0)) :: AAE,AAE1,AAE2,AAC,AAD,AAF,AZN,ADEXP,AREX33,AREX11,AREX13
	REAL(KIND=KIND(1.0D0)), DIMENSION(:,:), ALLOCATABLE :: H33,H11,H13,QC,QD,QF,HBS,HBD


	
	PRINT*,'ESTOY EN MEOCUP'
	
	PRINT*,'=======>>  SITIO=',IZI !SITIO PARA EL CUAL SE CALCULA LA FUNCION DE GREEN
	
	DO 1000 J=1,NBL1CUPSIS
!	PRINT*,'J=',J
	JBTX=IB1CUPSISX(J)
	JBTY=IB1CUPSISY(J)
	IF(((IBHUSED(JBTX)).EQ.1).AND.((IBHUSED(JBTY)).EQ.1))THEN
!	PRINT*,'PASE--> J=',J,JBTX,JBTY
	!*-> encontrar setores de O^{dagga} e O que conectam com ACUP
	IU1XX=IU1B1XCUPSIS(J)
	IU1YY=IU1B1YCUPSIS(J)

	ICHARXX=ICHARB1XCUPSIS(J)
	ICHARYY=ICHARB1YCUPSIS(J)

	DO IFIND=1,NBL1HSIS
	  IU1T=IU1B1SIS(IFIND)
	  ICHART=ICHARB1SIS(IFIND)

	  IF((IU1T==IU1YY).AND.(ICHART==ICHARYY)) JRIGHT=IFIND

	  IF((IU1T==IU1XX).AND.(ICHART==ICHARXX)) JLEFT=IFIND
	END DO
	
!*->	
	ALLOCATE(CG33(J)%op(IDIMB1CUPSIS(1,J),IDIMB1CUPSIS(2,J)) )
	CG33(J)%op=0.D0
	
	ALLOCATE(CG11(J)%op(IDIMB1CUPSIS(1,J),IDIMB1CUPSIS(2,J)) )
	CG11(J)%op=0.D0
	
	ALLOCATE(CG13(J)%op(IDIMB1CUPSIS(1,J),IDIMB1CUPSIS(2,J)) )
	CG13(J)%op=0.D0
	
	IF( ALLOCATED(H33) )  DEALLOCATE( H33 ) 
	ALLOCATE(H33(IDIMB1CUPSIS(1,J),IDIMB1CUPSIS(2,J)))
	H33=0.D0
	IF( ALLOCATED(H11) )  DEALLOCATE( H11 )
	ALLOCATE(H11(IDIMB1CUPSIS(1,J),IDIMB1CUPSIS(2,J)))
	H11=0.D0
! 	IF( ALLOCATED(H13) )  DEALLOCATE( H13 )
! 	ALLOCATE(H13(IDIMB1CUPSIS(1,J),IDIMB1CUPSIS(2,J)))
! 	H13=0.D0
	
	IF(IZI.NE.LAT)THEN
!	PRINT*,'LEFT-SIDE CHAIN',IZI
!CCCCC ---->     Cup x 1
	DO 500 IB1=1,NBL1CUP
!	ISITEM
	NS1=10000+1000*IZI+IB1
	READ(NS1)IEA,IEB
	ALLOCATE(HBS(IEA,IEB))
	READ(NS1)HBS
	REWIND(NS1)
	CLOSE(NS1)
	
!	IDUTEM
	NS1=20000+1000*IZI+IB1
	READ(NS1)IEA,IEB
	ALLOCATE(HBD(IEA,IEB))
	READ(NS1)HBD
	REWIND(NS1)
	CLOSE(NS1)
	
	DO 400 IB2=1,NBL2H

	IU1XX=IU1B1XCUP(IB1)+ IU1B2(IB2)
	ICHARXX=ICHARB1XCUP(IB1)+ ICHARB2(IB2)

	IF((IU1XX==IU1B1XCUPSIS(J)).AND.(ICHARXX==ICHARB1XCUPSIS(J))) THEN !250

	    IIB1=IB1CUPX(IB1)
	    JIB1=IB1CUPY(IB1)

	  DO 200 I1=1,IDIMB1CUP(1,IB1)
	  DO 200 I2=1,IDIMB2H(IB2)

	    IX11=IGSIS(IIB1,IB2)%op4(I1,I2)

	    DO 100 J1=1,IDIMB1CUP(2,IB1)
	    DO 100 J2=1,IDIMB2H(IB2)

	    IY11=IGSIS(JIB1,IB2)%op4(J1,J2)
!* CUPx1  
	AAE1=0.D0
	AAE2=0.D0
	IF(J2==I2)AAE1=HBS(I1,J1)
	IF(J2==I2)AAE2=HBD(I1,J1)

	H33(IX11,IY11)=H33(IX11,IY11)+AAE2

	H11(IX11,IY11)=H11(IX11,IY11)+AAE1
	
!	H13(IX11,IY11)=H13(IX11,IY11)+AAE
	    
100	    CONTINUE
200	  CONTINUE
250	END IF
400	CONTINUE
	DEALLOCATE(HBS)
	DEALLOCATE(HBD)
!	DEALLOCATE(CUPK(IB1)%op)
500	CONTINUE

	ELSE	!IZI
	
!	PRINT*,'RIGHT-SIDE CHAIN',IZI
	DO 900 IB1=1,NBL1H
	ISIG=(-1)**(ICHARB1(IB1))
	DO 900 IB2=1,NBL2CUP

	IU1XX=IU1B1(IB1)+IU1B2XCUP(IB2)
	ICHARXX=ICHARB1(IB1)+ICHARB2XCUP(IB2)

	IF((IU1XX==IU1B1XCUPSIS(J)).AND.(ICHARXX==ICHARB1XCUPSIS(J))) THEN !850

	    IIB2=IB2CUPX(IB2)
	    JIB2=IB2CUPY(IB2)

	  DO 800 I1=1,IDIMB1H(IB1)
	  DO 800 I2=1,IDIMB2CUP(1,IB2)

	    IX11=IGSIS(IB1,IIB2)%op4(I1,I2)

	    DO 700 J1=1,IDIMB1H(IB1)
	    DO 700 J2=1,IDIMB2CUP(2,IB2)

	    IY11=IGSIS(IB1,JIB2)%op4(J1,J2)
!* 1xCUP  
	AAE=0.D0
	IF(J1==I1)AAE=ISIG*CUPL(IB2)%op(I2,J2)
	
	IF(IIDUPLA2(IB2,I2,J2).EQ.1)THEN
	H33(IX11,IY11)=H33(IX11,IY11)+AAE
	END IF
	
	IF(IISIMPLE2(IB2,I2,J2).EQ.1)THEN
	H11(IX11,IY11)=H11(IX11,IY11)+AAE
	END IF
	
!	H13(IX11,IY11)=H13(IX11,IY11)+AAE
	
700	    CONTINUE
800	  CONTINUE
850	END IF
900	CONTINUE
	END IF		!IZI
	
! 	DO I1=1,IDIMB1CUPSIS(1,J)
! 	DO J1=1,IDIMB1CUPSIS(2,J)
! 	IF(H11(I1,J1)/=0.D0)PRINT*,J,I1,J1,H11(I1,J1)
! 	END DO
! 	END DO	
!CCCCCCCCCCCC   --------> ACUP_{B}=O^{+}ACUP_{Ba}O
	JJX=IDIMB1CUPSIS(1,J)
	JJ=IDIMB1HSIS(JRIGHT)
	JJYNEW=IDIMB1HSIS(JRIGHT)

	ALLOCATE(QC(JJX,JJYNEW) )
	QC=0.D0
	ALLOCATE(QD(JJX,JJYNEW) )
	QD=0.D0
!	ALLOCATE(QF(JJX,JJYNEW) )
!	QF=0.D0
	

	DO I1=1,JJX
	DO J1=1,JJYNEW
	  AAC=0.D0
	  AAD=0.D0
	  AAF=0.D0
	  DO K=1,JJ
	    AAC=AAC + H33(I1,K)*OS(JRIGHT)%op(K,J1)
	    AAD=AAD + H11(I1,K)*OS(JRIGHT)%op(K,J1)
!	    AAF=AAF + H13(I1,K)*OS(JRIGHT)%op(K,J1)
	  END DO
	  QC(I1,J1)=AAC
	  QD(I1,J1)=AAD
!	  QF(I1,J1)=AAF
	END DO
	END DO

	DEALLOCATE(H33)
	DEALLOCATE(H11)
!	DEALLOCATE(H13)

!* O^{+}([ACUP]O)
	JJ=IDIMB1HSIS(JLEFT)
	JJXNEW=IDIMB1HSIS(JLEFT)
	DO I1=1,JJXNEW
	DO J1=1,JJYNEW
	  AAC=0.D0
	  AAD=0.D0
	  AAF=0.D0
	  DO K=1,JJ
	    AAC=AAC + (OS(JLEFT)%op(K,I1))*QC(K,J1)
	    AAD=AAD + (OS(JLEFT)%op(K,I1))*QD(K,J1)
!	    AAF=AAF + (OS(JLEFT)%op(K,I1))*QF(K,J1)
	  END DO
	  CG33(J)%op(I1,J1)=AAC
	  CG11(J)%op(I1,J1)=AAD
!	  CG13(J)%op(I1,J1)=AAF
	END DO
	END DO
	DEALLOCATE(QC)
	DEALLOCATE(QD)
!	DEALLOCATE(QF)
	
! 	PRINT*,'%%%%%%%%%%%%     RESULT ROTATION   %%%%%%%%%'
! 	DO I1=1,IDIMB1CUPSIS(1,J)
! 	DO J1=1,IDIMB1CUPSIS(2,J)
! 	IF(CG33(J)%op(I1,J1)/=0.D0)PRINT*,J,I1,J1,CG33(J)%op(I1,J1)
! 	END DO
! 	END DO

    END IF ! Solo los bloques usados
1000	END DO		!#BLOQUES-CUPSIS
    

!----> Calculating the partition function of the system
	PRINT*,'BETA=',BETA,'EMIN=',EMIN
	AZN=0.D0
	DO J=1,NBL1HSIS
	IF((IBHUSED(J)).EQ.1)THEN
!	PRINT*,'Partition function',J
!    IF (ICHARB1SIS(J) .NE. LFIN)THEN
	DO M=1,IDIMB1HSIS(J)
	AZN=AZN+DEXP(-BETA*(WPHOS(J)%op3(M)-EMIN))
!	PRINT*,'AZ_LIMIT',J,M,WPHOS(J)%op3(M),EMIN,DEXP(-BETA*(WPHOS(J)%op3(M)-EMIN)),AZN
	END DO
!	ELSE
!	AZN=AZN+DEXP(-BETA*(WPHOS(J)%op3(1)-EMIN))
!	PRINT*,'AZ_LIMIT',J,1,WPHOS(J)%op3(1),EMIN,DEXP(-BETA*(WPHOS(J)%op3(1)-EMIN)),AZN	
!	END IF
	END IF
	END DO
	PRINT*,'partition function AZ=',AZN
	
        P=0
        PP=0
        PPP=0

	DO 2000 J=1,NBL1CUPSIS
!	J=6
!	PRINT*,'J=',J
	JBTX=IB1CUPSISX(J)
	JBTY=IB1CUPSISY(J)
	IF(((IBHUSED(JBTX)).EQ.1).AND.((IBHUSED(JBTY)).EQ.1))THEN
!	PRINT*,'PASE--> J=',J,JBTX,JBTY
	NBX=IB1CUPSISX(J)
	NBY=IB1CUPSISY(J)
!	PRINT*,J,IB1CUPSISX(J),IB1CUPSISY(J)
!	PRINT*,IDIMB1CUPSIS(1,J),IDIMB1CUPSIS(2,J)
	DO I1=1,IDIMB1CUPSIS(1,J)
	DO J1=1,IDIMB1CUPSIS(2,J)
	ADEXP=DEXP(-BETA*(WPHOS(NBX)%op3(I1)-EMIN))+DEXP(-BETA*(WPHOS(NBY)%op3(J1)-EMIN))
	ADEXP=ADEXP/AZN
	CG13(J)%op(I1,J1)=(CG33(J)%op(I1,J1)*CG11(J)%op(I1,J1))
	AREX33=ADEXP*((CG33(J)%op(I1,J1))**2)
	CG33(J)%op(I1,J1)=AREX33
        IF(ABS(CG33(J)%op(I1,J1)) .GT. 1.D-12) THEN
        P=P+1
        END IF
	AREX11=ADEXP*((CG11(J)%op(I1,J1))**2)
	CG11(J)%op(I1,J1)=AREX11
        IF(ABS(CG11(J)%op(I1,J1)) .GT. 1.D-12) THEN
        PP=PP+1
        END IF
!	AREX13=ADEXP*((CG13(J)%op(I1,J1))**2)
	AREX13=ADEXP*CG13(J)%op(I1,J1)
	CG13(J)%op(I1,J1)=AREX13
        IF(ABS(CG13(J)%op(I1,J1)) .GT. 1.D-12) THEN
        PPP=PPP+1
        END IF
!	PRINT*,I1,J1,CG33(J)%op(I1,J1),CG11(J)%op(I1,J1),CG13(J)%op(I1,J1)
!	PRINT*,I1,J1,WPHOS(NBX)%op3(I1),WPHOS(NBY)%op3(J1),ADEXP
	END DO
	END DO
	
	END IF
2000	END DO		!#BLOQUES-CUPSIS
	
        IF(ALLOCATED(V)) DEALLOCATE(V)
        ALLOCATE(V(P,3))
        V=0.D0
        IF(ALLOCATED(VV)) DEALLOCATE(VV)
        ALLOCATE(VV(PP,3))
        VV=0.D0
        IF(ALLOCATED(VVV)) DEALLOCATE(VVV)
        ALLOCATE(VVV(PPP,3))
        VVV=0.D0

        P=0
        PP=0
        PPP=0

        DO 6000 J=1,NBL1CUPSIS
        JBTX=IB1CUPSISX(J)
        JBTY=IB1CUPSISY(J)
        IF(((IBHUSED(JBTX)).EQ.1).AND.((IBHUSED(JBTY)).EQ.1))THEN
        NBX=IB1CUPSISX(J)
        NBY=IB1CUPSISY(J)
        DO I1=1,IDIMB1CUPSIS(1,J)
        DO J1=1,IDIMB1CUPSIS(2,J)
        IF(ABS(CG33(J)%op(I1,J1)) .GT. 1.D-12) THEN
        V(P+1,1)=J
        V(P+1,2)=I1
        V(P+1,3)=J1
        P=P+1
        END IF
        IF(ABS(CG11(J)%op(I1,J1)) .GT. 1.D-12) THEN
        VV(PP+1,1)=J
        VV(PP+1,2)=I1
        VV(PP+1,3)=J1
        PP=PP+1
        END IF
        IF(ABS(CG13(J)%op(I1,J1)) .GT. 1.D-12) THEN
        VVV(PPP+1,1)=J
        VVV(PPP+1,2)=I1
        VVV(PPP+1,3)=J1
        PPP=PPP+1
        END IF
        END DO
        END DO
        END IF
6000    END DO

	PRINT*,'O.K.'
	
	
	END SUBROUTINE MEOCUP
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC6CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE MEOCDOWN
	USE PRINCIPAL
	USE AMODESP1
	USE AMODESP2
	USE AMODIGSIS
	USE AMODO
	USE AMODHSIS
	USE AMODOLO
	
	IMPLICIT NONE
	INTEGER :: J,JJX,JJY,IU1XX,IU1YY,ICHARXX,ICHARYY,IFIND,IU1T,IEA,IEB
	INTEGER :: ICHART,JRIGHT,JLEFT,IB1,ISIG,IB2,M,NBX,NBY
	INTEGER :: IIB2,JIB2,I1,I2,IX11,J1,J2,IY11,JJ,JJYNEW,J11
	INTEGER :: K,JJXNEW,I11,IIB1,JIB1,NS1,JBTX,JBTY
	REAL(KIND=KIND(1.0D0)) :: AAE,AAE1,AAE2,AAC,AAD,AAF,AZN,ADEXP,AREX44,AREX22,AREX24
	REAL(KIND=KIND(1.0D0)), DIMENSION(:,:), ALLOCATABLE :: H44,H22,H24,QC,QD,QF,HBS,HBD


	
	PRINT*,'ESTOY EN MEOCDOWN'
	
	PRINT*,'=======>>  SITIO=',IZI !SITIO PARA EL CUAL SE CALCULA LA FUNCION DE GREEN
	
	DO 1000 J=1,NBL1CDOWNSIS
!	PRINT*,'J=',J
	JBTX=IB1CDOWNSISX(J)
	JBTY=IB1CDOWNSISY(J)
	IF(((IBHUSED(JBTX)).EQ.1).AND.((IBHUSED(JBTY)).EQ.1))THEN
!	PRINT*,'PASE--> J=',J,JBTX,JBTY
!*-> encontrar setores de O^{dagga} e O que conectam com ACDOWN
	IU1XX=IU1B1XCDOWNSIS(J)
	IU1YY=IU1B1YCDOWNSIS(J)

	ICHARXX=ICHARB1XCDOWNSIS(J)
	ICHARYY=ICHARB1YCDOWNSIS(J)

	DO IFIND=1,NBL1HSIS
	  IU1T=IU1B1SIS(IFIND)
	  ICHART=ICHARB1SIS(IFIND)

	  IF((IU1T==IU1YY).AND.(ICHART==ICHARYY)) JRIGHT=IFIND

	  IF((IU1T==IU1XX).AND.(ICHART==ICHARXX)) JLEFT=IFIND
	END DO
	
!*->	
	ALLOCATE(CG44(J)%op(IDIMB1CDOWNSIS(1,J),IDIMB1CDOWNSIS(2,J)) )
	CG44(J)%op=0.D0
	
	ALLOCATE(CG22(J)%op(IDIMB1CDOWNSIS(1,J),IDIMB1CDOWNSIS(2,J)) )
	CG22(J)%op=0.D0
	
	ALLOCATE(CG24(J)%op(IDIMB1CDOWNSIS(1,J),IDIMB1CDOWNSIS(2,J)) )
	CG24(J)%op=0.D0
	
	IF( ALLOCATED(H44) )  DEALLOCATE( H44 ) 
	ALLOCATE(H44(IDIMB1CDOWNSIS(1,J),IDIMB1CDOWNSIS(2,J)))
	H44=0.D0
	IF( ALLOCATED(H22) )  DEALLOCATE( H22 )
	ALLOCATE(H22(IDIMB1CDOWNSIS(1,J),IDIMB1CDOWNSIS(2,J)))
	H22=0.D0
! 	IF( ALLOCATED(H24) )  DEALLOCATE( H24 )
! 	ALLOCATE(H24(IDIMB1CDOWNSIS(1,J),IDIMB1CDOWNSIS(2,J)))
! 	H24=0.D0
	
	IF(IZI.NE.LAT)THEN
!	PRINT*,'LEFT-SIDE CHAIN',IZI
!CCCCC ---->     CDOWN x 1
	DO 500 IB1=1,NBL1CDOWN
!	ISITEM
	NS1=30000+1000*IZI+IB1
	READ(NS1)IEA,IEB
	ALLOCATE(HBS(IEA,IEB))
	READ(NS1)HBS
	REWIND(NS1)
	CLOSE(NS1)
	
!	IDUTEM
	NS1=40000+1000*IZI+IB1
	READ(NS1)IEA,IEB
	ALLOCATE(HBD(IEA,IEB))
	READ(NS1)HBD
	REWIND(NS1)
	CLOSE(NS1)
	
	DO 400 IB2=1,NBL2H

	IU1XX=IU1B1XCDOWN(IB1)+ IU1B2(IB2)
	ICHARXX=ICHARB1XCDOWN(IB1)+ ICHARB2(IB2)

	IF((IU1XX==IU1B1XCDOWNSIS(J)).AND.(ICHARXX==ICHARB1XCDOWNSIS(J))) THEN !250

	    IIB1=IB1CDOWNX(IB1)
	    JIB1=IB1CDOWNY(IB1)

	  DO 200 I1=1,IDIMB1CDOWN(1,IB1)
	  DO 200 I2=1,IDIMB2H(IB2)

	    IX11=IGSIS(IIB1,IB2)%op4(I1,I2)

	    DO 100 J1=1,IDIMB1CDOWN(2,IB1)
	    DO 100 J2=1,IDIMB2H(IB2)

	    IY11=IGSIS(JIB1,IB2)%op4(J1,J2)
!* CDOWNx1  
	AAE1=0.D0
	AAE2=0.D0
	IF(J2==I2)AAE1=HBS(I1,J1)
	IF(J2==I2)AAE2=HBD(I1,J1)

	H44(IX11,IY11)=H44(IX11,IY11)+AAE2

	H22(IX11,IY11)=H22(IX11,IY11)+AAE1
	
!	H24(IX11,IY11)=H24(IX11,IY11)+AAE
	    
100	    CONTINUE
200	  CONTINUE
250	END IF
400	CONTINUE
	DEALLOCATE(HBS)
	DEALLOCATE(HBD)
!	DEALLOCATE(CUPK(IB1)%op)
500	CONTINUE

	ELSE	!IZI
!!      CERAR H44??	
!	PRINT*,'RIGHT-SIDE CHAIN',IZI
	DO 900 IB1=1,NBL1H
	ISIG=(-1)**(ICHARB1(IB1))
	DO 900 IB2=1,NBL2CDOWN

	IU1XX=IU1B1(IB1)+IU1B2XCDOWN(IB2)
	ICHARXX=ICHARB1(IB1)+ICHARB2XCDOWN(IB2)

	IF((IU1XX==IU1B1XCDOWNSIS(J)).AND.(ICHARXX==ICHARB1XCDOWNSIS(J))) THEN !850

	    IIB2=IB2CDOWNX(IB2)
	    JIB2=IB2CDOWNY(IB2)

	  DO 800 I1=1,IDIMB1H(IB1)
	  DO 800 I2=1,IDIMB2CDOWN(1,IB2)

	    IX11=IGSIS(IB1,IIB2)%op4(I1,I2)

	    DO 700 J1=1,IDIMB1H(IB1)
	    DO 700 J2=1,IDIMB2CDOWN(2,IB2)

	    IY11=IGSIS(IB1,JIB2)%op4(J1,J2)
!* 1xCDOWN  
	AAE=0.D0
	IF(J1==I1)AAE=ISIG*CDOWNL(IB2)%op(I2,J2)
	
	IF(FFDUPLA2(IB2,I2,J2).EQ.1)THEN
	H44(IX11,IY11)=H44(IX11,IY11)+AAE
	END IF
	
	IF(FFSIMPLE2(IB2,I2,J2).EQ.1)THEN
	H22(IX11,IY11)=H22(IX11,IY11)+AAE
	END IF
	
!	H24(IX11,IY11)=H24(IX11,IY11)+AAE
	
700	    CONTINUE
800	  CONTINUE
850	END IF
900	CONTINUE
	END IF		!IZI
	
! 	DO I1=1,IDIMB1CDOWNSIS(1,J)
! 	DO J1=1,IDIMB1CDOWNSIS(2,J)
! 	IF(H22(I1,J1)/=0.D0)PRINT*,J,I1,J1,H22(I1,J1)
! 	END DO
! 	END DO	
!CCCCCCCCCCCC   --------> ACDOWN_{B}=O^{+}ACDOWN_{Ba}O
	JJX=IDIMB1CDOWNSIS(1,J)
	JJ=IDIMB1HSIS(JRIGHT)
	JJYNEW=IDIMB1HSIS(JRIGHT)

	ALLOCATE(QC(JJX,JJYNEW) )
	QC=0.D0
	ALLOCATE(QD(JJX,JJYNEW) )
	QD=0.D0
!	ALLOCATE(QF(JJX,JJYNEW) )
!	QF=0.D0
	

	DO I1=1,JJX
	DO J1=1,JJYNEW
	  AAC=0.D0
	  AAD=0.D0
	  AAF=0.D0
	  DO K=1,JJ
	    AAC=AAC + H44(I1,K)*OS(JRIGHT)%op(K,J1)
	    AAD=AAD + H22(I1,K)*OS(JRIGHT)%op(K,J1)
!	    AAF=AAF + H24(I1,K)*OS(JRIGHT)%op(K,J1)
	  END DO
	  QC(I1,J1)=AAC
	  QD(I1,J1)=AAD
!	  QF(I1,J1)=AAF
	END DO
	END DO

	DEALLOCATE(H44)
	DEALLOCATE(H22)
!	DEALLOCATE(H24)

!* O^{+}([ACDOWN]O)
	JJ=IDIMB1HSIS(JLEFT)
	JJXNEW=IDIMB1HSIS(JLEFT)
	DO I1=1,JJXNEW
	DO J1=1,JJYNEW
	  AAC=0.D0
	  AAD=0.D0
	  AAF=0.D0
	  DO K=1,JJ
	    AAC=AAC + (OS(JLEFT)%op(K,I1))*QC(K,J1)
	    AAD=AAD + (OS(JLEFT)%op(K,I1))*QD(K,J1)
!	    AAF=AAF + (OS(JLEFT)%op(K,I1))*QF(K,J1)
	  END DO
	  CG44(J)%op(I1,J1)=AAC
	  CG22(J)%op(I1,J1)=AAD
!	  CG24(J)%op(I1,J1)=AAF
	END DO
	END DO
	DEALLOCATE(QC)
	DEALLOCATE(QD)
!	DEALLOCATE(QF)
	
! 	PRINT*,'%%%%%%%%%%%%     RESULT ROTATION   %%%%%%%%%'
! 	DO I1=1,IDIMB1CDOWNSIS(1,J)
! 	DO J1=1,IDIMB1CDOWNSIS(2,J)
! 	IF(CG44(J)%op(I1,J1)/=0.D0)PRINT*,J,I1,J1,CG44(J)%op(I1,J1)
! 	END DO
! 	END DO

     END IF ! Solo los bloques usados
1000	END DO		!#BLOQUES-CDOWNSIS


!----> Calculating the partition function of the system
	PRINT*,'BETA=',BETA,'EMIN=',EMIN
	AZN=0.D0
	DO J=1,NBL1HSIS
	IF((IBHUSED(J)).EQ.1)THEN
!	PRINT*,'Partition function',J
!    IF (ICHARB1SIS(J) .NE. LFIN)THEN
	DO M=1,IDIMB1HSIS(J)
	AZN=AZN+DEXP(-BETA*(WPHOS(J)%op3(M)-EMIN))
!       PRINT*,'AZ_LIMIT',J,M,WPHOS(J)%op3(M),EMIN,DEXP(-BETA*(WPHOS(J)%op3(M)-EMIN)),AZN
	END DO
!        ELSE
!        AZN=AZN+DEXP(-BETA*(WPHOS(J)%op3(1)-EMIN))
!       PRINT*,'AZ_LIMIT',J,1,WPHOS(J)%op3(1),EMIN,DEXP(-BETA*(WPHOS(J)%op3(1)-EMIN)),AZN
!        END IF
	END IF
	END DO
        PRINT*,'partition function AZ=',AZN

        L=0
        LL=0
        LLL=0

	DO 2000 J=1,NBL1CDOWNSIS
!	J=6
!	PRINT*,'J=',J
	JBTX=IB1CDOWNSISX(J)
	JBTY=IB1CDOWNSISY(J)
	IF(((IBHUSED(JBTX)).EQ.1).AND.((IBHUSED(JBTY)).EQ.1))THEN
!	PRINT*,'PASE--> J=',J,JBTX,JBTY
	NBX=IB1CDOWNSISX(J)
	NBY=IB1CDOWNSISY(J)
!	PRINT*,J,IB1CDOWNSISX(J),IB1CDOWNSISY(J)
!	PRINT*,IDIMB1CDOWNSIS(1,J),IDIMB1CDOWNSIS(2,J)
	DO I1=1,IDIMB1CDOWNSIS(1,J)
	DO J1=1,IDIMB1CDOWNSIS(2,J)
	ADEXP=DEXP(-BETA*(WPHOS(NBX)%op3(I1)-EMIN))+DEXP(-BETA*(WPHOS(NBY)%op3(J1)-EMIN))
	ADEXP=ADEXP/AZN
	CG24(J)%op(I1,J1)=(CG44(J)%op(I1,J1)*CG22(J)%op(I1,J1))
	AREX44=ADEXP*((CG44(J)%op(I1,J1))**2)
	CG44(J)%op(I1,J1)=AREX44
        IF(ABS(CG44(J)%op(I1,J1)) .GT. 1.D-12) THEN
        L=L+1
        END IF
	AREX22=ADEXP*((CG22(J)%op(I1,J1))**2)
	CG22(J)%op(I1,J1)=AREX22
        IF(ABS(CG22(J)%op(I1,J1)) .GT. 1.D-12) THEN
        LL=LL+1
        END IF
!	AREX24=ADEXP*((CG24(J)%op(I1,J1))**2)
	AREX24=ADEXP*CG24(J)%op(I1,J1)
	CG24(J)%op(I1,J1)=AREX24
        IF(ABS(CG24(J)%op(I1,J1)) .GT. 1.D-12) THEN
        LLL=LLL+1
        END IF
!	PRINT*,I1,J1,CG44(J)%op(I1,J1),CG22(J)%op(I1,J1),CG24(J)%op(I1,J1)
!	PRINT*,I1,J1,WPHOS(NBX)%op3(I1),WPHOS(NBY)%op3(J1),ADEXP
	END DO
	END DO
	
	END IF
2000	END DO		!#BLOQUES-CDOWNSIS
	
        IF(ALLOCATED(S)) DEALLOCATE(S)
        ALLOCATE(S(L,3))
        S=0.D0
        IF(ALLOCATED(SS)) DEALLOCATE(SS)
        ALLOCATE(SS(LL,3))
        SS=0.D0
        IF(ALLOCATED(SSS)) DEALLOCATE(SSS)
        ALLOCATE(SSS(LLL,3))
        SSS=0.D0

        L=0
        LL=0
        LLL=0

        DO 6000 J=1,NBL1CDOWNSIS
        JBTX=IB1CDOWNSISX(J)
        JBTY=IB1CDOWNSISY(J)
        IF(((IBHUSED(JBTX)).EQ.1).AND.((IBHUSED(JBTY)).EQ.1))THEN
        NBX=IB1CDOWNSISX(J)
        NBY=IB1CDOWNSISY(J)
        DO I1=1,IDIMB1CDOWNSIS(1,J)
        DO J1=1,IDIMB1CDOWNSIS(2,J)
        IF(ABS(CG44(J)%op(I1,J1)) .GT. 1.D-12) THEN
        S(L+1,1)=J
        S(L+1,2)=I1
        S(L+1,3)=J1
        L=L+1
        END IF
        IF(ABS(CG22(J)%op(I1,J1)) .GT. 1.D-12) THEN
        SS(LL+1,1)=J
        SS(LL+1,2)=I1
        SS(LL+1,3)=J1
        LL=LL+1
        END IF
        IF(ABS(CG24(J)%op(I1,J1)) .GT. 1.D-12) THEN
        SSS(LLL+1,1)=J
        SSS(LLL+1,2)=I1
        SSS(LLL+1,3)=J1
        LLL=LLL+1
        END IF
        END DO
        END DO
        END IF
6000    END DO

	PRINT*,'O.K.'
	
	
	END SUBROUTINE MEOCDOWN	
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE LEE
	USE PRINCIPAL
	USE AMODESP1
	INTEGER :: J,IIAA,IEA,IEB


	DO 10 J=1,NBL1H
	 ALLOCATE(HB1(J)%op(IDIMB1H(J),IDIMB1H(J)))
	 HB1(J)%op=0.d0
	 READ(9100+J) IIAA
	 READ(9100+J) HB1(J)%op
	 REWIND(9100+J)
	 CLOSE(9100+J)
10	CONTINUE

	DO 15 J=1,NBL1H
	 ALLOCATE(ANN(J)%op(IDIMB1H(J),IDIMB1H(J)))
	 ANN(J)%op=0.D0
	 READ(9700+J) IIAA
	 READ(9700+J) ANN(J)%op
	 REWIND(9700+J)
	 CLOSE(9700+J)
15	CONTINUE

	DO 20 J=1,NBL1CUP
	 ALLOCATE(ACUP(J)%op(IDIMB1CUP(1,J),IDIMB1CUP(2,J)))
	 ACUP(J)%op=0.D0
	 READ(9300+J) IEA,IEB
	 READ(9300+J) ACUP(J)%op
	 REWIND(9300+J)
	 CLOSE(9300+J)
20	CONTINUE


	DO 30 J=1,NBL1CDOWN
	 ALLOCATE(ACDOWN(J)%op(IDIMB1CDOWN(1,J),IDIMB1CDOWN(2,J)))
	 ACDOWN(J)%op=0.D0
	 READ(9500+J) IEA,IEB
	 READ(9500+J) ACDOWN(J)%op
	 REWIND(9500+J)
	 CLOSE(9500+J)
30	CONTINUE
	
	DO 40 J=1,NBL1CUP
	 ALLOCATE(BCUP(J)%op(IDIMB1CUP(1,J),IDIMB1CUP(2,J)))
	 BCUP(J)%op=0.D0
	 READ(9600+J) IEA,IEB
	 READ(9600+J) BCUP(J)%op
	 REWIND(9600+J)
	 CLOSE(9600+J)
40	CONTINUE


	DO 50 J=1,NBL1CDOWN
	 ALLOCATE(BCDOWN(J)%op(IDIMB1CDOWN(1,J),IDIMB1CDOWN(2,J)))
	 BCDOWN(J)%op=0.D0
	 READ(9800+J) IEA,IEB
	 READ(9800+J) BCDOWN(J)%op
	 REWIND(9800+J)
	 CLOSE(9800+J)
50	CONTINUE


	END SUBROUTINE LEE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE TRANSFNEW
	USE PRINCIPAL
	INTEGER :: J

	NBL1H=NBL1HSIS
	DO J=1,NBL1H
	   ICHARB1(J)=ICHARB1SIS(J)
	   IU1B1(J)=IU1B1SIS(J)
	   IDIMB1H(J)=IDIMB1HSIS(J)
	END DO

	NBL1CUP=NBL1CUPSIS
	DO J=1,NBL1CUP
	  ICHARB1XCUP(J)=ICHARB1XCUPSIS(J)
	  ICHARB1YCUP(J)=ICHARB1YCUPSIS(J)
	  IU1B1XCUP(J)=IU1B1XCUPSIS(J)
	  IU1B1YCUP(J)=IU1B1YCUPSIS(J)  
	  IDIMB1CUP(1,J)=IDIMB1CUPSIS(1,J)
	  IDIMB1CUP(2,J)=IDIMB1CUPSIS(2,J)
	END DO

	NBL1CDOWN=NBL1CDOWNSIS
	DO J=1,NBL1CDOWN
	  ICHARB1XCDOWN(J)=ICHARB1XCDOWNSIS(J)
	  ICHARB1YCDOWN(J)=ICHARB1YCDOWNSIS(J)
	  IU1B1XCDOWN(J)=IU1B1XCDOWNSIS(J)
	  IU1B1YCDOWN(J)=IU1B1YCDOWNSIS(J)  
	  IDIMB1CDOWN(1,J)=IDIMB1CDOWNSIS(1,J)
	  IDIMB1CDOWN(2,J)=IDIMB1CDOWNSIS(2,J)
	END DO
	
	
	END SUBROUTINE TRANSFNEW
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE DEALLOC_IGSIS
	USE PRINCIPAL
	USE AMODIGSIS
	INTEGER :: J,IB1,IB2,IU1,ICHART

	DO J=1,NBL1HSIS

	DO 200 IB1=1,NBL1H
	DO 200 IB2=1,NBL2H

	IU1=IU1B1(IB1)+IU1B2(IB2)
	ICHART=ICHARB1(IB1)+ICHARB2(IB2)

	IF((IU1==IU1B1SIS(J)).AND.(ICHART==ICHARB1SIS(J))) THEN	!100

	 DEALLOCATE( IGSIS(IB1,IB2)%op4  )
100	END IF
200	CONTINUE
	END DO

	END SUBROUTINE DEALLOC_IGSIS
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE TCUPES1
	USE PRINCIPAL
	USE AMODESP1
	USE AMODESP2
	USE AMODIGSIS
	USE AMODO
	USE AMODOLO
	USE AMODHSIS 

	IMPLICIT NONE
	INTEGER :: J,JJX,JJY,IU1XX,IU1YY,ICHARXX,ICHARYY,IFIND,IU1T
	INTEGER :: ICHART,JRIGHT,JLEFT,IB1,ISIG,IB2,IZ,IIB1,JIB1
	INTEGER :: IIB2,JIB2,I1,I2,IX11,J1,J2,IY11,JJ,JJYNEW,J11
	INTEGER :: K,JJXNEW,I11,NS1,IEA,IEB
	REAL(KIND=KIND(1.0D0)) :: AAE1,AAE2,AAC
	REAL(KIND=KIND(1.0D0)), DIMENSION(:,:), ALLOCATABLE :: HBS,HBD

	DO 1000 IZ=1,LAT-1
	
	
	DO 900 J=1,NBL1CUPSIS
	
	ALLOCATE(CUPS(J)%op(IDIMB1CUPSIS(1,J),IDIMB1CUPSIS(2,J)))
	CUPS(J)%op=0.D0
	ALLOCATE(CUPD(J)%op(IDIMB1CUPSIS(1,J),IDIMB1CUPSIS(2,J)))
	CUPD(J)%op=0.D0
!CCCCC ---->     Cup x 1
	DO 500 IB1=1,NBL1CUP
!	ISITEM
	NS1=10000+1000*IZ+IB1
	READ(NS1)IEA,IEB
	ALLOCATE(HBS(IEA,IEB))
	READ(NS1)HBS
	REWIND(NS1)
	CLOSE(NS1)
	
!	IDUTEM
	NS1=20000+1000*IZ+IB1
	READ(NS1)IEA,IEB
	ALLOCATE(HBD(IEA,IEB))
	READ(NS1)HBD
	REWIND(NS1)
	CLOSE(NS1)
	
	DO 400 IB2=1,NBL2H

	IU1XX=IU1B1XCUP(IB1)+ IU1B2(IB2)
	ICHARXX=ICHARB1XCUP(IB1)+ ICHARB2(IB2)

	IF((IU1XX==IU1B1XCUPSIS(J)).AND.(ICHARXX==ICHARB1XCUPSIS(J))) THEN !250

	    IIB1=IB1CUPX(IB1)
	    JIB1=IB1CUPY(IB1)

	  DO 200 I1=1,IDIMB1CUP(1,IB1)
	  DO 200 I2=1,IDIMB2H(IB2)

	    IX11=IGSIS(IIB1,IB2)%op4(I1,I2)

	    DO 100 J1=1,IDIMB1CUP(2,IB1)
	    DO 100 J2=1,IDIMB2H(IB2)

	    IY11=IGSIS(JIB1,IB2)%op4(J1,J2)
!* CUPx1  
	AAE1=0.D0
	AAE2=0.D0
	IF(J2==I2)AAE1=HBS(I1,J1)
	IF(J2==I2)AAE2=HBD(I1,J1)

	CUPS(J)%op(IX11,IY11)=CUPS(J)%op(IX11,IY11)+AAE1
	CUPD(J)%op(IX11,IY11)=CUPD(J)%op(IX11,IY11)+AAE2
	    
100	    CONTINUE
200	  CONTINUE
250	END IF
400	CONTINUE
	DEALLOCATE(HBS)
	DEALLOCATE(HBD)
500	CONTINUE
900	END DO
	

	DO J=1,NBL1CUPSIS
!	ISITEM
	NS1=10000+1000*IZ+J
	WRITE(NS1)IDIMB1CUPSIS(1,J),IDIMB1CUPSIS(2,J)
	WRITE(NS1)CUPS(J)%op
	REWIND(NS1)
	CLOSE(NS1)
	DEALLOCATE(CUPS(J)%op)
	
!	IDUTEM
	NS1=20000+1000*IZ+J
	WRITE(NS1)IDIMB1CUPSIS(1,J),IDIMB1CUPSIS(2,J)
	WRITE(NS1)CUPD(J)%op
	REWIND(NS1)
	CLOSE(NS1)
	DEALLOCATE(CUPD(J)%op)
		
	END DO

1000	END DO	!sitios

	END SUBROUTINE TCUPES1
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE TCDOWNES1
	USE PRINCIPAL
	USE AMODESP1
	USE AMODESP2
	USE AMODIGSIS
	USE AMODO
	USE AMODOLO
	USE AMODHSIS 

	IMPLICIT NONE
	INTEGER :: J,JJX,JJY,IU1XX,IU1YY,ICHARXX,ICHARYY,IFIND,IU1T
	INTEGER :: ICHART,JRIGHT,JLEFT,IB1,ISIG,IB2,IZ,IIB1,JIB1
	INTEGER :: IIB2,JIB2,I1,I2,IX11,J1,J2,IY11,JJ,JJYNEW,J11
	INTEGER :: K,JJXNEW,I11,NS1,IEA,IEB
	REAL(KIND=KIND(1.0D0)) :: AAE1,AAE2,AAC
	REAL(KIND=KIND(1.0D0)), DIMENSION(:,:), ALLOCATABLE :: HBS,HBD

	DO 1000 IZ=1,LAT-1
	
	
	DO 900 J=1,NBL1CDOWNSIS
	
	ALLOCATE(CDOWNS(J)%op(IDIMB1CDOWNSIS(1,J),IDIMB1CDOWNSIS(2,J)))
	CDOWNS(J)%op=0.D0
	ALLOCATE(CDOWND(J)%op(IDIMB1CDOWNSIS(1,J),IDIMB1CDOWNSIS(2,J)))
	CDOWND(J)%op=0.D0
!CCCCC ---->     CDOWN x 1
	DO 500 IB1=1,NBL1CDOWN
!	ISITEM
	NS1=30000+1000*IZ+IB1
	READ(NS1)IEA,IEB
	ALLOCATE(HBS(IEA,IEB))
	READ(NS1)HBS
	REWIND(NS1)
	CLOSE(NS1)
	
!	IDUTEM
	NS1=40000+1000*IZ+IB1
	READ(NS1)IEA,IEB
	ALLOCATE(HBD(IEA,IEB))
	READ(NS1)HBD
	REWIND(NS1)
	CLOSE(NS1)
	
	DO 400 IB2=1,NBL2H

	IU1XX=IU1B1XCDOWN(IB1)+ IU1B2(IB2)
	ICHARXX=ICHARB1XCDOWN(IB1)+ ICHARB2(IB2)

	IF((IU1XX==IU1B1XCDOWNSIS(J)).AND.(ICHARXX==ICHARB1XCDOWNSIS(J))) THEN !250

	    IIB1=IB1CDOWNX(IB1)
	    JIB1=IB1CDOWNY(IB1)

	  DO 200 I1=1,IDIMB1CDOWN(1,IB1)
	  DO 200 I2=1,IDIMB2H(IB2)

	    IX11=IGSIS(IIB1,IB2)%op4(I1,I2)

	    DO 100 J1=1,IDIMB1CDOWN(2,IB1)
	    DO 100 J2=1,IDIMB2H(IB2)

	    IY11=IGSIS(JIB1,IB2)%op4(J1,J2)
!* CDOWNx1  
	AAE1=0.D0
	AAE2=0.D0
	IF(J2==I2)AAE1=HBS(I1,J1)
	IF(J2==I2)AAE2=HBD(I1,J1)

	CDOWNS(J)%op(IX11,IY11)=CDOWNS(J)%op(IX11,IY11)+AAE1
	CDOWND(J)%op(IX11,IY11)=CDOWND(J)%op(IX11,IY11)+AAE2
	    
100	    CONTINUE
200	  CONTINUE
250	END IF
400	CONTINUE
	DEALLOCATE(HBS)
	DEALLOCATE(HBD)
500	CONTINUE
900	END DO
	

	DO J=1,NBL1CDOWNSIS
!	ISITEM
	NS1=30000+1000*IZ+J
	WRITE(NS1)IDIMB1CDOWNSIS(1,J),IDIMB1CDOWNSIS(2,J)
	WRITE(NS1)CDOWNS(J)%op
	REWIND(NS1)
	CLOSE(NS1)
	DEALLOCATE(CDOWNS(J)%op)
	
!	IDUTEM
	NS1=40000+1000*IZ+J
	WRITE(NS1)IDIMB1CDOWNSIS(1,J),IDIMB1CDOWNSIS(2,J)
	WRITE(NS1)CDOWND(J)%op
	REWIND(NS1)
	CLOSE(NS1)
	DEALLOCATE(CDOWND(J)%op)
		
	END DO

1000	END DO	!sitios

	END SUBROUTINE TCDOWNES1	
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE TCUPSISQ
	USE PRINCIPAL
	USE AMODESP1
	USE AMODESP2
	USE AMODIGSIS
	USE AMODO
	USE AMODOLO
	USE AMODHSIS

	IMPLICIT NONE
	INTEGER :: J,JJX,JJY,IU1XX,IU1YY,ICHARXX,ICHARYY,IFIND,IU1T
	INTEGER :: ICHART,JRIGHT,JLEFT,IB1,ISIG,IB2
	INTEGER :: IIB1,JIB1,I1,I2,IX11,J1,J2,IY11,JJ,JJYNEW,J11
	INTEGER :: K,JJXNEW,I11,NS1
	REAL(KIND=KIND(1.0D0)) :: AAE,AAC

	REAL(KIND=KIND(1.0D0)), DIMENSION(:,:), ALLOCATABLE :: HBB,Q,HBS,HBD

!*-> encontrar setores de O^{dagga} e O que conectam com ACUP

	DO 1000 J=1,NBL1CUPSIS

	ALLOCATE(HBB(IDIMB1CUPSIS(1,J),IDIMB1CUPSIS(2,J)) )
	HBB=0.D0

!CCCCC ---->     Cup x1
	DO 900 IB1=1,NBL1CUP

!	ISIG=(-1)**(ICHARB1(IB1))

	DO 900 IB2=1,NBL2H

	IU1XX=IU1B1XCUP(IB1)+IU1B2(IB2)
	ICHARXX=ICHARB1XCUP(IB1)+ICHARB2(IB2)

	IF((IU1XX==IU1B1XCUPSIS(J)).AND.(ICHARXX==ICHARB1XCUPSIS(J))) THEN !850

	    IIB1=IB1CUPX(IB1)
	    JIB1=IB1CUPY(IB1)

	  DO 800 I1=1,IDIMB1CUP(1,IB1)
	  DO 800 I2=1,IDIMB2H(IB2)

	    IX11=IGSIS(IIB1,IB2)%op4(I1,I2)

	    DO 700 J1=1,IDIMB1CUP(2,IB1)
	    DO 700 J2=1,IDIMB2H(IB2)

	    IY11=IGSIS(JIB1,IB2)%op4(J1,J2)
!* CUPx1
	AAE=0.D0
	IF(J2==I2)AAE=BCUP(IB1)%op(I1,J1)

	    HBB(IX11,IY11)=HBB(IX11,IY11)+AAE


700	    CONTINUE
800	  CONTINUE
850	END IF
900	CONTINUE

	WRITE(9600+J) IDIMB1CUPSIS(1,J),IDIMB1CUPSIS(2,J)
	WRITE(9600+J) HBB
	REWIND(9600+J)
	CLOSE(9600+J)
	DEALLOCATE(HBB)

1000	END DO

	END SUBROUTINE TCUPSISQ
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE TCUPSIS
	USE PRINCIPAL
	USE AMODESP1
	USE AMODESP2
	USE AMODIGSIS
	USE AMODO
	USE AMODOLO
	USE AMODHSIS 

	IMPLICIT NONE
	INTEGER :: J,JJX,JJY,IU1XX,IU1YY,ICHARXX,ICHARYY,IFIND,IU1T
	INTEGER :: ICHART,JRIGHT,JLEFT,IB1,ISIG,IB2
	INTEGER :: IIB2,JIB2,I1,I2,IX11,J1,J2,IY11,JJ,JJYNEW,J11
	INTEGER :: K,JJXNEW,I11,NS1
	REAL(KIND=KIND(1.0D0)) :: AAE,AAC

	REAL(KIND=KIND(1.0D0)), DIMENSION(:,:), ALLOCATABLE :: HBB,Q,HBS,HBD

!*-> encontrar setores de O^{dagga} e O que conectam com ACUP

	DO 1000 J=1,NBL1CUPSIS

	ALLOCATE(HBB(IDIMB1CUPSIS(1,J),IDIMB1CUPSIS(2,J)) )
	HBB=0.D0
	ALLOCATE(HBS(IDIMB1CUPSIS(1,J),IDIMB1CUPSIS(2,J)) )
	HBS=0.D0
	ALLOCATE(HBD(IDIMB1CUPSIS(1,J),IDIMB1CUPSIS(2,J)) )
	HBD=0.D0
	
!CCCCC ---->     1 x Cup2
	DO 900 IB1=1,NBL1H

	ISIG=(-1)**(ICHARB1(IB1))

	DO 900 IB2=1,NBL2CUP

	IU1XX=IU1B1(IB1)+IU1B2XCUP(IB2)
	ICHARXX=ICHARB1(IB1)+ICHARB2XCUP(IB2)

	IF((IU1XX==IU1B1XCUPSIS(J)).AND.(ICHARXX==ICHARB1XCUPSIS(J))) THEN !850

	    IIB2=IB2CUPX(IB2)
	    JIB2=IB2CUPY(IB2)

	  DO 800 I1=1,IDIMB1H(IB1)
	  DO 800 I2=1,IDIMB2CUP(1,IB2)

	    IX11=IGSIS(IB1,IIB2)%op4(I1,I2)

	    DO 700 J1=1,IDIMB1H(IB1)
	    DO 700 J2=1,IDIMB2CUP(2,IB2)

	    IY11=IGSIS(IB1,JIB2)%op4(J1,J2)
!* 1xCUP  
	AAE=0.D0
	IF(J1==I1)AAE=ISIG*CUPL(IB2)%op(I2,J2)
	      
	    HBB(IX11,IY11)=HBB(IX11,IY11)+AAE
	    
	    
	IF(IIDUPLA2(IB2,I2,J2).EQ.1)THEN
	    HBD(IX11,IY11)=HBD(IX11,IY11)+AAE
	END IF
	
	IF(IISIMPLE2(IB2,I2,J2).EQ.1)THEN
	    HBS(IX11,IY11)=HBS(IX11,IY11)+AAE
	END IF	
700	    CONTINUE
800	  CONTINUE
850	END IF
900	CONTINUE

	WRITE(9300+J) IDIMB1CUPSIS(1,J),IDIMB1CUPSIS(2,J)
	WRITE(9300+J) HBB
	REWIND(9300+J)
	CLOSE(9300+J)
	DEALLOCATE(HBB)
	
!	ISITEM
	NS1=10000+1000*LAT+J
	WRITE(NS1)IDIMB1CUPSIS(1,J),IDIMB1CUPSIS(2,J)
	WRITE(NS1)HBS
	REWIND(NS1)
	CLOSE(NS1)
	DEALLOCATE(HBS)
	
!	IDUTEM
	NS1=20000+1000*LAT+J
	WRITE(NS1)IDIMB1CUPSIS(1,J),IDIMB1CUPSIS(2,J)
	WRITE(NS1)HBD
	REWIND(NS1)
	CLOSE(NS1)
	DEALLOCATE(HBD)

1000	END DO

	END SUBROUTINE TCUPSIS
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE TCDOWNSISQ
	USE PRINCIPAL
	USE AMODESP1
	USE AMODESP2
	USE AMODIGSIS
	USE AMODO
	USE AMODOLO

	IMPLICIT NONE
	INTEGER :: J,JJX,JJY,IU1XX,IU1YY,ICHARXX,ICHARYY,IFIND,IU1T
	INTEGER :: ICHART,JRIGHT,JLEFT,IB1,ISIG,IB2,IIB1,JIB1,I1,I2
	INTEGER :: IX11,J1,J2,IY11,JJ,JJYNEW,J11,K,JJXNEW,I11,NS1
	REAL(KIND=KIND(1.0D0)) ::AAE,AAC


	REAL(KIND=KIND(1.0D0)), DIMENSION(:,:), ALLOCATABLE :: HBB,Q,HBS,HBD

!*-> encontrar setores de O^{dagga} e O que conectam com ACDOWN
	DO 1000 J=1,NBL1CDOWNSIS

	ALLOCATE(HBB(IDIMB1CDOWNSIS(1,J),IDIMB1CDOWNSIS(2,J)) )
	HBB=0.D0
!CCCCC ---->     Cdown x 1
	DO 500 IB1=1,NBL1CDOWN

!	ISIG=(-1)**(ICHARB1(IB1))

	DO 500 IB2=1,NBL2H

	IU1XX=IU1B1XCDOWN(IB1)+IU1B2(IB2)
	ICHARXX=ICHARB1XCDOWN(IB1)+ICHARB2(IB2)

	IF((IU1XX==IU1B1XCDOWNSIS(J)).AND.(ICHARXX==ICHARB1XCDOWNSIS(J))) THEN !450

	    IIB1=IB1CDOWNX(IB1)
	    JIB1=IB1CDOWNY(IB1)

	  DO 200 I1=1,IDIMB1CDOWN(1,IB1)
	  DO 200 I2=1,IDIMB2H(IB2)

	    IX11=IGSIS(IIB1,IB2)%op4(I1,I2)

	    DO 100 J1=1,IDIMB1CDOWN(2,IB1)
	    DO 100 J2=1,IDIMB2H(IB2)

	    IY11=IGSIS(JIB1,IB2)%op4(J1,J2)
!* CDOWNx1
	AAE=0.D0
	IF(J2==I2)AAE=BCDOWN(IB1)%op(I1,J1)

	    HBB(IX11,IY11)=HBB(IX11,IY11)+AAE

100	    CONTINUE
200	  CONTINUE
450	END IF
500	CONTINUE

	WRITE(9800+J) IDIMB1CDOWNSIS(1,J),IDIMB1CDOWNSIS(2,J)
	WRITE(9800+J) HBB
	REWIND(9800+J)
	CLOSE(9800+J)
	DEALLOCATE(HBB)


1000	END DO

	END SUBROUTINE TCDOWNSISQ
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE TCDOWNSIS
	USE PRINCIPAL
	USE AMODESP1
	USE AMODESP2
	USE AMODIGSIS
	USE AMODO
	USE AMODOLO

	IMPLICIT NONE
	INTEGER :: J,JJX,JJY,IU1XX,IU1YY,ICHARXX,ICHARYY,IFIND,IU1T
	INTEGER :: ICHART,JRIGHT,JLEFT,IB1,ISIG,IB2,IIB2,JIB2,I1,I2
	INTEGER :: IX11,J1,J2,IY11,JJ,JJYNEW,J11,K,JJXNEW,I11,NS1
	REAL(KIND=KIND(1.0D0)) ::AAE,AAC


	REAL(KIND=KIND(1.0D0)), DIMENSION(:,:), ALLOCATABLE :: HBB,Q,HBS,HBD

!*-> encontrar setores de O^{dagga} e O que conectam com ACDOWN
	DO 1000 J=1,NBL1CDOWNSIS

	ALLOCATE(HBB(IDIMB1CDOWNSIS(1,J),IDIMB1CDOWNSIS(2,J)) )
	HBB=0.D0
	ALLOCATE(HBS(IDIMB1CDOWNSIS(1,J),IDIMB1CDOWNSIS(2,J)) )
	HBS=0.D0
	ALLOCATE(HBD(IDIMB1CDOWNSIS(1,J),IDIMB1CDOWNSIS(2,J)) )
	HBD=0.D0
!CCCCC ---->     1 x Cdown2
	DO 500 IB1=1,NBL1H

	ISIG=(-1)**(ICHARB1(IB1))

	DO 500 IB2=1,NBL2CDOWN

	IU1XX=IU1B1(IB1)+IU1B2XCDOWN(IB2)
	ICHARXX=ICHARB1(IB1)+ICHARB2XCDOWN(IB2)

	IF((IU1XX==IU1B1XCDOWNSIS(J)).AND.(ICHARXX==ICHARB1XCDOWNSIS(J))) THEN !450

	    IIB2=IB2CDOWNX(IB2)
	    JIB2=IB2CDOWNY(IB2)

	  DO 200 I1=1,IDIMB1H(IB1)
	  DO 200 I2=1,IDIMB2CDOWN(1,IB2)

	    IX11=IGSIS(IB1,IIB2)%op4(I1,I2)

	    DO 100 J1=1,IDIMB1H(IB1)
	    DO 100 J2=1,IDIMB2CDOWN(2,IB2)

	    IY11=IGSIS(IB1,JIB2)%op4(J1,J2)
!* 1xCDOWN
	AAE=0.D0
	IF(J1==I1)AAE=ISIG*CDOWNL(IB2)%op(I2,J2)

	    HBB(IX11,IY11)=HBB(IX11,IY11)+AAE
	    
	IF(FFDUPLA2(IB2,I2,J2).EQ.1)THEN
	    HBD(IX11,IY11)=HBD(IX11,IY11)+AAE
	END IF
	
	IF(FFSIMPLE2(IB2,I2,J2).EQ.1)THEN
	    HBS(IX11,IY11)=HBS(IX11,IY11)+AAE
	END IF		    

100	    CONTINUE
200	  CONTINUE
450	END IF
500	CONTINUE

	WRITE(9500+J) IDIMB1CDOWNSIS(1,J),IDIMB1CDOWNSIS(2,J)
	WRITE(9500+J) HBB
	REWIND(9500+J)
	CLOSE(9500+J)
	DEALLOCATE(HBB)
	
!	ISITEM
	NS1=30000+1000*LAT+J
	WRITE(NS1)IDIMB1CDOWNSIS(1,J),IDIMB1CDOWNSIS(2,J)
	WRITE(NS1)HBS
	REWIND(NS1)
	CLOSE(NS1)
	DEALLOCATE(HBS)
	
!	IDUTEM
	NS1=40000+1000*LAT+J
	WRITE(NS1)IDIMB1CDOWNSIS(1,J),IDIMB1CDOWNSIS(2,J)
	WRITE(NS1)HBD
	REWIND(NS1)
	CLOSE(NS1)
	DEALLOCATE(HBD)	
	
	
1000	END DO

	END SUBROUTINE TCDOWNSIS
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE THSIS
	USE PRINCIPAL
	USE AMODESP1
	USE AMODESP2
	USE AMODOLO
	USE AMODIGSIS
	USE AMODO
	USE AMODHSIS

	IMPLICIT NONE
	INTEGER :: J,IU1XSIS,ICHARXSIS,JJ,IB1,IB2,IU1,ICHART,I1,I2
	INTEGER :: II,J1,J2,ISIG,IU1XX,ICHARXX,IIB1,IIB2,JIB1,JIB2
	INTEGER :: IX11,IY12,IY11,IX12,JJNEW,J11,K,I11,M,INFO
        INTEGER :: LWORK,IL,IU,MM,LIWORK
	REAL(KIND=KIND(1.0D0)) :: AAE1,AAE2,AAE3,AAE,AAC,AAF
        REAL(KIND=KIND(1.0D0)) :: DLAMCH, VL, VU


	CHARACTER     JOBZ,UPLO
        INTEGER, DIMENSION(1) :: DUMMY2
        INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
        INTEGER, DIMENSION(:), ALLOCATABLE :: ISUPPZ
        REAL(KIND=KIND(1.0D0)), DIMENSION(1) :: DUMMY
	REAL(KIND=KIND(1.0D0)), DIMENSION(:), ALLOCATABLE :: W
	REAL(KIND=KIND(1.0D0)), DIMENSION(:), ALLOCATABLE :: work
	REAL(KIND=KIND(1.0D0)), DIMENSION(:,:), ALLOCATABLE :: A
        REAL(KIND=KIND(1.0D0)), DIMENSION(:,:), ALLOCATABLE :: Z

!	WRITE(44,*)'#############---->LAT=',LAT    
	EMIN=0.D0		!Energía minima
    PRINT*,'LFIN=',LFIN
	
	DO 1000 J=1,NBL1HSIS
	IU1XSIS=IU1B1SIS(J)
	ICHARXSIS=ICHARB1SIS(J)
	JJ=IDIMB1HSIS(J)
	IBHUSED(J)=0 ! IBHUSED(J)=0 ONLY CONSIDERS BLOCKS THAT CONNECT TO THE GROUND STATE (CAN ONLY BE USED AT LOW TEMPERATURES AND HALF-FILLING, I.E. T=1/BETAI~0, AMUI=0 AND ACOPnI=-0.5D0*ACOPUI-AMUI); IBHUSED(J)=1 CONSIDERS ALL BLOCKS

        ALLOCATE(WPHOS(J)%op3(JJ))
	ALLOCATE(OS(J)%op(JJ,JJ))

	IF( ALLOCATED(A) )  DEALLOCATE(A) 
	ALLOCATE (A(JJ,JJ))
	A=0.D0
	
	
!CCCCC ---->     HB1x1    y   1x HB2
	DO 50 IB1=1,NBL1H
	DO 50 IB2=1,NBL2H

	IU1=IU1B1(IB1)+IU1B2(IB2)
	ICHART=ICHARB1(IB1)+ICHARB2(IB2)

	IF((IU1==IU1XSIS).AND.(ICHART==ICHARXSIS)) THEN	!40
	
	  DO 20 I1=1,IDIMB1H(IB1)
	  DO 20 I2=1,IDIMB2H(IB2)

	    II=IGSIS(IB1,IB2)%op4(I1,I2)

	  DO 10 J1=1,IDIMB1H(IB1)
	  DO 10 J2=1,IDIMB2H(IB2)

	    JJ=IGSIS(IB1,IB2)%op4(J1,J2)

	AAE1=0.D0
	IF(J2==I2)AAE1=HB1(IB1)%op(I1,J1)

	AAE2=0.D0
	IF(J1==I1)AAE2=X2(IB2)%op(I2,J2)

	AAE3=ACOPV*ANN(IB1)%op(I1,J1)*XLN(IB2)%op(I2,J2)

	    A(II,JJ)=A(II,JJ)+AAE1+AAE2+AAE3
10	  CONTINUE
20	  CONTINUE  
40	END IF
50	CONTINUE
!CCCCC ---->     Cup1^{+}xCup2
	DO 300 IB1=1,NBL1CUP

	ISIG=(-1)**(ICHARB1YCUP(IB1))

	DO 300 IB2=1,NBL2CUP

	IU1XX=IU1B1XCUP(IB1)+IU1B2YCUP(IB2)
	ICHARXX=ICHARB1XCUP(IB1)+ICHARB2YCUP(IB2)
	
	IF((IU1XX==IU1XSIS).AND.(ICHARXX==ICHARXSIS)) THEN   !250

	    IIB1=IB1CUPX(IB1)
	    IIB2=IB2CUPY(IB2)

	    JIB1=IB1CUPY(IB1)
	    JIB2=IB2CUPX(IB2)

	  DO 200 I1=1,IDIMB1CUP(1,IB1)
	  DO 200 I2=1,IDIMB2CUP(2,IB2)

	    IX11=IGSIS(IIB1,IIB2)%op4(I1,I2)
	    IY12=IGSIS(IIB1,IIB2)%op4(I1,I2)

	    DO 100 J1=1,IDIMB1CUP(2,IB1)
	    DO 100 J2=1,IDIMB2CUP(1,IB2)

	    IY11=IGSIS(JIB1,JIB2)%op4(J1,J2)
	    IX12=IGSIS(JIB1,JIB2)%op4(J1,J2)

	AAE=-ISIG*ACt(LBL)*ACUP(IB1)%op(I1,J1)*CUPL(IB2)%op(J2,I2)

	AAF=0.D0
	IF(LAT.EQ.LFIN)THEN
	AAF=-ISIG*ACt(LFIN)*BCUP(IB1)%op(I1,J1)*CUPL(IB2)%op(J2,I2)
	END IF

	      A(IX11,IY11)=A(IX11,IY11)+AAE+AAF
	      A(IX12,IY12)=A(IX12,IY12)+AAE+AAF
100	  CONTINUE  
200	  CONTINUE
250	END IF
300	CONTINUE
!CCCCC ---->     Cdown1^{+}xCdown2
	DO 600 IB1=1,NBL1CDOWN

	ISIG=(-1)**(ICHARB1YCDOWN(IB1))

	DO 600 IB2=1,NBL2CDOWN

	IU1XX=IU1B1XCDOWN(IB1)+IU1B2YCDOWN(IB2)
	ICHARXX=ICHARB1XCDOWN(IB1)+ICHARB2YCDOWN(IB2)

	IF((IU1XX==IU1XSIS).AND.(ICHARXX==ICHARXSIS)) THEN   !550

	    IIB1=IB1CDOWNX(IB1)
	    IIB2=IB2CDOWNY(IB2)

	    JIB1=IB1CDOWNY(IB1)
	    JIB2=IB2CDOWNX(IB2)

	  DO 500 I1=1,IDIMB1CDOWN(1,IB1)
	  DO 500 I2=1,IDIMB2CDOWN(2,IB2)
	    
	    IX11=IGSIS(IIB1,IIB2)%op4(I1,I2)
	    IY12=IGSIS(IIB1,IIB2)%op4(I1,I2)

	    DO 400 J1=1,IDIMB1CDOWN(2,IB1)
	    DO 400 J2=1,IDIMB2CDOWN(1,IB2)

	    IY11=IGSIS(JIB1,JIB2)%op4(J1,J2)
	    IX12=IGSIS(JIB1,JIB2)%op4(J1,J2)

	AAE=-ISIG*ACt(LBL)*ACDOWN(IB1)%op(I1,J1)*CDOWNL(IB2)%op(J2,I2)

	AAF=0.D0
	IF(LAT.EQ.LFIN)THEN
	AAF=-ISIG*ACt(LFIN)*BCDOWN(IB1)%op(I1,J1)*CDOWNL(IB2)%op(J2,I2)
	END IF

	      A(IX11,IY11)=A(IX11,IY11)+AAE+AAF
	      A(IX12,IY12)=A(IX12,IY12)+AAE+AAF
400	  CONTINUE
500	  CONTINUE
550	END IF
600	CONTINUE

	ALLOCATE(HSIS(1)%op(IDIMB1HSIS(J),IDIMB1HSIS(J)))
	HSIS(1)%op=0.D0
	
	DO I1=1,JJ
	DO I2=1,JJ
	  HSIS(1)%op(I1,I2)=A(I1,I2)
!	IF(A(I1,I2).NE.0.D0)print*,J,I1,I2,A(I1,I2)
	END DO
	END DO
	
	WRITE(9100+J)IDIMB1HSIS(J)
	WRITE(9100+J) HSIS(1)%op
	REWIND(9100+J)
	CLOSE(9100+J)

	DEALLOCATE(HSIS(1)%op)
	
	IF(MOD(LFIN,2).EQ.0)THEN
	 IF((ICHARXSIS.EQ.LFIN).AND.(IU1XSIS.EQ.0))IBHUSED(J)=1
	 IF((ICHARXSIS.EQ.(LFIN-1)).AND.(IU1XSIS.EQ.(-1)))IBHUSED(J)=1
	 IF((ICHARXSIS.EQ.(LFIN+1)).AND.(IU1XSIS.EQ.(1)))IBHUSED(J)=1
	 ! DOWN
	 IF(IDOWN.EQ.1)THEN
	 IF((ICHARXSIS.EQ.(LFIN-1)).AND.(IU1XSIS.EQ.(1)))IBHUSED(J)=1
	 IF((ICHARXSIS.EQ.(LFIN+1)).AND.(IU1XSIS.EQ.(-1)))IBHUSED(J)=1
	 END IF
	ELSE
	 IF((ICHARXSIS.EQ.LFIN).AND.(IU1XSIS.EQ.1))IBHUSED(J)=1
	 IF((ICHARXSIS.EQ.(LFIN-1)).AND.(IU1XSIS.EQ.(0)))IBHUSED(J)=1
	 IF((ICHARXSIS.EQ.(LFIN+1)).AND.(IU1XSIS.EQ.(2)))IBHUSED(J)=1
	 ! DOWN
	 IF(IDOWN.EQ.1)THEN
	 IF((ICHARXSIS.EQ.(LFIN-1)).AND.(IU1XSIS.EQ.(2)))IBHUSED(J)=1
	 IF((ICHARXSIS.EQ.(LFIN+1)).AND.(IU1XSIS.EQ.(0)))IBHUSED(J)=1
	 END IF
     IF((ICHARXSIS.EQ.LFIN).AND.(IU1XSIS.EQ.(-1)))IBHUSED(J)=1
     IF((ICHARXSIS.EQ.(LFIN-1)).AND.(IU1XSIS.EQ.(-2)))IBHUSED(J)=1
     IF((ICHARXSIS.EQ.(LFIN+1)).AND.(IU1XSIS.EQ.(0)))IBHUSED(J)=1
     ! DOWN
     IF(IDOWN.EQ.1)THEN
     IF((ICHARXSIS.EQ.(LFIN-1)).AND.(IU1XSIS.EQ.(0)))IBHUSED(J)=1
     IF((ICHARXSIS.EQ.(LFIN+1)).AND.(IU1XSIS.EQ.(-2)))IBHUSED(J)=1
     END IF
	END IF
	
	
!CCCCCCCCCCCCCCCCCCCCCCCC  DIAGONALIZO CADA BLOQUE   CCCCC
    IF(LAT.EQ.LFIN)THEN     !999
    IF((IBHUSED(J)).EQ.1)THEN
    
!    PRINT*,'===========================================>>'   
!700	write(44,*)'===========================================>>'    
!	PRINT*,'BLOCO=',J,'CHARGE=',ICHARXSIS,'SPIN=',IU1XSIS,'DIM=',IDIMB1HSIS(J)
!    write(44,*)'BLOCO=',J,'CHARGE=',ICHARXSIS,'SPIN=',IU1XSIS,'DIM=',IDIMB1HSIS(J)
    
	ALLOCATE (W(JJ))
	W=0.D0

        LWORK=-1
        LIWORK=-1
        IL=1
        IU=1
!        IF (ICHARXSIS .NE. LFIN)THEN
        IF( ALLOCATED(Z) )  DEALLOCATE(Z)
        ALLOCATE (Z(JJ,JJ))
        Z=0.D0
        IF( ALLOCATED(ISUPPZ) )  DEALLOCATE(ISUPPZ)
        ALLOCATE (ISUPPZ(2*JJ))
        ISUPPZ=0.D0
        CALL DSYEVR('V','A','U',JJ,A,JJ,VL,VU,IL,IU,2*DLAMCH('S'),M,W,Z,JJ,ISUPPZ,DUMMY,LWORK,DUMMY2,LIWORK,INFO)
!        ELSE
!        IF( ALLOCATED(Z) )  DEALLOCATE(Z)
!        ALLOCATE (Z(JJ,IU-IL+1))
!        Z=0.D0
!        IF( ALLOCATED(ISUPPZ) )  DEALLOCATE(ISUPPZ)
!        ALLOCATE (ISUPPZ(2*(IU-IL+1)))
!        ISUPPZ=0.D0
!	CALL DSYEVR('V','I','U',JJ,A,JJ,VL,VU,IL,IU,2*DLAMCH('S'),M,W,Z,JJ,ISUPPZ,DUMMY,LWORK,DUMMY2,LIWORK,INFO)
!        END IF
        LWORK=INT(DUMMY(1))
        ALLOCATE (WORK(LWORK))
        LIWORK=DUMMY2(1)
        ALLOCATE (IWORK(LIWORK))
!        PRINT*, 'LWORK=',LWORK
!        PRINT*, 'LIWORK=',LIWORK
!        PRINT*, 'Z SHAPE=',SHAPE(Z)
!        IF (ICHARXSIS .NE. LFIN)THEN
        CALL DSYEVR('V','A','U',JJ,A,JJ,VL,VU,IL,IU,2*DLAMCH('S'),M,W,Z,JJ,ISUPPZ,WORK,LWORK,IWORK,LIWORK,INFO)
!        ELSE
!        CALL DSYEVR('V','I','U',JJ,A,JJ,VL,VU,IL,IU,2*DLAMCH('S'),M,W,Z,JJ,ISUPPZ,WORK,LWORK,IWORK,LIWORK,INFO)
!        END IF

	DEALLOCATE (A,work)

	IF(INFO/=0) THEN
	  WRITE(334,*) 'INFO dif zero',info
	  STOP
	END IF

        DEALLOCATE (IWORK,ISUPPZ)

        MM=M
!        WRITE(*,'(/A,I10)')' The total number of eigenvalues found:', MM

! 	ALLOCATE(WPHOS(J)%op3(JJ))
! 	ALLOCATE(OS(J)%op(JJ,JJ))
	DO M=1,MM
	IF(EMIN>W(M))EMIN=W(M)
	  WPHOS(J)%op3(M)=W(M)
!	PRINT*,'BLOCO',J,'AUTOVALOR->',M,'<E>=',WPHOS(J)%op3(M)
!        write(44,*)'AUTOVALOR->',M,'<E>=',WPHOS(J)%op3(M)	  
!	PRINT*,'AUTOVETOR->',M  
!	 WRITE(45,*)ICHARB1SIS(J),WPHOS(J)%op3(M)
!	 WRITE(44,*)'AUTOVETOR->',M
	  DO I1=1,JJ
	    OS(J)%op(I1,M)=Z(I1,M)
!	  PRINT*, I1,OS(J)%op(I1,M) 
!	  WRITE(44,*)I1,OS(J)%op(I1,M) 
! !	OJO	  
! 	    OS(J)%op(I1,M)=0.D0
! 	    IF(I1.EQ.M)THEN
! 	      OS(J)%op(I1,M)=1.D0
! 	    END IF
	  END DO
!	  WRITE(44,*)'-------->'
	END DO
	

	DEALLOCATE (W,Z)

!	PRINT*,'BOLQUE=',J,'CHARGE=',ICHARXSIS,'SPIN=',IU1XSIS

! 	DO M=1,JJ
! 	PRINT*,'I=',M,'<E>',WPHOS(J)%op3(M)
!         write(12,*)'I=',M,'<E>',WPHOS(J)%op3(M)
! 	END DO
800 END IF
	ELSE 
	 DEALLOCATE(A)
	END IF
!	PRINT*,J,IBHUSED(J)
1000	END DO


!	PRINT*,'EMIN=',EMIN
	
	END SUBROUTINE THSIS
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE TNSIS
	USE PRINCIPAL
	USE AMODESP1
	USE AMODESP2
	USE AMODOLO
	USE AMODIGSIS
	USE AMODO
	USE AMODHSIS

	IMPLICIT NONE
	INTEGER :: J,IU1XSIS,ICHARXSIS,JJ,IB1,IB2,IU1,ICHART,I1,I2
	INTEGER :: II,J1,J2,ISIG,IU1XX,ICHARXX,IIB1,IIB2,JIB1,JIB2
	INTEGER :: IX11,IY12,IY11,IX12,JJNEW,J11,K,I11,M,INFO
	REAL(KIND=KIND(1.0D0)) :: AAE1,AAE2,AAE3,AAE,AAC


	REAL(KIND=KIND(1.0D0)), DIMENSION(:,:), ALLOCATABLE :: HBB

	DO 1000 J=1,NBL1HSIS
	IU1XSIS=IU1B1SIS(J)
	ICHARXSIS=ICHARB1SIS(J)
    ALLOCATE(HBB(IDIMB1HSIS(J),IDIMB1HSIS(J)) )
	HBB=0.D0

!CCCCC ---->     HB1x1    y   1x HB2
	DO 50 IB1=1,NBL1H
	DO 50 IB2=1,NBL2H

	IU1=IU1B1(IB1)+IU1B2(IB2)
	ICHART=ICHARB1(IB1)+ICHARB2(IB2)

	IF((IU1==IU1XSIS).AND.(ICHART==ICHARXSIS)) THEN	!40

	  DO 20 I1=1,IDIMB1H(IB1)
	  DO 20 I2=1,IDIMB2H(IB2)

	    II=IGSIS(IB1,IB2)%op4(I1,I2)

	  DO 10 J1=1,IDIMB1H(IB1)
	  DO 10 J2=1,IDIMB2H(IB2)

	    JJ=IGSIS(IB1,IB2)%op4(J1,J2)

	AAE=0.D0
	IF(J1==I1)AAE=XLN(IB2)%op(I2,J2)

	    HBB(II,JJ)=HBB(II,JJ)+AAE
10	  CONTINUE
20	  CONTINUE
40	END IF
50	CONTINUE


	WRITE(9700+J)IDIMB1HSIS(J)
	WRITE(9700+J) HBB
	REWIND(9700+J)
	CLOSE(9700+J)

	DEALLOCATE(HBB)


1000	END DO

	END SUBROUTINE TNSIS
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE CONECTIONS
	USE PRINCIPAL
	IMPLICIT NONE 
	INTEGER ::J,IU1XX,ICHARXX,IU1YY,ICHARYY,IU1T,ICHART,IFIND


!*---------------------------- CUP ---------------------------------------
!* ESPACIO 1
	DO J=1,NBL1CUP
          IU1XX=IU1B1XCUP(J)
          ICHARXX=ICHARB1XCUP(J)
          IU1YY=IU1B1YCUP(J)
          ICHARYY=ICHARB1YCUP(J)
          DO IFIND=1,NBL1H
            IU1T=IU1B1(IFIND)
            ICHART=ICHARB1(IFIND)
            IF(IU1T.EQ.IU1XX.AND.ICHART.EQ.ICHARXX) IB1CUPX(J)=IFIND
	    IF(IU1T.EQ.IU1YY.AND.ICHART.EQ.ICHARYY) IB1CUPY(J)=IFIND
          END DO
	END DO


!* ESPACIO 2
	DO J=1,NBL2CUP
          IU1XX=IU1B2XCUP(J)
          ICHARXX=ICHARB2XCUP(J)
          IU1YY=IU1B2YCUP(J)
          ICHARYY=ICHARB2YCUP(J)
          DO IFIND=1,NBL2H
            IU1T=IU1B2(IFIND)
            ICHART=ICHARB2(IFIND)
            IF(IU1T.EQ.IU1XX.AND.ICHART.EQ.ICHARXX) IB2CUPX(J)=IFIND
	    IF(IU1T.EQ.IU1YY.AND.ICHART.EQ.ICHARYY) IB2CUPY(J)=IFIND
          END DO
	END DO
	
!*---------------------------- CDOWN --------------------------------------
!* ESPACIO 1
	DO J=1,NBL1CDOWN
          IU1XX=IU1B1XCDOWN(J)
          ICHARXX=ICHARB1XCDOWN(J)
          IU1YY=IU1B1YCDOWN(J)
          ICHARYY=ICHARB1YCDOWN(J)
          DO IFIND=1,NBL1H
            IU1T=IU1B1(IFIND)
            ICHART=ICHARB1(IFIND)
            IF(IU1T.EQ.IU1XX.AND.ICHART.EQ.ICHARXX) IB1CDOWNX(J)=IFIND
	    IF(IU1T.EQ.IU1YY.AND.ICHART.EQ.ICHARYY) IB1CDOWNY(J)=IFIND
          END DO
	END DO

!* ESPACIO 2
	DO J=1,NBL2CDOWN
          IU1XX=IU1B2XCDOWN(J)
          ICHARXX=ICHARB2XCDOWN(J)
          IU1YY=IU1B2YCDOWN(J)
          ICHARYY=ICHARB2YCDOWN(J)
          DO IFIND=1,NBL2H
            IU1T=IU1B2(IFIND)
            ICHART=ICHARB2(IFIND)
            IF(IU1T.EQ.IU1XX.AND.ICHART.EQ.ICHARXX) IB2CDOWNX(J)=IFIND
	    IF(IU1T.EQ.IU1YY.AND.ICHART.EQ.ICHARYY) IB2CDOWNY(J)=IFIND
          END DO
	END DO

!*-------->

	END SUBROUTINE CONECTIONS
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE ESTRUCTURA
	USE PRINCIPAL
	USE AMODIGSIS
	
	IMPLICIT NONE
	INTEGER :: I1,I2,IU1,ICHART,J,I,IIA,IIB,IIC,IU1A,ICHARA,IB1,IB2
	INTEGER :: ICONTSEC,IU1XX,ICHARXX,IU1YY,ICHARYY,IFIND,IU1T


!CCCCCCCCCCCCCCCCCCCCCCC  -----> SISTEMA
	ICONT=0
	DO 10 I1=1,NBL1H
	DO 10 I2=1,NBL2H
	   IU1=IU1B1(I1)+IU1B2(I2)
	   ICHART=ICHARB1(I1)+ICHARB2(I2)
	    DO J=1,ICONT
	     IF((IU1==IU1B1SIS(J)).AND.(ICHART==ICHARB1SIS(J))) THEN
	       IDIMB1HSIS(J)=IDIMB1HSIS(J)+IDIMB1H(I1)*IDIMB2H(I2)
! 	       IF(IDIMB1HSIS(J).GT.IDMAX) THEN
! 	        WRITE(334,*) ' IDMAX < IDSIS ',IDIMB1HSIS(J)
! !1	        STOP
! 	       END IF
	       GOTO 10
	     END IF
	   END DO

05	   IF(J>ICONT) THEN
	      ICONT=ICONT+1
	      IF(ICONT>NBMAX) THEN
	        WRITE(334,*) 'NBMAX nao e suficiente'
!1	        STOP
	      END IF
	      IU1B1SIS(ICONT)=IU1
	      ICHARB1SIS(ICONT)=ICHART
	      IDIMB1HSIS(ICONT)=IDIMB1H(I1)*IDIMB2H(I2)
	      NBL1HSIS=ICONT	     
	   END IF
10	CONTINUE

!*--> ordenado setores
!*1 iu1
	DO J=1,NBL1HSIS
	  DO I=J+1,NBL1HSIS
	  IF(IU1B1SIS(J)>IU1B1SIS(I)) THEN
	    IIA=IU1B1SIS(J)
	    IIB=ICHARB1SIS(J)
	    IIC=IDIMB1HSIS(J)
	    
	    IU1B1SIS(J)=IU1B1SIS(I)
	    ICHARB1SIS(J)=ICHARB1SIS(I)
	    IDIMB1HSIS(J)=IDIMB1HSIS(I)
	    
	    IU1B1SIS(I)=IIA
	    ICHARB1SIS(I)=IIB
	    IDIMB1HSIS(I)=IIC
	  END IF
	  END DO
	END DO  
!*2 ichar
	DO J=1,NBL1HSIS
	  DO I=J+1,NBL1HSIS
	  IF(IU1B1SIS(J)==IU1B1SIS(I)) THEN
	  IF(ICHARB1SIS(I)<ICHARB1SIS(J)) THEN
	    IIA=IU1B1SIS(J)
	    IIB=ICHARB1SIS(J)
	    IIC=IDIMB1HSIS(J)

	    IU1B1SIS(J)=IU1B1SIS(I)
	    ICHARB1SIS(J)=ICHARB1SIS(I)
	    IDIMB1HSIS(J)=IDIMB1HSIS(I)

	    IU1B1SIS(I)=IIA
	    ICHARB1SIS(I)=IIB
	    IDIMB1HSIS(I)=IIC
	  END IF
	  END IF
	  END DO
	END DO

! 	dimension maior bloco
	NBBIG=1
	DO J=1,NBL1HSIS
!	PRINT*,J,ICHARB1SIS(J),IU1B1SIS(J)
	IF(IDIMB1HSIS(J)>NBBIG)NBBIG=IDIMB1HSIS(J)
	END DO
!	PRINT*,'NBBIG=',NBL1HSIS,NBBIG
	



	
	
!*--> bloques/sectores/dimensiones de CUP= Cup^{+}
	NBL1CUPSIS=0
	DO J=1,NBL1HSIS
!* CUP implica aumentar la carga en 1 y el spin en 1

	    IU1=IU1B1SIS(J)+1
	    ICHART=ICHARB1SIS(J)+1

	    DO I=1,NBL1HSIS
	      IU1A=IU1B1SIS(I)
	      ICHARA=ICHARB1SIS(I)

	      IF((IU1==IU1A).AND.(ICHART==ICHARA)) THEN
	        NBL1CUPSIS=NBL1CUPSIS+1

	        IU1B1XCUPSIS(NBL1CUPSIS)=IU1A
	        IU1B1YCUPSIS(NBL1CUPSIS)=IU1-1

	        ICHARB1XCUPSIS(NBL1CUPSIS)=ICHARA
	        ICHARB1YCUPSIS(NBL1CUPSIS)=ICHART-1

	        IDIMB1CUPSIS(1,NBL1CUPSIS)=IDIMB1HSIS(I)
	        IDIMB1CUPSIS(2,NBL1CUPSIS)=IDIMB1HSIS(J)
	        
	        IXHCUPSIS(I)=J
	        IYHCUPSIS(J)=I
	      END IF
	      END DO
	END DO
! CONECT  CUPSIS-->HSIS
!    PRINT*,'NBL1CUPSIS=',NBL1CUPSIS
    DO J=1,NBL1CUPSIS
          IU1XX=IU1B1XCUPSIS(J)
          ICHARXX=ICHARB1XCUPSIS(J)
          IU1YY=IU1B1YCUPSIS(J)
          ICHARYY=ICHARB1YCUPSIS(J)
          DO IFIND=1,NBL1HSIS
            IU1T=IU1B1SIS(IFIND)
            ICHART=ICHARB1SIS(IFIND)
            IF(IU1T.EQ.IU1XX.AND.ICHART.EQ.ICHARXX) IB1CUPSISX(J)=IFIND
	    IF(IU1T.EQ.IU1YY.AND.ICHART.EQ.ICHARYY) IB1CUPSISY(J)=IFIND           
         END DO
!        PRINT*,J,IB1CUPSISX(J),IB1CUPSISY(J)
	END DO	
	
	
!-----> Conections  Cup sist------Hsist
	!* ESPACIO 1
!	PRINT*,'CONECTCUPSIS'
        DO IFIND=1,NBL1HSIS
            IU1T=IU1B1SIS(IFIND)
            ICHART=ICHARB1SIS(IFIND)
          DO J=1,NBL1CUPSIS
          IU1XX=IU1B1XCUPSIS(J)
          ICHARXX=ICHARB1XCUPSIS(J)
          IU1YY=IU1B1YCUPSIS(J)
          ICHARYY=ICHARB1YCUPSIS(J)
            IF(IU1T.EQ.IU1XX.AND.ICHART.EQ.ICHARXX) INVHCUPSISX(IFIND)=J
!	    IF(IU1T.EQ.IU1YY.AND.ICHART.EQ.ICHARYY) IB1CUPSISY(J)=IFIND
          END DO
!        PRINT*,J,IB1CUPSISX(J),IB1CUPSISY(J)
	END DO
         
!         PRINT*,'OJO---->'
! 	DO J=1,NBL1HSIS
! 	PRINT*,J,IXHCUPSIS(J),IYHCUPSIS(J),INVHCUPSISX(J)
!         END DO

	
	
	
!*--> bloques/sectores/dimensiones de CDOWN= Cdown^{+}
	NBL1CDOWNSIS=0
	DO J=1,NBL1HSIS
!* CDOWN implica aumentar la carga en 1
	    IU1=IU1B1SIS(J)-1
	    ICHART=ICHARB1SIS(J)+1

	    DO I=1,NBL1HSIS
	      IU1A=IU1B1SIS(I)
	      ICHARA=ICHARB1SIS(I)

	      IF((IU1==IU1A).AND.(ICHART==ICHARA)) THEN
	        NBL1CDOWNSIS=NBL1CDOWNSIS+1
	        IU1B1XCDOWNSIS(NBL1CDOWNSIS)=IU1A
	        IU1B1YCDOWNSIS(NBL1CDOWNSIS)=IU1+1


	        ICHARB1XCDOWNSIS(NBL1CDOWNSIS)=ICHARA
	        ICHARB1YCDOWNSIS(NBL1CDOWNSIS)=ICHART-1

	        IDIMB1CDOWNSIS(1,NBL1CDOWNSIS)=IDIMB1HSIS(I)
	        IDIMB1CDOWNSIS(2,NBL1CDOWNSIS)=IDIMB1HSIS(J)
	        
	        IXHCDOWNSIS(I)=J
	        IYHCDOWNSIS(J)=I
	      END IF
	      END DO
	END DO
	
! CONECT  CDOWNSIS-->HSIS	
!    PRINT*,'NBL1CDOWNSIS=',NBL1CDOWNSIS
        DO J=1,NBL1CDOWNSIS
          IU1XX=IU1B1XCDOWNSIS(J)
          ICHARXX=ICHARB1XCDOWNSIS(J)
          IU1YY=IU1B1YCDOWNSIS(J)
          ICHARYY=ICHARB1YCDOWNSIS(J)
          DO IFIND=1,NBL1HSIS
            IU1T=IU1B1SIS(IFIND)
            ICHART=ICHARB1SIS(IFIND)
            IF(IU1T.EQ.IU1XX.AND.ICHART.EQ.ICHARXX) IB1CDOWNSISX(J)=IFIND
	    IF(IU1T.EQ.IU1YY.AND.ICHART.EQ.ICHARYY) IB1CDOWNSISY(J)=IFIND           
         END DO
!        PRINT*,J,IB1CDOWNSISX(J),IB1CDOWNSISY(J)
	END DO		
	
!-----> Conections  CDOWN sist------Hsist
	!* ESPACIO 1
!	PRINT*,'CONECTCDOWNSIS'
        DO IFIND=1,NBL1HSIS
            IU1T=IU1B1SIS(IFIND)
            ICHART=ICHARB1SIS(IFIND)
          DO J=1,NBL1CDOWNSIS
          IU1XX=IU1B1XCDOWNSIS(J)
          ICHARXX=ICHARB1XCDOWNSIS(J)
          IU1YY=IU1B1YCDOWNSIS(J)
          ICHARYY=ICHARB1YCDOWNSIS(J)
            IF(IU1T.EQ.IU1XX.AND.ICHART.EQ.ICHARXX) INVHCDOWNSISX(IFIND)=J
!	    IF(IU1T.EQ.IU1YY.AND.ICHART.EQ.ICHARYY) IB1CUPSISY(J)=IFIND
          END DO
!        PRINT*,J,IB1CUPSISX(J),IB1CUPSISY(J)
	END DO	
	
!         PRINT*,'OJO---->'
! 	DO J=1,NBL1HSIS
! 	PRINT*,J,IXHCDOWNSIS(J),IYHCDOWNSIS(J),INVHCDOWNSISX(J)
!         END DO
	

!CCCCCCCCC ---------> IGSIS
	DO 100 J=1,NBL1HSIS

!	ALLOCATE( INVSIS1A(J)%op5(IDIMB1HSIS(J))  )
!	ALLOCATE( INVSIS1B(J)%op5(IDIMB1HSIS(J))  )
!	ALLOCATE( INVSIS2A(J)%op5(IDIMB1HSIS(J))  )
!	ALLOCATE( INVSIS2B(J)%op5(IDIMB1HSIS(J))  )

	ICONTSEC=0

	DO 50 IB1=1,NBL1H
	DO 50 IB2=1,NBL2H

	IU1=IU1B1(IB1)+IU1B2(IB2)
	ICHART=ICHARB1(IB1)+ICHARB2(IB2)

	IF((IU1==IU1B1SIS(J)).AND.(ICHART==ICHARB1SIS(J))) THEN	!40

	ALLOCATE( IGSIS(IB1,IB2)%op4(IDIMB1H(IB1),IDIMB2H(IB2))  )

	  DO 20 I1=1,IDIMB1H(IB1)
	  DO 20 I2=1,IDIMB2H(IB2)

	  ICONTSEC=ICONTSEC+1

	  IGSIS(IB1,IB2)%op4(I1,I2)=ICONTSEC

!	  INVSIS1A(J)%op5(ICONTSEC)=I1
!	  INVSIS1B(J)%op5(ICONTSEC)=IB1

!	  INVSIS2A(J)%op5(ICONTSEC)=I2
!	  INVSIS2B(J)%op5(ICONTSEC)=IB2

20	CONTINUE
40	END IF
50	CONTINUE
100	CONTINUE



	
	END SUBROUTINE ESTRUCTURA
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE DEALLOC_ESP1
	USE PRINCIPAL
	USE AMODESP1

	IMPLICIT NONE
	INTEGER :: J


	DO J=1,NBL1H
	  DEALLOCATE(HB1(J)%op)
	  DEALLOCATE(ANN(J)%op)
	END DO

	DO J=1,NBL1CUP
	  DEALLOCATE(ACUP(J)%op)
	  DEALLOCATE(BCUP(J)%op)
	END DO

	DO J=1,NBL1CDOWN
	  DEALLOCATE(ACDOWN(J)%op)
	  DEALLOCATE(BCDOWN(J)%op)
	END DO


	END SUBROUTINE DEALLOC_ESP1
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE MULTI
	USE PRINCIPAL
	USE AMODOLO
	USE AMODESP1
	USE AMODESP2
	USE AMODHSIS
	USE AMODIGSIS
	
	IMPLICIT NONE
	INTEGER ::J,II,IJ,IT
	
	IF(LBL.EQ.1)THEN
!* ESPACIO 1
	DO J=1,NBL1H
	  ALLOCATE(HB1(J)%op(IDIMB1H(J),IDIMB1H(J)))
	  HB1(J)%op=0.D0

	 DO II=1,IDIMB1H(J)
	 DO IJ=1,IDIMB1H(J)
	HB1(J)%op(II,IJ)=ACOPU*XLU(J)%op(II,IJ) + (ACOPn+ACS(LBL))*XLN(J)%op(II,IJ) - ACOPH*XLH(J)%op(II,IJ)

	 END DO
	 END DO
	END DO
	END IF
	
	
!* ESPACIO 2
	DO J=1,NBL2H
	  ALLOCATE(X2(J)%op(IDIMB2H(J),IDIMB2H(J)))
	  DO II=1,IDIMB2H(J)
	 DO IJ=1,IDIMB2H(J)
	X2(J)%op(II,IJ)=ACOPU*XLU(J)%op(II,IJ) + (ACOPn+ACS(LBL+1))*XLN(J)%op(II,IJ) - ACOPH*XLH(J)%op(II,IJ)

	 END DO
	 END DO
	END DO
	
	
		
	END SUBROUTINE MULTI
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	SUBROUTINE VALIN
	USE PRINCIPAL
	USE AMODOLO
	USE AMODESP1
	USE AMODESP2
	USE AMODHSIS
	USE AMODIGSIS
	
	IMPLICIT NONE

	INTEGER :: J, I,II,IJ,NS1
	INTEGER, DIMENSION(:,:), ALLOCATABLE ::IDU,IS


!*-------------------------- Aqui se define H -------------------
	DO J=1,2
	DO II=1,1
	DO IJ=1,1
	IIDUPLA2(J,II,IJ)=0
	IISIMPLE2(J,II,IJ)=0
	END DO
	END DO
	END DO
	
	

	NBL1H=4
	NBL2H=4

	DO J=1,NBL1H
	 IDIMB1H(J)=1
	END DO

	DO J=1,NBL1H
	 IDIMB2H(J)=IDIMB1H(J)
	END DO

!* ESPACIO 1
	ICHARB1(1)=1
	IU1B1(1)=-1
	ICHARB1(2)=0
	IU1B1(2)=0
	ICHARB1(3)=2
	IU1B1(3)=0
	ICHARB1(4)=1
	IU1B1(4)=1

!* ESPACIO 2-3-4
	DO J=1,NBL1H
	  ICHARB2(J)=ICHARB1(J)
	  IU1B2(J)=IU1B1(J)
	END DO	

!*-------------- Aqui se define los Operadores locales -------
!* ESPACIO
	DO J=1,NBL1H
	  ALLOCATE(XLU(J)%op(IDIMB1H(J),IDIMB1H(J)))
	  XLU(J)%op=0.d0
	END DO

	XLU(3)%op(1,1)=1.0D0

	DO J=1,NBL1H
	  ALLOCATE(XLN(J)%op(IDIMB1H(J),IDIMB1H(J)))
	  XLN(J)%op=0.d0
	  ALLOCATE(ANN(J)%op(IDIMB1H(J),IDIMB1H(J)))
	  ANN(J)%op=0.d0
	END DO

	XLN(1)%op(1,1)=1.0D0
	XLN(3)%op(1,1)=2.0D0
	XLN(4)%op(1,1)=1.0D0
	
	ANN(1)%op(1,1)=1.0D0
	ANN(3)%op(1,1)=2.0D0
	ANN(4)%op(1,1)=1.0D0
	
        DO J=1,NBL1H
          ALLOCATE(XLH(J)%op(IDIMB1H(J),IDIMB1H(J)))
          XLH(J)%op=0.d0
        END DO

        XLH(1)%op(1,1)=-1.0D0
        XLH(4)%op(1,1)=1.0D0
	
! * ESPACIO 1
! 	DO J=1,NBL1H
! 	  ALLOCATE(HB1(J)%op(IDIMB1H(J),IDIMB1H(J)))
! 	  HB1(J)%op=0.D0
! 
! 	 DO II=1,IDIMB1H(J)
! 	 DO IJ=1,IDIMB1H(J)
! 	HB1(J)%op(II,IJ)=ACOPU*XLU(J)%op(II,IJ) + ACOPn*XLN(J)%op(II,IJ)
! 
! 	 END DO
! 	 END DO
! 	END DO
! * ESPACIO 2
! 	DO J=1,NBL2H
! 	  ALLOCATE(X2(J)%op(IDIMB2H(J),IDIMB2H(J)))
! 	  X2(J)%op=HB1(J)%op
! 	END DO
!*-------------------------- Aqui se define CUP  ---------------------

	NBL1CUP=2
	NBL2CUP=2
!* IDIMB1CUP(i,k) -> o bloco k do oper. CUP: possui dimensao x = idimb1CUP(1,k)
!*                                               dimensao y = idimb1CUP(2,k)
!* i=1-> x, i=2-> y
!* Espacio 1 
	IDIMB1CUP(1,1)=1
	IDIMB1CUP(2,1)=1
	IDIMB1CUP(1,2)=1
	IDIMB1CUP(2,2)=1
!* Espacio 2
	DO J=1,NBL1CUP
	DO I=1,2
	IDIMB2CUP(I,J)=IDIMB1CUP(I,J)
	END DO
	END DO

!* ISIMB1XCUP(k) -> direccion x, k:bloco 
!* Espacio 1
	ICHARB1XCUP(1)=2
	IU1B1XCUP(1)=0
	ICHARB1YCUP(1)=1
	IU1B1YCUP(1)=-1
	ICHARB1XCUP(2)=1
	IU1B1XCUP(2)=1
	ICHARB1YCUP(2)=0
	IU1B1YCUP(2)=0


	DO J=1,NBL1CUP
	ICHARB2XCUP(J)=ICHARB1XCUP(J)
	IU1B2XCUP(J)=IU1B1XCUP(J)
	ICHARB2YCUP(J)=ICHARB1YCUP(J)
	IU1B2YCUP(J)=IU1B1YCUP(J)
	END DO

!* ESPACIO
	DO J=1,NBL1CUP
	  ALLOCATE(CUPL(J)%op(IDIMB1CUP(1,J),IDIMB1CUP(2,J)))
	  CUPL(J)%op=0.D0
	END DO

	CUPL(1)%op(1,1)=1.D0

	CUPL(2)%op(1,1)=1.D0

	DO J=1,NBL1CUP
	  ALLOCATE(ACUP(J)%op(IDIMB1CUP(1,J),IDIMB1CUP(2,J)))
	  ACUP(J)%op=CUPL(J)%op
	  ALLOCATE(BCUP(J)%op(IDIMB1CUP(1,J),IDIMB1CUP(2,J)))
	  BCUP(J)%op=CUPL(J)%op
	  
	  ALLOCATE(CUPS(J)%op(IDIMB1CUP(1,J),IDIMB1CUP(2,J)))
	  CUPS(J)%op=0.D0
	  ALLOCATE(CUPD(J)%op(IDIMB1CUP(1,J),IDIMB1CUP(2,J)))
	  CUPD(J)%op=0.D0
	END DO	

	
	
!	IIDUPLA(1,1,1,1) -> sitio 1,bloco 1, posicao 1,1 -> dupla ocpuacao
!	IIDUPLA(1,1,1,1)=1
!	IISIMPLE(1,2,1,1) -> sitio 1, bloco 2, posicao 1,1 -> ocpuacao simples
!	IISIMPLE(1,2,1,1)=1
	
	IIDUPLA2(1,1,1)=1
	IISIMPLE2(2,1,1)=1

	CUPS(2)%op(1,1)=1.D0

	CUPD(1)%op(1,1)=1.D0
	
	DO J=1,NBL1CUP
	NS1=10000+1000*1+J
	WRITE(NS1)IDIMB1CUP(1,J),IDIMB1CUP(2,J)
	WRITE(NS1)CUPS(J)%op
	REWIND(NS1)
	CLOSE(NS1)
	DEALLOCATE(CUPS(J)%op)
	END DO
	
	DO J=1,NBL1CUP
	NS1=20000+1000*1+J
	WRITE(NS1)IDIMB1CUP(1,J),IDIMB1CUP(2,J)
	WRITE(NS1)CUPD(J)%op
	REWIND(NS1)
	CLOSE(NS1)
	DEALLOCATE(CUPD(J)%op)
	END DO
	
	
!*-------------------------- Aqui se define CDOWN  ---------------------

	NBL1CDOWN=2
	NBL2CDOWN=2
!* IDIMB1CDOWN(i,k) -> o bloco k do oper. CDOWN: possui dimensao x = idimb1CDOWN(1,k)
!*                                               dimensao y = idimb1CDOWN(2,k)
!* i=1-> x, i=2-> y
!* Espacio 1
	IDIMB1CDOWN(1,1)=1
	IDIMB1CDOWN(2,1)=1
	IDIMB1CDOWN(1,2)=1
	IDIMB1CDOWN(2,2)=1

!* Espacio 2
	DO J=1,NBL1CDOWN
	DO I=1,2
	IDIMB2CDOWN(I,J)=IDIMB1CDOWN(I,J)
	END DO
	END DO
!* ISIMB1XCDOWN(k) -> direccion x, k:bloco
!* Espacio 1
	ICHARB1XCDOWN(1)=1
	IU1B1XCDOWN(1)=-1
	ICHARB1YCDOWN(1)=0
	IU1B1YCDOWN(1)=0
	ICHARB1XCDOWN(2)=2
	IU1B1XCDOWN(2)=0
	ICHARB1YCDOWN(2)=1
	IU1B1YCDOWN(2)=1

	DO J=1,NBL1CDOWN
	ICHARB2XCDOWN(J)=ICHARB1XCDOWN(J)
	IU1B2XCDOWN(J)=IU1B1XCDOWN(J)
	ICHARB2YCDOWN(J)=ICHARB1YCDOWN(J)
	IU1B2YCDOWN(J)=IU1B1YCDOWN(J)
	END DO

!* ESPACIO
	DO J=1,NBL1CDOWN
	  ALLOCATE(CDOWNL(J)%op(IDIMB1CDOWN(1,J),IDIMB1CDOWN(2,J)))
	  CDOWNL(J)%op=0.d0
	END DO

	CDOWNL(1)%op(1,1)=1.D0

	CDOWNL(2)%op(1,1)=-1.D0

	DO J=1,NBL1CDOWN
	  ALLOCATE(ACDOWN(J)%op(IDIMB1CDOWN(1,J),IDIMB1CDOWN(2,J)))
	  ACDOWN(J)%op=CDOWNL(J)%op
	  ALLOCATE(BCDOWN(J)%op(IDIMB1CDOWN(1,J),IDIMB1CDOWN(2,J)))
	  BCDOWN(J)%op=CDOWNL(J)%op
	  
	  ALLOCATE(CDOWNS(J)%op(IDIMB1CDOWN(1,J),IDIMB1CDOWN(2,J)))
	  CDOWNS(J)%op=0.D0
	  ALLOCATE(CDOWND(J)%op(IDIMB1CDOWN(1,J),IDIMB1CDOWN(2,J)))
	  CDOWND(J)%op=0.D0
	END DO	
	
	FFDUPLA2(2,1,1)=1
	FFSIMPLE2(1,1,1)=1
	
	CDOWNS(1)%op(1,1)=1.D0
	
	CDOWND(2)%op(1,1)=-1.D0
	
	DO J=1,NBL1CDOWN
	NS1=30000+1000*1+J
	WRITE(NS1)IDIMB1CDOWN(1,J),IDIMB1CDOWN(2,J)
	WRITE(NS1)CDOWNS(J)%op
	REWIND(NS1)
	CLOSE(NS1)
	DEALLOCATE(CDOWNS(J)%op)
	END DO
	
	DO J=1,NBL1CDOWN
	NS1=40000+1000*1+J
	WRITE(NS1)IDIMB1CDOWN(1,J),IDIMB1CDOWN(2,J)
	WRITE(NS1)CDOWND(J)%op
	REWIND(NS1)
	CLOSE(NS1)
	DEALLOCATE(CDOWND(J)%op)
	END DO
	
	END SUBROUTINE VALIN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


       Function Fmin(ax,bx,f,Tol)
       Implicit double precision (a-h,o-z)
       External f
        C=0.5d0*(3.0d0-dsqrt(5.0d0))
       Eps=1.0d0
10     Eps=Eps/2.0d0
       Tol1=1.0d0+Eps
       If(tol1.gt.1.0d0) go to 10
       Eps=dsqrt(Eps)
       A=ax
       B=bx
       V=A+C*(B-A)
       W=V
       X=V
       E=0.0d0
       FX=f(X)
       FV=FX
       FW=FX
20     XM=0.5d0*(A+B)
       TOL1=Eps*Dabs(X)+tol/3.0d0
       TOL2=2.0d0*TOL1


       If(Dabs(X-XM).le.(TOL2-0.5d0*(B-A))) go to 90

       If(Abs(E).le.Tol1) go to 40

       R=(X-W)*(FX-FV)
       Q=(X-V)*(FX-FW)
       P=(X-V)*Q-(X-W)*R
       Q=2.0d0*(Q-R)
       If(Q.gt.0.0d0) P=-P
       Q=Dabs(Q)
       R=E
       E=G
!
!   � Parabola es aceptable ?
!
30     If(dabs(P).ge.dabs((0.5d0*Q*R))) go to 40
       If(P.ge.Q*(A-X)) go to 40
       If(P.ge.Q*(B-X)) go to 40
!
!   El paso de la interpolacion parabolica
!
       G=P/Q
       U=X+G
!
!   F no se permite calcular demasiado cerca al AX o BX
!
       If((U-A).lt.Tol2) G=Dsign(Tol1,XM-X)
       If((B-U).lt.Tol2) G=Dsign(Tol1,XM-X)
       Go to 50
!
!   El paso de seccion de oro
!
40     If(X.ge.XM) E=A-X
       If(X.lt.XM) E=B-X
       G=C*E
!
!   F no se permite calcular demasiado cerca al X
!
50     If(Dabs(G).ge.Tol1) U=X+G
       If(Dabs(G).lt.Tol1) U=X+Dsign(Tol1,G)
       FU=f(u)
!
!  Asignar valores nuevos a los parametros A,B,V,W y X
!
       If(FU.gt.FX) go to 60
       If(U.ge.X) A=X
       If(U.lt.X) B=X
       V=W
       FV=FW
       W=X
       FW=FX
       X=U
       FX=FU
       go to 20
60     if(U.lt.X) A=U
       If(U.ge.X) B=U
       If(FU.le.FW) go to 70
       If(W.eq.X) go to 70
       If(FU.le.FV) go to 80
       If(V.eq.X) go to 80
       If(V.eq.W) go to 80
       go to 20
70     V=W
       FV=FW
       W=U
       FW=FU
       go to 20
80     V=U
       FV=FU
       go to 20
!
!   Fin del ciclo principal
!
90     Fmin=X
       Return
       End
