! Comments:
! 060907: Introduced the evaluation of derivatives of F using five-point
! 	derivative (used in the coefficients of the equations for wave
!	amplitudes).
! 060912: Introduced the evaluation of wave and particle energies 
!	(subroutine Energy).
! 060921: Introduces the possibility of choice of the resonance condition
! 	in the evaluation of the Ai and Dij coefficients, through the
!	parameter ADResCondition, which can be "Exact" or "Approx".
! 060923: Introduces the equation for ion-sound waves.
! 060925: Introduces the contribution of induced and spontaneous decay to
!       the equation for the evolution of L waves.
! 061113: Introduces the possibility of choice of the resonance condition
! 	in the evaluation of the Ldecay coefficients, through the
!	parameter LdecayResCondition, which can be "Exact" or "Approx".
! 061226: Introduces the contribution of induced and spontaneous decay to
!       the equation for the evolution of S waves.
! 070104: Introduces the symmetry in 'qx'. Uses the formulation appearing
!	in notes J-06, dated January 04, 2007.
! 070119: Introduces the symmetry in 'ux', and reformulated some of the
!	formulation. 
!	Introduces the possibility of choice of the resonance condition
! 	in the evaluation of the Sdecay coefficients, through the
!	parameter SdecayResCondition, which can be "Exact" or "Approx".
!	Uses the formulation appearing in notes J-06, dated 
!	January 19, 2007.
! 070123: Introduces the use of five-point derivative for the explicit
!	derivative in the ADI routine. Uses the formulation of notes J-06, 
!	dated January, 23, 2007. Introduces the possibility of choice between
! 	the three-point derivative (subroutine 'Adi3p') and five-point
!	derivative (subroutine 'Adi5p').
! 070203: Introduces the use of a Runge-Kutta procedure with adjustable step
!	(based on Numerical Recipes) for the evaluation of the wave equations. 
!	After each step done, uses the value of DTau obtained with
!	the Runge-Kutta for a step in the particle equation, using 'Adi3p'
!	or 'Adi5p'.
! 070215: Based on version 070203 (Runge-Kutta for wave equations). Uses the
!	formulation appearing in notes A-07, dated February 15, 2007 (the
!	integrations are performed over 'uz' and 'qz', looking for similarity
!	with the 1-D case). Separates the routine Coef_AD into Coef_A and
!	Coef_D. This saves time since the Ai are not modified along time
!	evolution. Evaluates the normalization of the 'Fe", which is saved
!	as a 6th. column in the file 'E_ratio.wt'. It is possible to 
!	re-normalize at each time step, but this possibility has been 
!	commented out (see 'Anorm commented out').
! 070313: Introduces the contribution of induced and spontaneous scattering
!	to the equation for the evolution of L waves (notes A-07, dated
!	March 21, 2007).
! 070323: Modifies the boundary conditions, introducing constant derivatives
!       of the distribution function at ux_lim and \pm uz_lim, instead of
!       constant value of the distribution. The new boundary condition allows
!	for change of the distribution at the boundaries (see notes A-07,
!	dated March 21, 2007).
!	A mistake in the saving of the boundary condition, for continued
!	calculation, was corrected in September 27, 2007.
! 070425: Introduces the possibility of occurrence of two beams in the
!	initial distribution function, one forward and another backward.
!	Corrects a mistake in the term for scattering of L waves,
!	introduced in version 070313.
! 070908: Introduces the evaluation of 1D projections of the electron
!       distribution and of the L and S spectra.
!       Corrects a mistake of previous versions: Fi is now saved properly, for
!       continued calculation. The mistake affected coefficients for S wave
!       spontaneous and induced emission.
! 071001: Introduces the choice between Runge-Kutta with fixed time step or
!	with adaptative time step, through parameter 'RKStep', which may be
!	'Fixed' or 'Adapt'.
! 071017: The integrals over Fi in the quasilinear term of the equation for
!	S waves is now performed analytically, instead of numerically (see
!	the difference between Eq. (54) and Eq. (55), in notes A07, version
!	17 October 2007. 
!	Another mistake was corrected: A factor "Muq" was removed from the
!	decay coefficient for S waves (this does not constitute a new
!	version). Another small improvement: A variable was created to
!	replace 'Zsq/Q' in the initialization of the 'S' spectra and in other
! 	places related to 'S' waves. This allows for analytical evaluation of
!	a quotient 'Q/Q', avoiding trouble for small Q.
! 071217: Some conditional branching were introduced, to avoid trouble at small
!	values of "Qx" and "Qz". These trouble appeared mainly in routine
!	'Coef_Swave' and in the scattering term in 'Coef_Lwave', and caused
!	'running errors' with the Fortran of the Alpha workstations,
!	explosive instabilities when using the option of compilation '-axN'
!	of the Intel Fortran, and explosive behavior when including the
!	scattering term, using the version of the compiler installed in the
!	cluster 'minuano'. After the corrections, some tests have shown the
!	same (well-behaved) results using different versions of compiler, and 
!	using or not the '-axN' option.
! 071221: Modified the procedure of evaluation of the some terms appearing
!	in the decay coefficients, in subroutines 'Aux_Coef_Lwave' and 
!	'Aux_Coef_Swave'. Introduced the subroutine 'Aitp2d' in order to 
!	interpolate over the wave spectra, in 2D. As a consequence of this
!	modification, variable 'ídif' has been eliminated from these routines
!	and from 'Coef_Lwave' and 'Coef_Swave'.

!---------------------------------------------------------------------
! List of modules, subroutines, and functions:

! MODULE Common_Arrays
! MODULE Common_Params
! MODULE Math_Constants
! MODULE Phys_Constants

! SUBROUTINE Definitions(RMiMe,AA,Geff)
! SUBROUTINE Init_Wave(Ulim,Qxi,Qxf,Qzi,Qzf,InitialLevel,Iw0,RatioNf,Uf,&
!     RTfTe,RatioNb,Ub,RTbTe,Tau,Dqx,Dqz,Dux,Duz,Anorm,Epart,Ewave,EppEw)
! SUBROUTINE Save_Results(Ulim,Ucrit,Qxi,Qxf,Qzi,Qzf,InitialLevel,Iw0,&
!     RatioNf,Uf,RTfTe,RatioNb,Ub,RTbTe,RTeTi,G,ADResCondition,&
!     LdecayResCondition,SdecayResCondition,Lemis,Ldecay,Lscat,&
!     Semis,Sdecay,Sscat,AdiRoutine,BoundaryCondition,RKStep,Tau,&
!     Dqx,Dqz,Dux,Duz,Anorm0,Epart0,Ewave0,EppEw0,Rn,Rp,Rw,Rs)
! SUBROUTINE Read_Results(Tau,Dqx,Dqz,Dux,Duz,Anorm0,Epart0,Ewave0,EppEw0,&
!     Rn,Rp,Rw,Rs)
! SUBROUTINE Coef_A
! SUBROUTINE Coef_D
! SUBROUTINE Coef_Lwave(Qx,Qz,iqx,kqz,sigma,Dfdux,Dfduz,Iwave,Iwavem,Coef)
! SUBROUTINE Aux_Coef_Lwave(sigma,SignCase,iqxp,Qx,Qz,Qxp,Qzp,Qxdif,Qzdif,&
!     Aux0,Aux1,Zlqp,Zlqdif,Iqp,Iqdif)
! REAL FUNCTION Exp_Sqr_Lwave(Aux0,A,B)
! SUBROUTINE Coef_Swave(Qx,Qz,iqx,kqz,sigma,Dfdux,Dfduz,Iwave,Coef)
! SUBROUTINE Aux_Coef_Swave(SignCase,iqxp,Qxp,Qzp,Qxdif,Qzdif,&
!     Aux1,Zlqp,Zlqdif,Iqp,Iqdif)
! SUBROUTINE Fnorm(Anorm)
! SUBROUTINE Energy(Epart,Ewave,EppEw)
! SUBROUTINE Output(WriteChoice)
! SUBROUTINE Adi3p(Ucrit,Dux,Duz,DTau,BoundaryCondition)
! SUBROUTINE Adi5p(Ucrit,Dux,Duz,DTau,BoundaryCondition)
! SUBROUTINE Tridag(If,k,A,B,C,D,V,Beta,Gamma)
! SUBROUTINE Cor_Ampli(Vf,n)
! SUBROUTINE Aitp1d(nx,Vx,Fx,Xp,Fp)
! SUBROUTINE Aitp2d(nx,ny,Vx,Vy,Fxy,Xp,Yp,Fp)
! SUBROUTINE Locate(Xx,N,X,J)
! SUBROUTINE Simpson(Vx,F,N,Res)
! SUBROUTINE Derivxy5p(nx,ny,Vx,Vy,Fxy,Dfdx,Dfdy)
! SUBROUTINE RKQC(Y,DYDX,N,X,HTRY,EPS,YSCAL,HDID,HNEXT)
! SUBROUTINE RK4(Y,DYDX,N,X,H,YOUT)

! REAL FUNCTION ZL(Qx,Qz)
! REAL FUNCTION ZS(Qx,Qz)

!---------------------------------------------------------------------

!---------------------------------------------------------------------
! Modules with definitions:
!---------------------------------------------------------------------

MODULE Common_Arrays
IMPLICIT NONE

INTEGER, PARAMETER :: nqx=51,nqz=51,nux=51,nuz=101
! nqx, nqz, nux, and nuz must be odd numbers.
REAL, DIMENSION(nqx) :: VQx
REAL, DIMENSION(nqz) :: VQz
REAL, DIMENSION(nqx,nqz) :: ILp,ILm,ISp,ISm
REAL, DIMENSION(nux) :: VUx
REAL, DIMENSION(nuz) :: VUz
REAL, DIMENSION(nux,nuz) :: Fe0,Fe,Fe_old,Fi
REAL, DIMENSION(nux,nuz) :: Ax,Az,Dxx,Dxz,Dzx,Dzz
REAL, DIMENSION(nux) :: Fepzm,Fepzp
REAL, DIMENSION(nuz) :: Fepxp

END MODULE Common_Arrays

MODULE Common_Params
IMPLICIT NONE

REAL :: G,RTeTi
CHARACTER(LEN=6) :: ADResCondition
CHARACTER(LEN=6) :: LdecayResCondition
CHARACTER(LEN=3) :: Lemis,Ldecay,Lscat
CHARACTER(LEN=6) :: SdecayResCondition
CHARACTER(LEN=3) :: Semis,Sdecay,Sscat

END MODULE Common_Params

MODULE Math_Constants
IMPLICIT none

REAL :: Pi= 3.1415926535897932384626433832795
REAL :: Sqtwo= 1.414213562
REAL :: Infinity= 1.0E+30
REAL :: Degree= 0.017453293
REAL :: EpsMin= 1.0E-16
REAL :: Qmin= 1.0E-4

COMPLEX :: Zi= (0.,1.)
COMPLEX :: Zzero= (0.,0.)

END MODULE Math_Constants

MODULE Phys_Constants
IMPLICIT none

REAL :: MeC2= 511.01E+3          ! electron mass (eV)
REAL :: MpC2= 938.50E+6          ! proton mass (eV)
REAL :: C_SI= 2.997925E+8        ! speed of light (m/s)
REAL :: C_cgs= 2.997925E+10      ! speed of light (cm/s)

END MODULE Phys_Constants

SUBROUTINE Definitions(RMiMe,AA,Geff)
USE Common_Arrays
USE Common_Params
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL :: RMpMe,RMeMp
REAL :: RMiMe,AA,Geff

RMpMe= MpC2/MeC2
RMeMp= MeC2/MpC2
RMiMe= RMpMe
AA= SQRT(1.+3./RTeTi)/SQRT(RMiMe)/SQRT(2.)
Geff= G/(2.*SQRT(2.)*(4.*Pi)**2)

END SUBROUTINE Definitions


!---------------------------------------------------------------------
! Main program:
!---------------------------------------------------------------------

PROGRAM WT_2D
USE Common_Arrays
USE Common_Params
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL :: Ulim,Ucrit,Qxi,Qxf,Qzi,Qzf,Iw0
REAL :: RatioNf,Uf,RTfTe
REAL :: RatioNb,Ub,RTbTe
REAL :: H,H1,Hmin,Hnext,Hdid
REAL :: Epsrk
REAL :: DTau,Tau,Tau1,Tau2,TauAdd
REAL :: Dqx,Dqz,Dux,Duz 
REAL :: Ewave0,Ewave,Epart0,Epart,EppEw0,EppEw
REAL :: Rn,Rp,Rw,Rs
REAL :: Anorm0,Anorm
REAL, DIMENSION(4*nqx*nqz) :: Ft,Dfdt,Yscal
INTEGER :: Iflag
INTEGER :: ind,nvar,nok,nbad,i,k
INTEGER :: it,Nitera 
CHARACTER(LEN=6) :: InitialLevel
CHARACTER(LEN=5) :: AdiRoutine
CHARACTER(LEN=5) :: RKStep 
CHARACTER(LEN=10) :: BoundaryCondition

OPEN(1,FILE='Start.wt')
READ(1,*) DTau
READ(1,*) Nitera       ! Useful only in case of RKStep="Adapt"
READ(1,*) TauAdd
CLOSE(1)

OPEN(2,FILE='E_ratio.wt')
OPEN(1,FILE='Ini.wt')
READ(1,*) Iflag
READ(1,*) Ulim
READ(1,*) Ucrit
READ(1,*) Qxi 
READ(1,*) Qxf 
READ(1,*) Qzi 
READ(1,*) Qzf 
READ(1,*) InitialLevel 
READ(1,*) Iw0 
READ(1,*) RatioNf 
READ(1,*) Uf 
READ(1,*) RTfTe
READ(1,*) RatioNb 
READ(1,*) Ub 
READ(1,*) RTbTe
READ(1,*) RTeTi 
READ(1,*) G
READ(1,*) ADResCondition
READ(1,*) Lemis
READ(1,*) Ldecay
READ(1,*) Lscat
READ(1,*) Semis
READ(1,*) Sdecay
READ(1,*) Sscat
READ(1,*) LdecayResCondition
READ(1,*) SdecayResCondition
READ(1,*) AdiRoutine
READ(1,*) BoundaryCondition
READ(1,*) RKStep

IF(Qxi<Qmin .OR. Qzi<Qmin) THEN
 OPEN(98,FILE='Warning_Main_Qi.wt')
 WRITE(98,*) ' Qxi= ',Qxi,'  Qzi= ',Qzi
 WRITE(98,*) ' Qxi and Qzi can not be too small !!'
 CLOSE(98)
 IF(Qxi<Qmin) THEN
  Qxi= Qmin
 ELSE
 END IF
 IF(Qzi<Qmin) THEN
  Qzi= Qmin
 ELSE
 END IF
ELSE
END IF

SELECT CASE(IFlag)

 CASE(0)
  CALL Init_Wave(Ulim,Qxi,Qxf,Qzi,Qzf,InitialLevel,Iw0,RatioNf,Uf,&
     RTfTe,RatioNb,Ub,RTbTe,Tau,Dqx,Dqz,Dux,Duz,Anorm0,Epart0,Ewave0,EppEw0)
  CALL Output("Fe0")
!  CALL Output("Fi ")
  CALL Output("IL0")
  CALL Output("IS0")
  it= 0
  Rn= 1.     ! Anorm/Anorm0
  Rp= 1.     ! Epart/Epart0
  Rw= 1.     ! Ewave/Ewave0
  Rs= 1.     ! EppEw/EppEw0
  WRITE(2,2001) it,Tau,Rn,Rp,Rw,Rs
2001 FORMAT(1x,i5,5(1x,e13.5e3))

 CASE(1)
  CALL Read_Results(Tau,Dqx,Dqz,Dux,Duz,Anorm0,Epart0,Ewave0,EppEw0,&
       Rn,Rp,Rw,Rs)
  it= 0
  WRITE(2,2001) it,Tau,Rn,Rp,Rw,Rs

 CASE DEFAULT
  OPEN(98,FILE='Warning_Main.wt')
  WRITE(98,*) ' Iflag= ',Iflag
  WRITE(98,*) ' Iflag must be 0 or 1 !!'
  CLOSE(98)
  STOP

END SELECT

! Start time evolution:

CALL Coef_A

! Procedure based on subroutine ODEINT (Numerical Recipes):
OPEN(99,FILE='Runge-Kutta.wt')
nvar= 4*nqx*nqz
ind= 0
DO i= 1,nqx
 DO k= 1,nqz
  ind= ind+1
  Ft(ind)= ILp(i,k)
 END DO
END DO
DO i= 1,nqx
 DO k= 1,nqz
  ind= ind+1
  Ft(ind)= ILm(i,k)
 END DO
END DO
DO i= 1,nqx
 DO k= 1,nqz
  ind= ind+1
  Ft(ind)= ISp(i,k)
 END DO
END DO
DO i= 1,nqx
 DO k= 1,nqz
  ind= ind+1
  Ft(ind)= ISm(i,k)
 END DO
END DO

SELECT CASE(RKStep)

 CASE("Fixed")
  H= DTau

  Tau1= Tau
  Tau2= Tau1+TauAdd
  Nitera= 1+(Tau2-Tau1)/DTau

  !***     Time evolution:    ***
  DO it= 1,Nitera
   CALL Derivs(Tau,Ft,Dfdt)
   IF( (Tau+H-Tau2)*(Tau+H-Tau1) > 0. ) H=Tau2-Tau
   CALL RK4(Ft,Dfdt,nvar,Tau,H,Ft)
   Tau= Tau+H
   CALL Cor_Ampli(Ft,nvar)

   ind= 0
   DO i= 1,nqx
    DO k= 1,nqz
     ind= ind+1
     ILp(i,k)= Ft(ind)
    END DO
   END DO
   DO i= 1,nqx
    DO k= 1,nqz
     ind= ind+1
     ILm(i,k)= Ft(ind)
    END DO
   END DO
   DO i= 1,nqx
    DO k= 1,nqz
     ind= ind+1
     ISp(i,k)= Ft(ind)
    END DO
   END DO
   DO i= 1,nqx
    DO k= 1,nqz
     ind= ind+1
     ISm(i,k)= Ft(ind)
    END DO
   END DO

   SELECT CASE(AdiRoutine)
    CASE("Adi3p")
     CALL Adi3p(Ucrit,Dux,Duz,DTau,BoundaryCondition)
    CASE("Adi5p")
     CALL Adi5p(Ucrit,Dux,Duz,DTau,BoundaryCondition)
    CASE DEFAULT
    OPEN(98,FILE='Warning_Main_Program.wt')
    WRITE(98,*) ' AdiRoutine= ',AdiRoutine
    WRITE(98,*) ' AdiRoutine must be (Adi3p) or (Adi5p) !!'
    CLOSE(98)
    STOP
   END SELECT

   CALL Fnorm(Anorm)
   CALL Energy(Epart,Ewave,EppEw)
   Rn= Anorm/Anorm0
   Rp= Epart/Epart0
   Rw= Ewave/Ewave0
   Rs= EppEw/EppEw0
   WRITE(2,2001) it,Tau,Rn,Rp,Rw,Rs
   OPEN(98,FILE='status.wt')
   WRITE(98,*) ' it= ',it,'   Tau= ',Tau
   CLOSE(98)

   IF ( (Tau-Tau2)*(Tau2-Tau1) >= 0. ) GO TO 1
  END DO

 CASE("Adapt")
  H1= DTau
  Hmin=0.1
  Epsrk= 1.E-3

  Tau1= Tau
  Tau2= Tau1+TauAdd
  H= SIGN(H1,Tau2-Tau1)
  nok=0
  nbad=0

  !***     Time evolution:    ***
  DO it= 1,Nitera
   CALL Derivs(Tau,Ft,Dfdt)
   DO i= 1,nvar
    Yscal(i)= ABS(Ft(i)) + ABS(H*Dfdt(i)) + EpsMin
   END DO
   IF( (Tau+H-Tau2)*(Tau+H-Tau1) > 0. ) H=Tau2-Tau
   CALL RKQC(Ft,Dfdt,nvar,Tau,H,Epsrk,Yscal,Hdid,Hnext)
   IF( Hdid == H ) THEN
    nok=nok+1
   ELSE
    nbad=nbad+1
   END IF
   CALL Cor_Ampli(Ft,nvar)

   ind= 0
   DO i= 1,nqx
    DO k= 1,nqz
     ind= ind+1
     ILp(i,k)= Ft(ind)
    END DO
   END DO
   DO i= 1,nqx
    DO k= 1,nqz
     ind= ind+1
     ILm(i,k)= Ft(ind)
    END DO
   END DO
   DO i= 1,nqx
    DO k= 1,nqz
     ind= ind+1
     ISp(i,k)= Ft(ind)
    END DO
   END DO
   DO i= 1,nqx
    DO k= 1,nqz
     ind= ind+1
     ISm(i,k)= Ft(ind)
    END DO
   END DO
   DTau= Hdid
   
   SELECT CASE(AdiRoutine)
    CASE("Adi3p")
     CALL Adi3p(Ucrit,Dux,Duz,DTau,BoundaryCondition)
    CASE("Adi5p")
     CALL Adi5p(Ucrit,Dux,Duz,DTau,BoundaryCondition)
    CASE DEFAULT
     OPEN(98,FILE='Warning_Main_Program.wt')
     WRITE(98,*) ' AdiRoutine= ',AdiRoutine
     WRITE(98,*) ' AdiRoutine must be (Adi3p) or (Adi5p) !!'
     CLOSE(98)
     STOP
   END SELECT

   CALL Fnorm(Anorm)
   CALL Energy(Epart,Ewave,EppEw)
   Rn= Anorm/Anorm0
   Rp= Epart/Epart0
   Rw= Ewave/Ewave0
   Rs= EppEw/EppEw0
   WRITE(2,2001) it,Tau,Rn,Rp,Rw,Rs
   OPEN(98,FILE='status.wt')
   WRITE(98,*) ' it= ',it,'   Tau= ',Tau
   CLOSE(98)

   IF ( (Tau-Tau2)*(Tau2-Tau1) >= 0. ) THEN
    Iflag= 1
    H1= Hnext
    GO TO 1
   ELSE 
   END IF
   IF ( ABS(Hnext) < Hmin ) WRITE(99,*) &
     ' Stepsize smaller than minimum.'
   H= Hnext
  END DO

 CASE DEFAULT
  OPEN(98,FILE='Warning_Main.wt')
  WRITE(98,*) 
  WRITE(98,*) ' RKStep must be (Fixed) or (Adapt) !!'
  CLOSE(98)
  STOP

END SELECT

1 CONTINUE
WRITE(99,*)' '
WRITE(99,*)'   Iterations made= ',it,'  (requested= ',Nitera,')'
WRITE(99,*)' '

CLOSE(2)

CALL Save_Results(Ulim,Ucrit,Qxi,Qxf,Qzi,Qzf,InitialLevel,Iw0,&
     RatioNf,Uf,RTfTe,RatioNb,Ub,RTbTe,RTeTi,G,ADResCondition,&
     LdecayResCondition,SdecayResCondition,Lemis,Ldecay,Lscat,&
     Semis,Sdecay,Sscat,AdiRoutine,BoundaryCondition,RKStep,Tau,&
     Dqx,Dqz,Dux,Duz,Anorm0,Epart0,Ewave0,EppEw0,Rn,Rp,Rw,Rs)

CALL Output("Ai ")
CALL Output("Dij")
CALL Output("Fe ")
CALL Output("IL ")
CALL Output("IS ")
CALL Output("Fe1")
CALL Output("IL1")
CALL Output("IS1")

END PROGRAM WT_2D


!---------------------------------------------------------------------
! Functions and Subroutines:
!---------------------------------------------------------------------

SUBROUTINE Init_Wave(Ulim,Qxi,Qxf,Qzi,Qzf,InitialLevel,Iw0,RatioNf,Uf,&
     RTfTe,RatioNb,Ub,RTbTe,Tau,Dqx,Dqz,Dux,Duz,Anorm,Epart,Ewave,EppEw)
USE Common_Arrays
USE Common_Params
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL :: Ulim,Qxi,Qxf,Qzi,Qzf,Iw0
REAL :: RatioNf,Uf,RTfTe
REAL :: RatioNb,Ub,RTbTe
REAL :: RMiMe,AA,Geff,U0
REAL, INTENT(out) :: Tau
REAL :: Dqx,Dqz,Dux,Duz 
REAL :: Ux,Uz,U,U2,Qx,Qz,Q 
REAL :: ZL,Zlq,Zs,Zsq,Zsoq,Aux
REAL :: Anorm,Epart,Ewave,EppEw
REAL, DIMENSION(nux,nuz) :: Dfdux,Dfduz
INTEGER :: i,j,k,nu2
CHARACTER(LEN=6) :: InitialLevel

CALL Definitions(RMiMe,AA,Geff)
Tau= 0.      ! Initializes the time counter.
Dqx= (Qxf-Qxi)/(nqx-1)
DO i= 1,nqx
 VQx(i)= Qxi+(i-1)*Dqx
END DO
VQx(nqx)= Qxf
Dqz= (Qzf-Qzi)/(nqz-1)
DO i= 1,nqz
 VQz(i)= Qzi+(i-1)*Dqz
END DO
VQz(nqz)= Qzf

! Initializes the spectrum of L waves and S waves:
SELECT CASE(InitialLevel)

 CASE("Choose")
  ILp= Iw0
  ILm= Iw0
  ISp= Iw0
  ISm= Iw0

 CASE("Auto  ")
  DO i= 1,nqx
   Qx= VQx(i)
   DO j= 1,nqz
    Qz= VQz(j)
    Q= SQRT(Qx**2+Qz**2)
    Zlq= ZL(Qx,Qz)
    Zsq= ZS(Qx,Qz)
    Zsoq= AA/SQRT(1.+Q**2/2.)
    ILp(i,j)= Geff/2./Zlq**2
    ILm(i,j)= Geff/2./Zlq**2
    Aux= (EXP(-(Zsoq**2))+SQRT(RMiMe*RTeTi)*EXP(-RMiMe*RTeTi*(Zsoq)**2)) &
      / (EXP(-(Zsoq)**2)+RTeTi*SQRT(RMiMe*RTeTi)*EXP(-RMiMe*RTeTi*(Zsoq)**2))
    ISp(i,j)= Geff/2./Zlq/Zsq * Aux
    ISm(i,j)= Geff/2./Zlq/Zsq * Aux
   END DO
  END DO

 CASE DEFAULT
  OPEN(98,FILE='Warning_Init_Wave.wt')
  WRITE(98,*) ' InitialLevel= ',InitialLevel
  WRITE(98,*) ' InitialLevel must be (Choose) or (Auto  ) !!'
  CLOSE(98)
  STOP

END SELECT

Dux= (Ulim-0.)/(nux-1)
DO i= 1,nux
 VUx(i)= 0.+(i-1)*Dux
END DO
Duz= (Ulim-(-Ulim))/(nuz-1)
DO i= 1,nuz
 VUz(i)= -Ulim+(i-1)*Duz
END DO

! Symmetrization of the array for Uz:
nu2= (nuz-1)/2+1
DO i= 1,nu2
 VUz(nuz+1-i)= -VUz(i)
END DO
VUz(nu2)= 0.

! Initializes the electron and ion distribution functions:
U0= -(RatioNf*Uf+RatioNb*Ub)/(1.-RatioNf-RatioNb)
DO k= 1,nuz
 Uz= VUz(k)
 DO i= 1,nux
  Ux= VUx(i)
  U2= Ux**2+Uz**2
  U= SQRT(U2)
  Fe(i,k)= (1.-RatioNf-RatioNb)*EXP(-(Ux**2+(Uz-U0)**2))/(Pi) &
     + RatioNf/(Pi*RTfTe)*EXP(-((Ux**2))/3 -(((Uz-Uf)**2))/1) &
     + RatioNb/(Pi*RTbTe)*EXP(-(Ux**2+(Uz-Ub)**2)/RTbTe)
  Fi(i,k)= (RMiMe*RTeTi)*EXP(-RMiMe*RTeTi*U2)/(Pi)
 END DO
END DO
CALL Fnorm(Anorm)
CALL Energy(Epart,Ewave,EppEw)
Fe0= Fe

! Initializes the boundary conditions for the eletron distribution: 
CALL Derivxy5p(nux,nuz,VUx,VUz,Fe0,Dfdux,Dfduz)
DO k= 1,nuz
 Fepxp(k)= Dfdux(nux,k)
END DO
DO i= 1,nux
 Fepzp(i)= Dfduz(i,nuz)
 Fepzm(i)= Dfduz(i,1)
END DO

RETURN
END SUBROUTINE Init_Wave

SUBROUTINE Save_Results(Ulim,Ucrit,Qxi,Qxf,Qzi,Qzf,InitialLevel,Iw0,&
     RatioNf,Uf,RTfTe,RatioNb,Ub,RTbTe,RTeTi,G,ADResCondition,&
     LdecayResCondition,SdecayResCondition,Lemis,Ldecay,Lscat,&
     Semis,Sdecay,Sscat,AdiRoutine,BoundaryCondition,RKStep,Tau,&
     Dqx,Dqz,Dux,Duz,Anorm0,Epart0,Ewave0,EppEw0,Rn,Rp,Rw,Rs)
USE Common_Arrays
!USE Common_Params
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL :: Ulim,Ucrit,Qxi,Qxf,Qzi,Qzf,Iw0
REAL :: RatioNf,Uf,RTfTe
REAL :: RatioNb,Ub,RTbTe
REAL :: RTeTi,G
REAL :: Tau
REAL :: Dqx,Dqz,Dux,Duz 
REAL :: Anorm0,Ewave0,Epart0,EppEw0
REAL :: Rn,Rp,Rw,Rs
INTEGER :: i,k,iflag
CHARACTER(LEN=6) :: InitialLevel
CHARACTER(LEN=6) :: ADResCondition,LdecayResCondition,SdecayResCondition
CHARACTER(LEN=3) :: Lemis,Ldecay,Lscat,Semis,Sdecay,Sscat
CHARACTER(LEN=5) :: AdiRoutine
CHARACTER(LEN=10) :: BoundaryCondition
CHARACTER(LEN=5) :: RKStep

Iflag= 1
OPEN(1,FILE='Out.wt')
WRITE(1,*) Iflag
WRITE(1,*) Ulim
WRITE(1,*) Ucrit
WRITE(1,*) Qxi 
WRITE(1,*) Qxf 
WRITE(1,*) Qzi 
WRITE(1,*) Qzf 
WRITE(1,*) InitialLevel
WRITE(1,*) Iw0 
WRITE(1,*) RatioNf 
WRITE(1,*) Uf 
WRITE(1,*) RTfTe 
WRITE(1,*) RatioNb 
WRITE(1,*) Ub 
WRITE(1,*) RTbTe 
WRITE(1,*) RTeTi 
WRITE(1,*) G
WRITE(1,*) ADResCondition
WRITE(1,*) Lemis
WRITE(1,*) Ldecay
WRITE(1,*) Lscat
WRITE(1,*) Semis
WRITE(1,*) Sdecay
WRITE(1,*) Sscat
WRITE(1,*) LdecayResCondition
WRITE(1,*) SdecayResCondition
WRITE(1,*) AdiRoutine
WRITE(1,*) BoundaryCondition
WRITE(1,*) RKStep

DO i= 1,nux
 WRITE(1,*) VUx(i)
END DO
DO k= 1,nuz
 WRITE(1,*) VUz(k)
END DO
DO i= 1,nux
 DO k= 1,nuz
  WRITE(1,*) Fe0(i,k)
 END DO
END DO
DO i= 1,nux
 DO k= 1,nuz
  WRITE(1,*) Fe(i,k)
 END DO
END DO
DO i= 1,nux
 DO k= 1,nuz
  WRITE(1,*) Fi(i,k)
 END DO
END DO
DO i= 1,nqx
 WRITE(1,*) VQx(i)
END DO
DO k= 1,nqz
 WRITE(1,*) VQz(k)
END DO
DO i= 1,nqx
 DO k= 1,nqz
  WRITE(1,*) ILp(i,k)
 END DO
END DO
DO i= 1,nqx
 DO k= 1,nqz
  WRITE(1,*) ILm(i,k)
 END DO
END DO
DO i= 1,nqx
 DO k= 1,nqz
  WRITE(1,*) ISp(i,k)
 END DO
END DO
DO i= 1,nqx
 DO k= 1,nqz
  WRITE(1,*) ISm(i,k)
 END DO
END DO
DO i= 1,nux
 WRITE(1,*) Fepzp(i),Fepzm(i)
END DO
DO k= 1,nuz
 WRITE(1,*) Fepxp(k)
END DO
WRITE(1,*) nux
WRITE(1,*) nuz
WRITE(1,*) nqx
WRITE(1,*) nqz
WRITE(1,*) Dqx
WRITE(1,*) Dqz
WRITE(1,*) Dux
WRITE(1,*) Duz
WRITE(1,*) Anorm0
WRITE(1,*) Epart0
WRITE(1,*) Ewave0
WRITE(1,*) EppEw0
WRITE(1,*) Rn
WRITE(1,*) Rp
WRITE(1,*) Rw
WRITE(1,*) Rs
WRITE(1,*) Tau
CLOSE(1)

RETURN
END SUBROUTINE Save_Results 

SUBROUTINE Read_Results(Tau,Dqx,Dqz,Dux,Duz,Anorm0,Epart0,Ewave0,EppEw0,&
     Rn,Rp,Rw,Rs)
USE Common_Arrays
USE Common_Params
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL :: Tau
REAL :: Dqx,Dqz,Dux,Duz
REAL :: Anorm0,Ewave0,Epart0,EppEw0
REAL :: Rn,Rp,Rw,Rs
INTEGER :: nuxback,nuzback,nqxback,nqzback
INTEGER :: i,k

DO i= 1,nux
 READ(1,*) VUx(i)
END DO
DO k= 1,nuz
 READ(1,*) VUz(k)
END DO
DO i= 1,nux
 DO k= 1,nuz
  READ(1,*) Fe0(i,k)
 END DO
END DO
DO i= 1,nux
 DO k= 1,nuz
  READ(1,*) Fe(i,k)
 END DO
END DO
DO i= 1,nux
 DO k= 1,nuz
  READ(1,*) Fi(i,k)
 END DO
END DO
DO i= 1,nqx
 READ(1,*) VQx(i)
END DO
DO k= 1,nqz
 READ(1,*) VQz(k)
END DO
DO i= 1,nqx
 DO k= 1,nqz
  READ(1,*) ILp(i,k)
 END DO
END DO
DO i= 1,nqx
 DO k= 1,nqz
  READ(1,*) ILm(i,k)
 END DO
END DO
DO i= 1,nqx
 DO k= 1,nqz
  READ(1,*) ISp(i,k)
 END DO
END DO
DO i= 1,nqx
 DO k= 1,nqz
  READ(1,*) ISm(i,k)
 END DO
END DO
DO i= 1,nux
 READ(1,*) Fepzp(i),Fepzm(i)
END DO
DO k= 1,nuz
 READ(1,*) Fepxp(k)
END DO
READ(1,*) nuxback
READ(1,*) nuzback
READ(1,*) nqxback
READ(1,*) nqzback
READ(1,*) Dqx
READ(1,*) Dqz
READ(1,*) Dux
READ(1,*) Duz
READ(1,*) Anorm0
READ(1,*) Epart0
READ(1,*) Ewave0
READ(1,*) EppEw0
READ(1,*) Rn
READ(1,*) Rp
READ(1,*) Rw
READ(1,*) Rs
READ(1,*) Tau
CLOSE(1)

IF(nuxback/=nux .OR. nuzback/=nuz .OR. nqxback/=nqx .OR. nqzback/=nqz) THEN
 OPEN(98,FILE='Warning_Read_Results.wt')
 WRITE(98,*) ' nuxback= ',nuxback,' nux= ',nux
 WRITE(98,*) ' nuzback= ',nuzback,' nuz= ',nuz
 WRITE(98,*) ' nqxback= ',nqxback,' nqx= ',nqx
 WRITE(98,*) ' nqzback= ',nqzback,' nqz= ',nqz
 WRITE(98,*) ' The quantities nux, nuz, nqx, and nqz can not be modified !!'
 CLOSE(98)
 STOP
ELSE
END IF

RETURN
END SUBROUTINE Read_Results 


SUBROUTINE Coef_A
USE Common_Arrays
USE Common_Params
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL :: RMiMe,AA,Geff
REAL :: Ux,Uz,Qx,Qx2
REAL :: AbsUz
REAL :: Qp1,Qm1,Qp12,Qm12
REAL :: Qp2,Qm2,Qp22,Qm22
REAL :: Aux,Auxp,Auxm
REAL :: Auxx,Auxz
REAL, DIMENSION(nqx) :: VintAx,VintAz
INTEGER :: i,l,m,sigma

CALL Definitions(RMiMe,AA,Geff)
!   Initialization of the Ai:
Ax= 0.
Az= 0.

SELECT CASE(ADResCondition)

 CASE("Exact ")
  DO m= 1,nuz
   Uz= VUz(m)
   DO l= 1,nux
    Ux= VUx(l)
    DO sigma= -1,1,2
     VintAx= 0.
     VintAz= 0.
     DO i= 1,nqx
      Qx= VQx(i)
      Qx2= Qx**2
      Aux= Uz**2-3.+3.*sigma*Qx*Ux-(9./4.)*Qx**2
      IF (Aux>0.) THEN
       Auxp= SQRT(Aux)
       Qp1= sigma*2./3. * (Uz+Auxp)
       IF (Qp1>0.) THEN
        Qp12= Qp1**2
        VintAx(i)= VintAx(i)+(Qx*Ux+Qp1*Uz)/(Qx2+Qp12) * Qx
        VintAz(i)= VintAz(i)+(Qx*Ux+Qp1*Uz)/(Qx2+Qp12) * Qp1
       ELSE
        VintAx(i)= VintAx(i)+0.
        VintAz(i)= VintAz(i)+0.
       END IF
       Qm1= sigma*2./3. * (Uz-Auxp)
       IF (Qm1>0.) THEN
        Qm12= Qm1**2
        VintAx(i)= VintAx(i)+(Qx*Ux+Qm1*Uz)/(Qx2+Qm12) * Qx
        VintAz(i)= VintAz(i)+(Qx*Ux+Qm1*Uz)/(Qx2+Qm12) * Qm1
       ELSE
        VintAx(i)= VintAx(i)+0.
        VintAz(i)= VintAz(i)+0.
       END IF
       VintAx(i)= VintAx(i)/Auxp
       VintAz(i)= VintAz(i)/Auxp
      ELSE
       VintAx(i)= 0.
       VintAz(i)= 0.
      END IF
     END DO
     CALL Simpson(VQx,VintAx,nqx,Auxx)
     CALL Simpson(VQx,VintAz,nqx,Auxz)
     Ax(l,m)= Ax(l,m)+2.*Geff*Auxx
     Az(l,m)= Az(l,m)+2.*Geff*Auxz
    END DO
    DO sigma= -1,1,2
     VintAx= 0.
     VintAz= 0.
     DO i= 1,nqx
      Qx= VQx(i)
      Qx2= Qx**2
      Aux= Uz**2-3.-3.*sigma*Qx*Ux-(9./4.)*Qx**2
      IF (Aux>0.) THEN
       Auxm= SQRT(Aux)
       Qp2= sigma*2./3. * (Uz+Auxm)
       IF (Qp2>0.) THEN
        Qp22= Qp2**2
        VintAx(i)= VintAx(i)+(Qx*Ux+Qp2*Uz)/(Qx2+Qp22) * Qx
        VintAz(i)= VintAz(i)+(Qx*Ux+Qp2*Uz)/(Qx2+Qp22) * Qp2
       ELSE
        VintAx(i)= VintAx(i)+0.
        VintAz(i)= VintAz(i)+0.
       END IF
       Qm2= sigma*2./3. * (Uz-Auxm)
       IF (Qm2>0.) THEN
        Qm22= Qm2**2
        VintAx(i)= VintAx(i)+(Qx*Ux+Qm2*Uz)/(Qx2+Qm22) * Qx
        VintAz(i)= VintAz(i)+(Qx*Ux+Qm2*Uz)/(Qx2+Qm22) * Qm2
       ELSE
        VintAx(i)= VintAx(i)+0.
        VintAz(i)= VintAz(i)+0.
       END IF
       VintAx(i)= VintAx(i)/Auxm
       VintAz(i)= VintAz(i)/Auxm
      ELSE
       VintAx(i)= 0.
       VintAz(i)= 0.
      END IF
     END DO
     CALL Simpson(VQx,VintAx,nqx,Auxx)
     CALL Simpson(VQx,VintAz,nqx,Auxz)
     Ax(l,m)= Ax(l,m)+2.*Geff*Auxx
     Az(l,m)= Az(l,m)+2.*Geff*Auxz
    END DO
   END DO
  END DO

 CASE("Approx")
  DO m= 1,nuz
   Uz= VUz(m)
   AbsUz= ABS(Uz)
   IF (AbsUz>1.E-6) THEN
    DO l= 1,nux
     Ux= VUx(l)
     DO sigma= -1,1,2
      VintAx= 0.
      VintAz= 0.
      DO i= 1,nqx
       Qx= VQx(i)
       Qx2= Qx**2
       Qp1= (sigma-Qx*Ux)/Uz
       IF (Qp1>0.) THEN
        Qp12= Qp1**2
        VintAx(i)= VintAx(i)+(Qx*Ux+Qp1*Uz)/(Qx2+Qp12) * Qx
        VintAz(i)= VintAz(i)+(Qx*Ux+Qp1*Uz)/(Qx2+Qp12) * Qp1
       ELSE
        VintAx(i)= VintAx(i)+0.
        VintAz(i)= VintAz(i)+0.
       END IF
       Qp2= (sigma+Qx*Ux)/Uz
       IF (Qp2>0.) THEN
        Qp22= Qp2**2
        VintAx(i)= VintAx(i)+(-Qx*Ux+Qp2*Uz)/(Qx2+Qp22) * (-Qx)
        VintAz(i)= VintAz(i)+(-Qx*Ux+Qp2*Uz)/(Qx2+Qp22) * Qp2
       ELSE
        VintAx(i)= VintAx(i)+0.
        VintAz(i)= VintAz(i)+0.
       END IF
      END DO
      CALL Simpson(VQx,VintAx,nqx,Auxx)
      CALL Simpson(VQx,VintAz,nqx,Auxz)
      Ax(l,m)= Ax(l,m)+(2.*Geff/AbsUz)*Auxx
      Az(l,m)= Az(l,m)+(2.*Geff/AbsUz)*Auxz
     END DO
    END DO
   ELSE
    Ax(l,m)= 0.
    Az(l,m)= 0.
   END IF
  END DO

 CASE DEFAULT
  OPEN(98,FILE='Warning_Coef_A.wt')
  WRITE(98,*) ' ADResCondition= ',ADResCondition
  WRITE(98,*) ' ADResCondition must be (Exact ) or (Approx) !!'
  CLOSE(98)
  STOP

END SELECT

RETURN
END SUBROUTINE Coef_A

SUBROUTINE Coef_D
USE Common_Arrays
USE Common_Params
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL :: Ux,Uz,Qx,Qx2
REAL :: AbsUz,EsLqp1,EsLqp2,EsLqm1,EsLqm2,EsLq1,EsLq2
REAL :: Qp1,Qm1,Qp12,Qm12
REAL :: Qp2,Qm2,Qp22,Qm22
REAL :: Aux,Auxp,Auxm
REAL :: Auxxx,Auxxz,Auxzz
REAL, DIMENSION(nqz) :: Vauxq
REAL, DIMENSION(nqx) :: VintDxx,VintDxz,VintDzz
INTEGER :: i,l,m,sigma

!   Initialization of the Dij:
Dxx= 0.
Dxz= 0.
Dzx= 0.
Dzz= 0.

SELECT CASE(ADResCondition)

 CASE("Exact ")
  DO m= 1,nuz
   Uz= VUz(m)
   DO l= 1,nux
    Ux= VUx(l)
    DO sigma= -1,1,2
     VintDxx= 0.
     VintDxz= 0.
     VintDzz= 0.
     DO i= 1,nqx
      Qx= VQx(i)
      Qx2= Qx**2
      Aux= Uz**2-3.+3.*sigma*Qx*Ux-(9./4.)*Qx**2
      IF (Aux>0.) THEN
       Auxp= SQRT(Aux)
       Qp1= sigma*2./3. * (Uz+Auxp)
       IF (Qp1>0.) THEN
        Qp12= Qp1**2
        IF (sigma==1) THEN
         Vauxq(:)= ILp(i,:)
         CALL Aitp1d(nqz,VQz,Vauxq,Qp1,EsLqp1)
        ELSE
         Vauxq(:)= ILm(i,:)
         CALL Aitp1d(nqz,VQz,Vauxq,Qp1,EsLqp1)
        END IF
        VintDxx(i)= VintDxx(i)+Qx2/(Qx2+Qp12) * EsLqp1
        VintDxz(i)= VintDxz(i)+Qx*Qp1/(Qx2+Qp12) * EsLqp1
        VintDzz(i)= VintDzz(i)+Qp12/(Qx2+Qp12) * EsLqp1
       ELSE
        VintDxx(i)= VintDxx(i)+0.
        VintDxz(i)= VintDxz(i)+0.
        VintDzz(i)= VintDzz(i)+0.
       END IF
       Qm1= sigma*2./3. * (Uz-Auxp)
       IF (Qm1>0.) THEN
        Qm12= Qm1**2
        IF (sigma==1) THEN
         Vauxq(:)= ILp(i,:)
         CALL Aitp1d(nqz,VQz,Vauxq,Qm1,EsLqm1)
        ELSE
         Vauxq(:)= ILm(i,:)
         CALL Aitp1d(nqz,VQz,Vauxq,Qm1,EsLqm1)
        END IF
        VintDxx(i)= VintDxx(i)+Qx2/(Qx2+Qm12) * EsLqm1
        VintDxz(i)= VintDxz(i)+Qx*Qm1/(Qx2+Qm12) * EsLqm1
        VintDzz(i)= VintDzz(i)+Qm12/(Qx2+Qm12) * EsLqm1
       ELSE
        VintDxx(i)= VintDxx(i)+0.
        VintDxz(i)= VintDxz(i)+0.
        VintDzz(i)= VintDzz(i)+0.
       END IF
       VintDxx(i)= VintDxx(i)/Auxp
       VintDxz(i)= VintDxz(i)/Auxp
       VintDzz(i)= VintDzz(i)/Auxp
      ELSE
       VintDxx(i)= 0.
       VintDxz(i)= 0.
       VintDzz(i)= 0.
      END IF
     END DO
     CALL Simpson(VQx,VintDxx,nqx,Auxxx)
     CALL Simpson(VQx,VintDxz,nqx,Auxxz)
     CALL Simpson(VQx,VintDzz,nqx,Auxzz)
     Dxx(l,m)= Dxx(l,m)+2.*Auxxx
     Dxz(l,m)= Dxz(l,m)+2.*Auxxz
     Dzx(l,m)= Dzx(l,m)+2.*Auxxz
     Dzz(l,m)= Dzz(l,m)+2.*Auxzz
    END DO
    DO sigma= -1,1,2
     VintDxx= 0.
     VintDxz= 0.
     VintDzz= 0.
     DO i= 1,nqx
      Qx= VQx(i)
      Qx2= Qx**2
      Aux= Uz**2-3.-3.*sigma*Qx*Ux-(9./4.)*Qx**2
      IF (Aux>0.) THEN
       Auxm= SQRT(Aux)
       Qp2= sigma*2./3. * (Uz+Auxm)
       IF (Qp2>0.) THEN
        Qp22= Qp2**2
        IF (sigma==1) THEN
         Vauxq(:)= ILp(i,:)
         CALL Aitp1d(nqz,VQz,Vauxq,Qp2,EsLqp2)
        ELSE
         Vauxq(:)= ILm(i,:)
         CALL Aitp1d(nqz,VQz,Vauxq,Qp2,EsLqp2)
        END IF
        VintDxx(i)= VintDxx(i)+Qx2/(Qx2+Qp22) * EsLqp2
        VintDxz(i)= VintDxz(i)+Qx*Qp2/(Qx2+Qp22) * EsLqp2
        VintDzz(i)= VintDzz(i)+Qp22/(Qx2+Qp22) * EsLqp2
       ELSE
        VintDxx(i)= VintDxx(i)+0.
        VintDxz(i)= VintDxz(i)+0.
        VintDzz(i)= VintDzz(i)+0.
       END IF
       Qm2= sigma*2./3. * (Uz-Auxm)
       IF (Qm2>0.) THEN
        Qm22= Qm2**2
        IF (sigma==1) THEN
         Vauxq(:)= ILp(i,:)
         CALL Aitp1d(nqz,VQz,Vauxq,Qm2,EsLqm2)
        ELSE
         Vauxq(:)= ILm(i,:)
         CALL Aitp1d(nqz,VQz,Vauxq,Qm2,EsLqm2)
        END IF
        VintDxx(i)= VintDxx(i)+Qx2/(Qx2+Qm22) * EsLqm2
        VintDxz(i)= VintDxz(i)+Qx*Qm2/(Qx2+Qm22) * EsLqm2
        VintDzz(i)= VintDzz(i)+Qm22/(Qx2+Qm22) * EsLqm2
       ELSE
        VintDxx(i)= VintDxx(i)+0.
        VintDxz(i)= VintDxz(i)+0.
        VintDzz(i)= VintDzz(i)+0.
       END IF
       VintDxx(i)= VintDxx(i)/Auxm
       VintDxz(i)= VintDxz(i)/Auxm
       VintDzz(i)= VintDzz(i)/Auxm
      ELSE
       VintDxx(i)= 0.
       VintDxz(i)= 0.
       VintDzz(i)= 0.
      END IF
     END DO
     CALL Simpson(VQx,VintDxx,nqx,Auxxx)
     CALL Simpson(VQx,VintDxz,nqx,Auxxz)
     CALL Simpson(VQx,VintDzz,nqx,Auxzz)
     Dxx(l,m)= Dxx(l,m)+2.*Auxxx
     Dxz(l,m)= Dxz(l,m)+2.*Auxxz
     Dzx(l,m)= Dzx(l,m)+2.*Auxxz
     Dzz(l,m)= Dzz(l,m)+2.*Auxzz
    END DO
   END DO
  END DO

 CASE("Approx")
  DO m= 1,nuz
   Uz= VUz(m)
   AbsUz= ABS(Uz)
   IF (AbsUz>1.E-6) THEN
    DO l= 1,nux
     Ux= VUx(l)
     DO sigma= -1,1,2
      VintDxx= 0.
      VintDxz= 0.
      VintDzz= 0.
      DO i= 1,nqx
       Qx= VQx(i)
       Qx2= Qx**2
       Qp1= (sigma-Qx*Ux)/Uz
       IF (Qp1>0.) THEN
        Qp12= Qp1**2
        IF (sigma==1) THEN
         Vauxq(:)= ILp(i,:)
         CALL Aitp1d(nqz,VQz,Vauxq,Qp1,EsLq1)
        ELSE
         Vauxq(:)= ILm(i,:)
         CALL Aitp1d(nqz,VQz,Vauxq,Qp1,EsLq1)
        END IF
        VintDxx(i)= VintDxx(i)+Qx2/(Qx2+Qp12) * EsLq1
        VintDxz(i)= VintDxz(i)+Qx*Qp1/(Qx2+Qp12) * EsLq1
        VintDzz(i)= VintDzz(i)+Qp12/(Qx2+Qp12) * EsLq1
       ELSE
        VintDxx(i)= VintDxx(i)+0.
        VintDxz(i)= VintDxz(i)+0.
        VintDzz(i)= VintDzz(i)+0.
       END IF
       Qp2= (sigma+Qx*Ux)/Uz
       IF (Qp2>0.) THEN
        Qp22= Qp2**2
        IF (sigma==1) THEN
         Vauxq(:)= ILp(i,:)
         CALL Aitp1d(nqz,VQz,Vauxq,Qp2,EsLq2)
        ELSE
         Vauxq(:)= ILm(i,:)
         CALL Aitp1d(nqz,VQz,Vauxq,Qp2,EsLq2)
        END IF
        VintDxx(i)= VintDxx(i)+Qx2/(Qx2+Qp22) * EsLq2
        VintDxz(i)= VintDxz(i)+(-Qx*Qp2)/(Qx2+Qp22) * EsLq2
        VintDzz(i)= VintDzz(i)+Qp22/(Qx2+Qp22) * EsLq2
       ELSE
        VintDxx(i)= VintDxx(i)+0.
        VintDxz(i)= VintDxz(i)+0.
        VintDzz(i)= VintDzz(i)+0.
       END IF
      END DO
      CALL Simpson(VQx,VintDxx,nqx,Auxxx)
      CALL Simpson(VQx,VintDxz,nqx,Auxxz)
      CALL Simpson(VQx,VintDzz,nqx,Auxzz)
      Dxx(l,m)= Dxx(l,m)+(2./AbsUz)*Auxxx
      Dxz(l,m)= Dxz(l,m)+(2./AbsUz)*Auxxz
      Dzx(l,m)= Dzx(l,m)+(2./AbsUz)*Auxxz
      Dzz(l,m)= Dzz(l,m)+(2./AbsUz)*Auxzz
     END DO
    END DO
   ELSE
    Dxx(l,m)= 0.
    Dxz(l,m)= 0.
    Dzx(l,m)= 0.
    Dzz(l,m)= 0.
   END IF
  END DO

 CASE DEFAULT
  OPEN(98,FILE='Warning_Coef_D.wt')
  WRITE(98,*) ' ADResCondition= ',ADResCondition
  WRITE(98,*) ' ADResCondition must be (Exact ) or (Approx) !!'
  CLOSE(98)
  STOP

END SELECT

RETURN
END SUBROUTINE Coef_D

SUBROUTINE Coef_Lwave(Qx,Qz,iqx,kqz,sigma,Dfdux,Dfduz,Iwave,Iwavem,Coef)
USE Common_Arrays
USE Common_Params
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL :: RMiMe,AA,Geff
REAL :: Ux,Qx,Qz
REAL :: Qx2,Qz2,Q2,Ures
REAL :: ZL,Zlq,Zlqp,Zlqdif
REAL :: Exp_Sqr_Lwave
REAL :: Aux,Aux0,Aux1,Aux2,Aux3,Aux4,Aux5,Aux6
REAL :: Coef
REAL :: CoefE,CoefD,CoefS
REAL :: Qxp,Qzp,Qxp2,Qzp2,Qzpden
REAL :: Qd1
REAL :: Qxdif,Qzdif,Iqp,Iqdif
REAL, DIMENSION(nux) :: Vintux
REAL, DIMENSION(nuz) :: Vauxuz
REAL, DIMENSION(nux,nuz) :: Dfdux,Dfduz
REAL, DIMENSION(nqx,nqz) :: Iwave,Iwavem
REAL, DIMENSION(nqx) :: Vintqx
REAL, DIMENSION(nqz) :: Vintqz
REAL, DIMENSION(nqz) :: VauxqA,VauxqB
INTEGER :: i,k,l,sigma
INTEGER :: iqx,kqz
CHARACTER(LEN=2) :: SignCase 

CALL Definitions(RMiMe,AA,Geff)
Qx2= Qx**2
Qz2= Qz**2
Q2= Qx2+Qz2
Zlq= ZL(Qx,Qz)

IF(Lemis=="yes") THEN
 ! Contribution due to spontaneous and induced emission:
 DO l= 1,nux
  Ux= VUx(l)
  Vauxuz(:)= Geff*Fe(l,:)-(sigma*Zlq)*Iwave(iqx,kqz)*Qx*Dfdux(l,:) &
             +(sigma*Zlq)*Iwave(iqx,kqz)*Qz*Dfduz(l,:)
  Ures= (sigma*Zlq+Qx*Ux)/Qz
  CALL Aitp1d(nuz,VUz,Vauxuz,Ures,Aux)
  Vintux(l)= Aux
  Vauxuz(:)= Geff*Fe(l,:)+(sigma*Zlq)*Iwave(iqx,kqz)*Qx*Dfdux(l,:) &
             +(sigma*Zlq)*Iwave(iqx,kqz)*Qz*Dfduz(l,:)
  Ures= (sigma*Zlq-Qx*Ux)/Qz
  CALL Aitp1d(nuz,VUz,Vauxuz,Ures,Aux)
  Vintux(l)= Vintux(l) + Aux
 END DO
 CALL Simpson(VUx,Vintux,nux,Aux)
 CoefE= (Pi/Q2/(ABS(Qz))) * Aux
ELSE
 CoefE= 0.
END IF

If(Ldecay=="yes") THEN
 ! Contribution due to spontaneous and induced decay:
 VauxqA= 0.   ! Initializes the integrand of the Qzp integrals
 VauxqB= 0.

 SELECT CASE(LdecayResCondition)

  CASE("Exact ")
   OPEN(98,FILE='Warning_Coef_Lwave.wt')
   WRITE(98,*) ' LdecayResCondition= ',LdecayResCondition
   WRITE(98,*) ' The case (Exact ) is not yet developed for L waves !!'
   CLOSE(98)
   STOP

  CASE("Approx")
   DO i= 1,nqx
    Qxp= VQx(i)
    Aux= Qx**2+Qz**2-Qxp**2
    IF (Aux>0.) THEN
     Qd1= SQRT(Aux)
     Qzp= Qd1
     SignCase= "pm"
     Qxdif= Qx+Qxp
     Qzdif= Qz-Qzp
     Aux0= (Qx*Qxp-Qz*Qzp)**2
     CALL Aux_Coef_Lwave(sigma,SignCase,i,Qx,Qz,Qxp,Qzp,Qxdif,Qzdif,&
	Aux0,Aux1,Zlqp,Zlqdif,Iqp,Iqdif)
     VauxqA(i)= VauxqA(i)+Aux1*Iqdif*Iqp
     VauxqB(i)= VauxqB(i)+Aux1*(Iqdif)
     SignCase= "mm"
     Qxdif= Qx-Qxp
     Qzdif= Qz-Qzp
     Aux0= (Qx*Qxp+Qz*Qzp)**2
     CALL Aux_Coef_Lwave(sigma,SignCase,i,Qx,Qz,Qxp,Qzp,Qxdif,Qzdif,&
	Aux0,Aux1,Zlqp,Zlqdif,Iqp,Iqdif)
     VauxqA(i)= VauxqA(i)+Aux1*Iqdif*Iqp
     VauxqB(i)= VauxqB(i)+Aux1*(Iqdif)
     SignCase= "pp"
     Qxdif= Qx+Qxp
     Qzdif= Qz+Qzp
     Aux0= (Qx*Qxp+Qz*Qzp)**2
     CALL Aux_Coef_Lwave(sigma,SignCase,i,Qx,Qz,Qxp,Qzp,Qxdif,Qzdif,&
	Aux0,Aux1,Zlqp,Zlqdif,Iqp,Iqdif)
     VauxqA(i)= VauxqA(i)+Aux1*Iqdif*Iqp
     VauxqB(i)= VauxqB(i)+Aux1*(Iqdif)
     SignCase= "mp"
     Qxdif= Qx-Qxp
     Qzdif= Qz+Qzp
     Aux0= (Qx*Qxp-Qz*Qzp)**2
     CALL Aux_Coef_Lwave(sigma,SignCase,i,Qx,Qz,Qxp,Qzp,Qxdif,Qzdif,&
	Aux0,Aux1,Zlqp,Zlqdif,Iqp,Iqdif)
     VauxqA(i)= VauxqA(i)+Aux1*Iqdif*Iqp
     VauxqB(i)= VauxqB(i)+Aux1*(Iqdif)
     IF(Qzp<0.01) THEN
      Qzpden= 0.01
     ELSE
      Qzpden= Qzp
     END IF
     VauxqA(i)= VauxqA(i)/Qzpden
     VauxqB(i)= VauxqB(i)/Qzpden
    ELSE
     VauxqA(i)= VauxqA(i) + 0.
     VauxqB(i)= VauxqB(i) + 0.
    END IF
   END DO
   CALL Simpson(VQx,VauxqA,nqx,Aux)
   CoefD= Aux
   CALL Simpson(VQx,VauxqB,nqx,Aux)
   CoefD= CoefD - Aux*Iwave(iqx,kqz)
   CoefD= AA * (2.*SQRT(2.)/3.) * (Zlq/Q2)**2 * CoefD

  CASE DEFAULT
   OPEN(98,FILE='Warning_Coef_Lwave.wt')
   WRITE(98,*) ' LdecayResCondition= ',LdecayResCondition
   WRITE(98,*) ' LdecayResCondition must be (Exact ) or (Approx) !!'
   CLOSE(98)
   STOP

 END SELECT

ELSE
 CoefD= 0.
END IF

IF(Lscat=="yes") THEN
 ! Contribution due to spontaneous and induced scattering:

  Vintqx= 0.
  DO i= 1,nqx
   Qxp= VQx(i)
   Qxp2= Qxp**2
   Vintqz= 0.
   DO k= 1,nqz
    Qzp= VQz(k)
    Qzp2= Qzp**2
    Zlqp= ZL(Qxp,Qzp)
    IF(Zlqp==Zlq) THEN
     Aux0= 0.
    ELSE
     Aux0= RMiMe*RTeTi*(Zlq-Zlqp)**2
    END IF
    Aux1= Qx-Qxp
    Aux2= Qz-Qzp
    Aux3= Exp_Sqr_Lwave(Aux0,Aux1,Aux2) 
    Aux1= Qx+Qxp
    Aux2= Qz+Qzp
    Aux4= Exp_Sqr_Lwave(Aux0,Aux1,Aux2) 
    Aux1= Qx+Qxp
    Aux2= Qz-Qzp
    Aux5= Exp_Sqr_Lwave(Aux0,Aux1,Aux2) 
    Aux1= Qx-Qxp
    Aux2= Qz+Qzp
    Aux6= Exp_Sqr_Lwave(Aux0,Aux1,Aux2) 
    Vintqz(k)= (Qx*Qxp+Qz*Qzp)**2/Q2/(Qxp2+Qzp2) &
     * ( (Geff*(Zlqp*Iwave(iqx,kqz)-Zlq*Iwave(i,k)) &
     + Iwave(i,k)*Iwave(iqx,kqz)*(Zlq-Zlqp)*2.*RTeTi) * Aux3 &
!     * EXP(-Aux0/((Qx-Qxp)**2+(Qz-Qzp)**2)) &
!     /SQRT((Qx-Qxp)**2+(Qz-Qzp)**2) &
     + (Geff*(Zlqp*Iwave(iqx,kqz)-Zlq*Iwavem(i,k)) &
     + Iwavem(i,k)*Iwave(iqx,kqz)*(Zlq-Zlqp)*2.*RTeTi) * Aux4 ) &
!     * EXP(-Aux0/((Qx+Qxp)**2+(Qz+Qzp)**2)) &
!     /SQRT((Qx+Qxp)**2+(Qz+Qzp)**2) ) &
     + (Qx*Qxp-Qz*Qzp)**2/Q2/(Qxp2+Qzp2) &
     * ( (Geff*(Zlqp*Iwave(iqx,kqz)-Zlq*Iwave(i,k)) &
     + Iwave(i,k)*Iwave(iqx,kqz)*(Zlq-Zlqp)*2.*RTeTi) * Aux5 &
!     * EXP(-Aux0/((Qx+Qxp)**2+(Qz-Qzp)**2)) &
!     /SQRT((Qx+Qxp)**2+(Qz-Qzp)**2) &
     + (Geff*(Zlqp*Iwave(iqx,kqz)-Zlq*Iwavem(i,k)) &
     + Iwavem(i,k)*Iwave(iqx,kqz)*(Zlq-Zlqp)*2.*RTeTi) * Aux6 )
!     * EXP(-Aux0/((Qx-Qxp)**2+(Qz+Qzp)**2)) &
!     /SQRT((Qx-Qxp)**2+(Qz+Qzp)**2) )
   END DO
   CALL Simpson(VQz,Vintqz,nqz,Aux)
   Vintqx(i)= Aux
  END DO
  CALL Simpson(VQx,Vintqx,nqx,Aux)
 CoefS= -(Zlq)*SQRT(RMiMe*RTeTi/Pi)*Aux
ELSE
 CoefS= 0.
END IF

Coef= CoefE+CoefD+CoefS
RETURN
END SUBROUTINE Coef_Lwave

SUBROUTINE Aux_Coef_Lwave(sigma,SignCase,iqxp,Qx,Qz,Qxp,Qzp,Qxdif,Qzdif,&
     Aux0,Aux1,Zlqp,Zlqdif,Iqp,Iqdif)
USE Common_Arrays
USE Common_Params
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL :: Qx,Qz,Qxp,Qzp
REAL :: ZL,Zlqp,Zlqdif
REAL :: Qxdif,Qzdif,Iqp,Iqdif
REAL :: Qxdifaux,Qzdifaux
REAL :: Aux0,Aux1
REAL, DIMENSION(nqz) :: Vauxq
REAL, DIMENSION(nqx,nqz) :: ILsigma,ISsigma
INTEGER :: sigma,iqxp
CHARACTER(LEN=2) :: SignCase 

SELECT CASE(sigma)
 CASE(1)
  SELECT CASE(SignCase)
   CASE("pm")
    ILsigma= ILp
    ISsigma= ISp+ISm
    IF((Qx**2+Qz**2+Qx*Qxp-Qz*Qzp)>0.) THEN
     Aux1= SQRT(Qx**2+Qz**2+Qx*Qxp-Qz*Qzp)*Aux0
    ELSE
     Aux1= 0.
    END IF
   CASE("mm")
    ILsigma= ILp
    ISsigma= ISp+ISm
    IF((Qx**2+Qz**2-Qx*Qxp-Qz*Qzp)>0.) THEN
     Aux1= SQRT(Qx**2+Qz**2-Qx*Qxp-Qz*Qzp)*Aux0
    ELSE
     Aux1= 0.
    END IF
   CASE("pp")
    ILsigma= ILm
    ISsigma= ISp+ISm
    IF((Qx**2+Qz**2+Qx*Qxp+Qz*Qzp)>0.) THEN
     Aux1= SQRT(Qx**2+Qz**2+Qx*Qxp+Qz*Qzp)*Aux0
    ELSE
     Aux1= 0.
    END IF
   CASE("mp")
    ILsigma= ILm
    ISsigma= ISp+ISm
    IF((Qx**2+Qz**2-Qx*Qxp+Qz*Qzp)>0.) THEN
     Aux1= SQRT(Qx**2+Qz**2-Qx*Qxp+Qz*Qzp)*Aux0
    ELSE
     Aux1= 0.
    END IF
   CASE DEFAULT
    OPEN(98,FILE='Warning_Aux_Coef_Lwave.wt')
    WRITE(98,*) ' SignCase= ',SignCase
    WRITE(98,*) ' SignCase must be pp, pm, mp, mm !! '
    CLOSE(98)
    STOP
  END SELECT
 CASE(-1)
  SELECT CASE(SignCase)
   CASE("pm")
    ILsigma= ILm
    ISsigma= ISp+ISm
    IF((Qx**2+Qz**2+Qx*Qxp-Qz*Qzp)>0.) THEN
     Aux1= SQRT(Qx**2+Qz**2+Qx*Qxp-Qz*Qzp)*Aux0
    ELSE
     Aux1= 0.
    END IF
   CASE("mm")
    ILsigma= ILm
    ISsigma= ISp+ISm
    IF((Qx**2+Qz**2-Qx*Qxp-Qz*Qzp)>0.) THEN
     Aux1= SQRT(Qx**2+Qz**2-Qx*Qxp-Qz*Qzp)*Aux0
    ELSE
     Aux1= 0.
    END IF
   CASE("pp")
    ILsigma= ILp
    ISsigma= ISp+ISm
    IF((Qx**2+Qz**2+Qx*Qxp+Qz*Qzp)>0.) THEN
     Aux1= SQRT(Qx**2+Qz**2+Qx*Qxp+Qz*Qzp)*Aux0
    ELSE
     Aux1= 0.
    END IF
   CASE("mp")
    ILsigma= ILp
    ISsigma= ISp+ISm
    IF((Qx**2+Qz**2-Qx*Qxp+Qz*Qzp)>0.) THEN
     Aux1= SQRT(Qx**2+Qz**2-Qx*Qxp+Qz*Qzp)*Aux0
    ELSE
     Aux1= 0.
    END IF
   CASE DEFAULT
    OPEN(98,FILE='Warning_Aux_Coef_Lwave.wt')
    WRITE(98,*) ' SignCase= ',SignCase
    WRITE(98,*) ' SignCase must be pp, pm, mp, mm !! '
    CLOSE(98)
    STOP
  END SELECT
 CASE DEFAULT
  OPEN(98,FILE='Warning_Aux_Coef_Lwave.wt')
  WRITE(98,*) ' sigma= ',sigma
  WRITE(98,*) ' sigma must be 1 or -1 !!'
  CLOSE(98)
  STOP
END SELECT
Zlqp= ZL(Qxp,Qzp)
Zlqdif= ZL(Qxdif,Qzdif)
IF (Qxdif>=0.) THEN
 Qxdifaux= Qxdif
ELSE
 Qxdifaux= -Qxdif
END IF
IF (Qzdif>=0.) THEN
 Qzdifaux= Qzdif
ELSE
 Qzdifaux= -Qzdif
END IF
CALL Aitp2d(nqx,nqz,VQx,VQz,ISsigma,Qxdifaux,Qzdifaux,Iqdif)
Vauxq(:)= ILsigma(iqxp,:)
CALL Aitp1d(nqz,VQz,Vauxq,Qzp,Iqp)
RETURN
END SUBROUTINE Aux_Coef_Lwave

REAL FUNCTION Exp_Sqr_Lwave(Aux0,A,B)
USE Math_Constants
IMPLICIT NONE
REAL, INTENT(in) :: Aux0,A,B
REAL :: C
C= A**2+B**2
IF(C<=EpsMin) THEN
 Exp_Sqr_Lwave= 0.
ELSE
 Exp_Sqr_Lwave= EXP(-Aux0/C)/SQRT(C)
END IF
RETURN
END FUNCTION Exp_Sqr_Lwave

SUBROUTINE Coef_Swave(Qx,Qz,iqx,kqz,sigma,Dfdux,Dfduz,Iwave,Coef)
USE Common_Arrays
USE Common_Params
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL :: RMiMe,AA,Geff
REAL :: Ux,Qx,Qz,Q,Q2
REAL :: Zlq,Zsq,Zsoq,Ures,Aux,Aux1
REAL :: ZL,ZS
REAL :: CoefE,CoefD,Coef
REAL :: Zlqp,Zlqdif
REAL :: Qdp,Qdm
REAL :: Muq
REAL :: Qxp,Qzp
REAL :: Qxdif,Qzdif,Iqp,Iqdif
REAL, DIMENSION(nux) :: Vintux
REAL, DIMENSION(nuz) :: Vauxuz
REAL, DIMENSION(nux,nuz) :: Dfdux,Dfduz
REAL, DIMENSION(nqx,nqz) :: Iwave
REAL, DIMENSION(nqx) :: VauxqA,VauxqB
INTEGER :: i,l,sigma
INTEGER :: iqx,kqz
CHARACTER(LEN=2) :: SignCase 

CALL Definitions(RMiMe,AA,Geff)
Q2= Qx**2+Qz**2
Q= SQRT(Q2)
Muq= Q**3*AA/2.
Zlq= ZL(Qx,Qz)
Zsq= ZS(Qx,Qz)
Zsoq= AA/SQRT(1.+Q**2/2.)

IF(Semis=="yes") THEN
 ! Contribution due to spontaneous and induced emission:
 DO l= 1,nux
  Ux= VUx(l)
  Vauxuz(:)= Geff*Fe(l,:)+(sigma*Zlq)*Iwave(iqx,kqz) &
     * (-Qx*Dfdux(l,:)+Qz*Dfduz(l,:))
  Ures= (sigma*Zsq+Qx*Ux)/Qz
  CALL Aitp1d(nuz,VUz,Vauxuz,Ures,Aux)
  Vintux(l)= Aux
  Vauxuz(:)= Geff*Fe(l,:)+(sigma*Zlq)*Iwave(iqx,kqz) &
     * (Qx*Dfdux(l,:)+Qz*Dfduz(l,:))
  Ures= (sigma*Zsq-Qx*Ux)/Qz
  CALL Aitp1d(nuz,VUz,Vauxuz,Ures,Aux)
  Vintux(l)= Vintux(l) + Aux
 END DO
 CALL Simpson(VUx,Vintux,nux,Aux)
 CoefE= Pi * Q*AA/2./(ABS(Qz)) * Aux
 CoefE= CoefE + SQRT(Pi)*AA/2.*(Geff-2.*RTeTi*Zlq*Zsq*Iwave(iqx,kqz)) &
     * SQRT(RMiMe*RTeTi)*EXP(-RMiMe*RTeTi*(Zsoq)**2)
ELSE
 CoefE= 0.
END IF

If(Sdecay=="yes") THEN
 ! Contribution due to spontaneous and induced decay:
 VauxqA= 0.   ! Initializes the integrand of the Qzp integrals
 VauxqB= 0.

 SELECT CASE(SdecayResCondition)

  CASE("Exact ")
   OPEN(98,FILE='Warning_Coef_Swave.wt')
   WRITE(98,*) ' SdecayResCondition= ',SdecayResCondition
   WRITE(98,*) ' The case (Exact ) is not yet developed for S waves !!'
   CLOSE(98)
   STOP

  CASE("Approx")
   DO i= 1,nqx
    Qxp= VQx(i)
    Qdp= (Q2+2.*Qx*Qxp)/(2.*Qz)
    Qdm= (Q2-2.*Qx*Qxp)/(2.*Qz)
    Qzp= Qdp
    Qxdif= Qx-Qxp
    Qzdif= Qz+Qzp
    SignCase= "pp"
    CALL Aux_Coef_Swave(SignCase,i,Qxp,Qzp,Qxdif,Qzdif,&
	Aux1,Zlqp,Zlqdif,Iqp,Iqdif)
    VauxqA(i)= VauxqA(i)+Aux1*(sigma*Zlq*Iqp*Iqdif)
    VauxqB(i)= VauxqB(i)+Aux1*(-Zlqp*Iqdif+Zlqdif*Iqp)
    SignCase= "mm"
    CALL Aux_Coef_Swave(SignCase,i,Qxp,Qzp,Qxdif,Qzdif,&
	Aux1,Zlqp,Zlqdif,Iqp,Iqdif)
    VauxqA(i)= VauxqA(i)+Aux1*(sigma*Zlq*Iqp*Iqdif)
    VauxqB(i)= VauxqB(i)+Aux1*(Zlqp*Iqdif-Zlqdif*Iqp)
    IF (Qdm<0.) THEN
     Qzp= -Qdm
     Qxdif= Qx+Qxp
     IF(Qzp<=Qz) THEN
      Qzdif= Qz-Qzp
      SignCase= "pm"
      CALL Aux_Coef_Swave(SignCase,i,Qxp,Qzp,Qxdif,Qzdif,&
	Aux1,Zlqp,Zlqdif,Iqp,Iqdif)
      VauxqA(i)= VauxqA(i)+Aux1*(sigma*Zlq*Iqp*Iqdif)
      VauxqB(i)= VauxqB(i)+Aux1*(Zlqp*Iqdif-Zlqdif*Iqp)
      SignCase= "mp"
      CALL Aux_Coef_Swave(SignCase,i,Qxp,Qzp,Qxdif,Qzdif,&
	Aux1,Zlqp,Zlqdif,Iqp,Iqdif)
      VauxqA(i)= VauxqA(i)+Aux1*(sigma*Zlq*Iqp*Iqdif)
      VauxqB(i)= VauxqB(i)+Aux1*(-Zlqp*Iqdif+Zlqdif*Iqp)
     ELSE
      Qzdif= Qzp-Qz
      SignCase= "pp"
      CALL Aux_Coef_Swave(SignCase,i,Qxp,Qzp,Qxdif,Qzdif,&
	Aux1,Zlqp,Zlqdif,Iqp,Iqdif)
      VauxqA(i)= VauxqA(i)+Aux1*(sigma*Zlq*Iqp*Iqdif)
      VauxqB(i)= VauxqB(i)+Aux1*(Zlqp*Iqdif-Zlqdif*Iqp)
      SignCase= "mm"
      CALL Aux_Coef_Swave(SignCase,i,Qxp,Qzp,Qxdif,Qzdif,&
	Aux1,Zlqp,Zlqdif,Iqp,Iqdif)
      VauxqA(i)= VauxqA(i)+Aux1*(sigma*Zlq*Iqp*Iqdif)
      VauxqB(i)= VauxqB(i)+Aux1*(-Zlqp*Iqdif+Zlqdif*Iqp)
     END IF
    ELSE
     Qzp= Qdm
     Qxdif= Qx-Qxp
     IF(Qzp<=Qz) THEN
      Qzdif= Qz-Qzp
      SignCase= "pm"
      CALL Aux_Coef_Swave(SignCase,i,Qxp,Qzp,Qxdif,Qzdif,&
	Aux1,Zlqp,Zlqdif,Iqp,Iqdif)
      VauxqA(i)= VauxqA(i)+Aux1*(sigma*Zlq*Iqp*Iqdif)
      VauxqB(i)= VauxqB(i)+Aux1*(Zlqp*Iqdif-Zlqdif*Iqp)
      SignCase= "mp"
      CALL Aux_Coef_Swave(SignCase,i,Qxp,Qzp,Qxdif,Qzdif,&
	Aux1,Zlqp,Zlqdif,Iqp,Iqdif)
      VauxqA(i)= VauxqA(i)+Aux1*(sigma*Zlq*Iqp*Iqdif)
      VauxqB(i)= VauxqB(i)+Aux1*(-Zlqp*Iqdif+Zlqdif*Iqp)
     ELSE
      Qzdif= Qzp-Qz
      SignCase= "pp"
      CALL Aux_Coef_Swave(SignCase,i,Qxp,Qzp,Qxdif,Qzdif,&
	Aux1,Zlqp,Zlqdif,Iqp,Iqdif)
      VauxqA(i)= VauxqA(i)+Aux1*(sigma*Zlq*Iqp*Iqdif)
      VauxqB(i)= VauxqB(i)+Aux1*(Zlqp*Iqdif-Zlqdif*Iqp)
      SignCase= "mm"
      CALL Aux_Coef_Swave(SignCase,i,Qxp,Qzp,Qxdif,Qzdif,&
	Aux1,Zlqp,Zlqdif,Iqp,Iqdif)
      VauxqA(i)= VauxqA(i)+Aux1*(sigma*Zlq*Iqp*Iqdif)
      VauxqB(i)= VauxqB(i)+Aux1*(-Zlqp*Iqdif+Zlqdif*Iqp)
     END IF
    END IF
   END DO
   CALL Simpson(VQx,VauxqA,nqx,Aux)
   CoefD= Aux
   CALL Simpson(VQx,VauxqB,nqx,Aux)
   CoefD= CoefD - Aux*Iwave(iqx,kqz)
   CoefD= (AA/3.) * (sigma*Zlq) * Q/(ABS(Qz)) * CoefD

  CASE DEFAULT
   OPEN(98,FILE='Warning_Coef_Swave.wt')
   WRITE(98,*) ' SdecayResCondition= ',SdecayResCondition
   WRITE(98,*) ' SdecayResCondition must be (Exact ) or (Approx) !!'
   CLOSE(98)
   STOP

 END SELECT

ELSE
 CoefD= 0.
END IF

Coef= CoefE+CoefD
RETURN
END SUBROUTINE Coef_Swave

SUBROUTINE Aux_Coef_Swave(SignCase,iqxp,Qxp,Qzp,Qxdif,Qzdif,&
     Aux1,Zlqp,Zlqdif,Iqp,Iqdif)
USE Common_Arrays
USE Common_Params
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL :: Qxp,Qzp
REAL :: ZL,Zlqp,Zlqdif
REAL :: Qxdif,Qzdif,Iqp,Iqdif
REAL :: Qxdifaux,Qzdifaux
REAL :: Aux1
REAL, DIMENSION(nqz) :: Vauxq
REAL, DIMENSION(nqx,nqz) :: ILq,ILqdif
INTEGER :: iqxp
CHARACTER(LEN=2) :: SignCase 

SELECT CASE(SignCase)
 CASE("pp")
  ILqdif= ILp
  ILq= ILp
 CASE("pm")
  ILqdif= ILp
  ILq= ILm
 CASE("mp")
  ILqdif= ILm
  ILq= ILp
 CASE("mm")
  ILqdif= ILm
  ILq= ILm
 CASE DEFAULT
  OPEN(98,FILE='Warning_Aux_Coef_Swave.wt')
  WRITE(98,*) ' SignCase= ',SignCase
  WRITE(98,*) ' SignCase must be pp, pm, mp, mm !! '
  CLOSE(98)
  STOP
END SELECT
Zlqp= ZL(Qxp,Qzp)
Zlqdif= ZL(Qxdif,Qzdif)
IF(ABS(Qzp)>VQz(nqz) .AND. ABS(Qzdif)>VQz(nqz)) THEN
 Aux1= 0.
 Iqp= 0.
 Iqdif= 0.
ELSE
 Aux1= (Qxp**2*Qxdif**2+Qzp**2*Qzdif**2)/(Qxdif**2+Qzdif**2+EpsMin) &
      /(Qxp**2+Qzp**2)
 IF (Qxdif>=0.) THEN
  Qxdifaux= Qxdif
 ELSE
  Qxdifaux= -Qxdif
 END IF
 IF (Qzdif>=0.) THEN
  Qzdifaux= Qzdif
 ELSE
  Qzdifaux= -Qzdif
 END IF
 CALL Aitp2d(nqx,nqz,VQx,VQz,ILqdif,Qxdifaux,Qzdifaux,Iqdif)
 Vauxq(:)= ILq(iqxp,:)
 CALL Aitp1d(nqz,VQz,Vauxq,Qzp,Iqp)
END IF
RETURN
END SUBROUTINE Aux_Coef_Swave

REAL FUNCTION ZL(Qx,Qz)
IMPLICIT NONE
REAL, INTENT(in) :: Qx,Qz
REAL :: Q2

Q2= Qx**2+Qz**2
ZL= SQRT(1.+1.5*Q2)
IF (Qz < 0.) THEN
 ZL= -ZL
ELSE
END IF
RETURN
END FUNCTION ZL

REAL FUNCTION ZS(Qx,Qz)
IMPLICIT NONE
REAL, INTENT(in) :: Qx,Qz
REAL :: Q,Q2
REAL :: RMiMe,AA,Geff

CALL Definitions(RMiMe,AA,Geff)
Q2= Qx**2+Qz**2
Q= SQRT(Q2)
ZS= Q*AA/SQRT(1.+Q2/2.)
IF (Qz < 0.) THEN
 ZS= -ZS
ELSE
END IF
RETURN
END FUNCTION ZS

SUBROUTINE Fnorm(Anorm)
USE Common_Arrays
USE Common_Params
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
INTEGER :: l,m
REAL :: Anorm
REAL, DIMENSION(nux) :: Vintux
REAL, DIMENSION(nuz) :: Vintuz

DO l= 1,nux
 DO m= 1,nuz
  Vintuz(m)= Fe(l,m)
 END DO
 CALL Simpson(VUz,Vintuz,nuz,Anorm)
 Vintux(l)= Anorm
END DO
CALL Simpson(VUx,Vintux,nux,Anorm)
Anorm= 2.*Anorm  ! Symmetry in Ux.
! Anorm commented out:
! DO l= 1,nux
! DO m= 1,nuz
!  Fe(l,m)= Fe(l,m)/Anorm
! END DO
! END DO
RETURN
END SUBROUTINE Fnorm

SUBROUTINE Energy(Epart,Ewave,EppEw)
USE Common_Arrays
USE Common_Params
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
INTEGER :: i,k
REAL :: Ewave,Epart,EppEw
REAL :: Ux,Uz,U2,Aux
REAL :: Qx,Qz,Q,Muq
REAL :: RMiMe,AA,Geff
REAL, DIMENSION(nux) :: Vintux
REAL, DIMENSION(nuz) :: Vintuz
REAL, DIMENSION(nqx) :: Vintqx
REAL, DIMENSION(nqz) :: Vintqz

CALL Definitions(RMiMe,AA,Geff)

DO k= 1,nuz
 Uz= VUz(k)
 DO i= 1,nux
  Ux= Vux(i)
  U2= Ux**2+Uz**2
  Vintux(i)= U2*Fe(i,k)
 END DO
 CALL Simpson(VUx,Vintux,nux,Aux)
 Vintuz(k)= Aux
END DO
CALL Simpson(VUz,Vintuz,nuz,Aux)
Epart= Aux/2.

DO i= 1,nqx
 Qx= VQx(i)
 DO k= 1,nqz
  Qz= VQz(k) 
  Q= SQRT(Qx**2+Qz**2)
  Muq= Q**3*AA/2.
  Vintqz(k)= ILp(i,k)+ILm(i,k)+Muq*(ISp(i,k)+ISm(i,k))
 END DO
 CALL Simpson(VQz,Vintqz,nqz,Aux)
 Vintqx(i)= Aux
END DO
CALL Simpson(VQx,Vintqx,nqx,Aux)
Ewave= Aux

EppEw= Epart+Ewave

RETURN
END SUBROUTINE Energy

SUBROUTINE Output(WriteChoice)
USE Common_Arrays
USE Common_Params
USE Math_Constants
IMPLICIT NONE
INTEGER :: i,k
INTEGER, PARAMETER :: Step=2
REAL :: Aux
REAL :: Qx,Qz,Q,Muq
REAL :: RMiMe,AA,Geff
REAL, DIMENSION(nqx) :: Vintqx
REAL, DIMENSION(nqz) :: Iqzm,Iqzp,Iqzm2,Iqzp2
REAL, DIMENSION(nux) :: Vintux
REAL, DIMENSION(nuz) :: Feuz
CHARACTER(LEN=3) :: WriteChoice

CALL Definitions(RMiMe,AA,Geff)

SELECT CASE(WriteChoice)

 CASE("Fe0")
  OPEN(1,FILE='Fe0.wt')
  DO i= 1,nux,Step
   WRITE(1,*)' '
   DO k= 1,nuz,Step
    IF(ABS(Fe0(nux+1-i,k))<=EpsMin) THEN
     Aux= 0.E0
    ELSE
     Aux= Fe0(nux+1-i,k)
    END IF
    WRITE(1,*) -VUx(nux+1-i),VUz(k),Aux
   END DO
  END DO
  DO i= 1,nux,Step
   WRITE(1,*)' '
   DO k= 1,nuz,Step
    IF(ABS(Fe0(i,k))<=EpsMin) THEN
     Aux= 0.E0
    ELSE
     Aux= Fe0(i,k)
    END IF
    WRITE(1,*) VUx(i),VUz(k),Aux
   END DO
  END DO
  CLOSE(1)

 CASE("Fe ")
  OPEN(1,FILE='Fe.wt')
  DO i= 1,nux,Step
   WRITE(1,*)' '
   DO k= 1,nuz,Step
    IF(ABS(Fe(nux+1-i,k))<=EpsMin) THEN
     Aux= 0.E0
    ELSE
     Aux= Fe(nux+1-i,k)
    END IF
    WRITE(1,*) -VUx(nux+1-i),VUz(k),Aux
   END DO
  END DO
  DO i= 1,nux,Step
   WRITE(1,*)' '
   DO k= 1,nuz,Step
    IF(ABS(Fe(i,k))<=EpsMin) THEN
     Aux= 0.E0
    ELSE
     Aux= Fe(i,k)
    END IF
    WRITE(1,*) VUx(i),VUz(k),Aux
   END DO
  END DO
  CLOSE(1)

 CASE("Ai ")
  OPEN(1,FILE='Ai.wt')
  DO i= 1,nux,Step
   WRITE(1,*)' '
   DO k= 1,nuz,Step
    WRITE(1,*) -VUx(nux+1-i),VUz(k),Ax(nux+1-i,k),Az(nux+1-i,k)
   END DO
  END DO
  DO i= 1,nux,Step
   WRITE(1,*)' '
   DO k= 1,nuz,Step
    WRITE(1,*) VUx(i),VUz(k),Ax(i,k),Az(i,k)
   END DO
  END DO
  CLOSE(1)

 CASE("Dij")
  OPEN(1,FILE='Dij.wt')
  DO i= 1,nux,Step
   WRITE(1,*)' '
   DO k= 1,nuz,Step
    WRITE(1,*) -VUx(nux+1-i),VUz(k),Dxx(nux+1-i,k),Dxz(nux+1-i,k),&
	Dzz(nux+1-i,k)
   END DO
  END DO
  DO i= 1,nux,Step
   WRITE(1,*)' '
   DO k= 1,nuz,Step
    WRITE(1,*) VUx(i),VUz(k),Dxx(i,k),Dxz(i,k),Dzz(i,k)
   END DO
  END DO
  CLOSE(1)

 CASE("IL ")
  OPEN(1,FILE='IL.wt')
  DO i= 1,nqx,Step
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    WRITE(1,*) -VQx(nqx+1-i),-VQz(nqz+1-k),ILm(nqx+1-i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    WRITE(1,*) -VQx(nqx+1-i),VQz(k),ILp(nqx+1-i,k)
   END DO
  END DO
  DO i= 1,nqx,Step
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    WRITE(1,*) VQx(i),-VQz(nqz+1-k),ILm(i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    WRITE(1,*) VQx(i),VQz(k),ILp(i,k)
   END DO
  END DO
  CLOSE(1)

 CASE("IS ")
  OPEN(1,FILE='IS.wt')
  DO i= 1,nqx,Step
   Qx= -VQx(nqx+1-i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    Q= SQRT(Qx**2+Qz**2)
    Muq= Q**3*AA/2.
    WRITE(1,*) Qx,-Qz,ISm(nqx+1-i,nqz+1-k),Muq*ISm(nqx+1-i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    Q= SQRT(Qx**2+Qz**2)
    Muq= Q**3*AA/2.
    WRITE(1,*) Qx,Qz,ISp(nqx+1-i,k),Muq*ISp(nqx+1-i,k)
   END DO
  END DO
  DO i= 1,nqx,Step
   Qx= VQx(i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    Q= SQRT(Qx**2+Qz**2)
    Muq= Q**3*AA/2.
    WRITE(1,*) Qx,-Qz,ISm(i,nqz+1-k),Muq*ISm(i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    Q= SQRT(Qx**2+Qz**2)
    Muq= Q**3*AA/2.
    WRITE(1,*) Qx,Qz,ISp(i,k),Muq*ISp(i,k)
   END DO
  END DO
  CLOSE(1)

 CASE("IL0")
  OPEN(1,FILE='IL0.wt')
  DO i= 1,nqx,Step
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    WRITE(1,*) -VQx(nqx+1-i),-VQz(nqz+1-k),ILm(nqx+1-i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    WRITE(1,*) -VQx(nqx+1-i),VQz(k),ILp(nqx+1-i,k)
   END DO
  END DO
  DO i= 1,nqx,Step
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    WRITE(1,*) VQx(i),-VQz(nqz+1-k),ILm(i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    WRITE(1,*) VQx(i),VQz(k),ILp(i,k)
   END DO
  END DO
  CLOSE(1)
  OPEN(1,FILE='IL01D.wt')
  DO k= 1,nqz
   DO i= 1,nqx
    Vintqx(i)= ILm(i,k)
   END DO
   CALL Simpson(VQx,Vintqx,nqx,Aux)
   Iqzm(k)= Aux
   DO i= 1,nqx
    Vintqx(i)= ILp(i,k)
   END DO
   CALL Simpson(VQx,Vintqx,nqx,Aux)
   Iqzp(k)= Aux
  END DO
  DO k= 1,nqz,Step
   WRITE(1,*) -VQz(nqz+1-k),Iqzm(nqz+1-k)
  END DO
  DO k= 1,nqz,Step
   WRITE(1,*) VQz(k),Iqzp(k)
  END DO
  CLOSE(1)

 CASE("IS0")
  OPEN(1,FILE='IS0.wt')
  DO i= 1,nqx,Step
   Qx= -VQx(nqx+1-i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    Q= SQRT(Qx**2+Qz**2)
    Muq= Q**3*AA/2.
    WRITE(1,*) Qx,-Qz,ISm(nqx+1-i,nqz+1-k),Muq*ISm(nqx+1-i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    Q= SQRT(Qx**2+Qz**2)
    Muq= Q**3*AA/2.
    WRITE(1,*) Qx,Qz,ISp(nqx+1-i,k),Muq*ISp(nqx+1-i,k)
   END DO
  END DO
  DO i= 1,nqx,Step
   Qx= VQx(i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    Q= SQRT(Qx**2+Qz**2)
    Muq= Q**3*AA/2.
    WRITE(1,*) Qx,-Qz,ISm(i,nqz+1-k),Muq*ISm(i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    Q= SQRT(Qx**2+Qz**2)
    Muq= Q**3*AA/2.
    WRITE(1,*) Qx,Qz,ISp(i,k),Muq*ISp(i,k)
   END DO
  END DO
  CLOSE(1)
  OPEN(1,FILE='IS01D.wt')
  DO k= 1,nqz
   Qz= VQz(k)
   DO i= 1,nqx
    Vintqx(i)= ISm(i,k)
   END DO
   CALL Simpson(VQx,Vintqx,nqx,Aux)
   Iqzm(k)= Aux
   DO i= 1,nqx
    Qx= VQx(i)
    Q= SQRT(Qx**2+Qz**2)
    Muq= Q**3*AA/2.
    Vintqx(i)= Muq*ISm(i,k)
   END DO
   CALL Simpson(VQx,Vintqx,nqx,Aux)
   Iqzm2(k)= Aux
   DO i= 1,nqx
    Vintqx(i)= ISp(i,k)
   END DO
   CALL Simpson(VQx,Vintqx,nqx,Aux)
   Iqzp(k)= Aux
   DO i= 1,nqx
    Qx= VQx(i)
    Q= SQRT(Qx**2+Qz**2)
    Muq= Q**3*AA/2.
    Vintqx(i)= Muq*ISp(i,k)
   END DO
   CALL Simpson(VQx,Vintqx,nqx,Aux)
   Iqzp2(k)= Aux
  END DO
  DO k= 1,nqz,Step
   WRITE(1,*) -VQz(nqz+1-k),Iqzm(nqz+1-k),Iqzm2(nqz+1-k)
  END DO
  DO k= 1,nqz,Step
   WRITE(1,*) VQz(k),Iqzp(k),Iqzp2(k)
  END DO
  CLOSE(1)

 CASE("Fe1")
  OPEN(1,FILE='Fe1D.wt')
  DO k= 1,nuz
   DO i= 1,nux
    Vintux(i)= Fe(i,k)
   END DO
   CALL Simpson(VUx,Vintux,nux,Aux)
   Feuz(k)= Aux
  END DO
  DO k= 1,nuz,Step
   WRITE(1,*) VUz(k),Feuz(k)
  END DO
  CLOSE(1)

 CASE("IL1")
  OPEN(1,FILE='IL1D.wt')
  DO k= 1,nqz
   DO i= 1,nqx
    Vintqx(i)= ILm(i,k)
   END DO
   CALL Simpson(VQx,Vintqx,nqx,Aux)
   Iqzm(k)= Aux
   DO i= 1,nqx
    Vintqx(i)= ILp(i,k)
   END DO
   CALL Simpson(VQx,Vintqx,nqx,Aux)
   Iqzp(k)= Aux
  END DO
  DO k= 1,nqz,Step
   WRITE(1,*) -VQz(nqz+1-k),Iqzm(nqz+1-k)
  END DO
  DO k= 1,nqz,Step
   WRITE(1,*) VQz(k),Iqzp(k)
  END DO
  CLOSE(1)

 CASE("IS1")
  OPEN(1,FILE='IS1D.wt')
  DO k= 1,nqz
   Qz= VQz(k)
   DO i= 1,nqx
    Vintqx(i)= ISm(i,k)
   END DO
   CALL Simpson(VQx,Vintqx,nqx,Aux)
   Iqzm(k)= Aux
   DO i= 1,nqx
    Qx= VQx(i)
    Q= SQRT(Qx**2+Qz**2)
    Muq= Q**3*AA/2.
    Vintqx(i)= Muq*ISm(i,k)
   END DO
   CALL Simpson(VQx,Vintqx,nqx,Aux)
   Iqzm2(k)= Aux
   DO i= 1,nqx
    Vintqx(i)= ISp(i,k)
   END DO
   CALL Simpson(VQx,Vintqx,nqx,Aux)
   Iqzp(k)= Aux
   DO i= 1,nqx
    Qx= VQx(i)
    Q= SQRT(Qx**2+Qz**2)
    Muq= Q**3*AA/2.
    Vintqx(i)= Muq*ISp(i,k)
   END DO
   CALL Simpson(VQx,Vintqx,nqx,Aux)
   Iqzp2(k)= Aux
  END DO
  DO k= 1,nqz,Step
   WRITE(1,*) -VQz(nqz+1-k),Iqzm(nqz+1-k),Iqzm2(nqz+1-k)
  END DO
  DO k= 1,nqz,Step
   WRITE(1,*) VQz(k),Iqzp(k),Iqzp2(k)
  END DO
  CLOSE(1)

  CASE DEFAULT

END SELECT
RETURN

END SUBROUTINE Output


!---------------------------------------------------------------------
! Mathematical Routines:
!---------------------------------------------------------------------

SUBROUTINE Adi5p(Ucrit,Dux,Duz,DTau,BoundaryCondition)

! See Notes J-06 (also G-05).
! Fe_old(nux,nuz): Keeps the values of the function at the start of each
! 		 sub-interval.

USE Common_Arrays
USE Common_Params
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL :: Ux,Uz
REAL :: Auxx,Auxz,Auxxx,Auxzz,Auxxz
REAL :: DTau,Dux,Duz,U2,Ucrit
REAL, DIMENSION(nux) :: Alpha1,Beta1,Gamma1,Psi1
REAL, DIMENSION(nuz) :: Alpha2,Beta2,Gamma2,Psi2
REAL, DIMENSION(nux) :: Faux1
REAL, DIMENSION(nuz) :: Faux2
REAL, DIMENSION(nux) :: Baux1,Gaux1
REAL, DIMENSION(nuz) :: Baux2,Gaux2
REAL :: Aux
INTEGER :: i,k
INTEGER :: ip1,im1,kp1,km1,ip2,im2,kp2,km2
CHARACTER(LEN=10) :: BoundaryCondition

! First sub-interval:

CALL Coef_D

Fe_old= Fe
Auxx= DTau/4./Dux
Auxz= DTau/4./Duz
Auxxx= DTau/4./Dux/Dux
Auxxz= DTau/8./Dux/Duz
Auxzz= DTau/4./Duz/Duz

DO k= 2,nuz-1
 kp1= k+1
 km1= k-1
 kp2= k+2
 km2= k-2
 DO i= 2,nux-1
  ip1= i+1
  im1= i-1
  Alpha1(i)= Ax(im1,k)*Auxx - (Dxx(im1,k)+Dxx(i,k))*Auxxx
  Beta1(i)= 1. + (Dxx(im1,k)+2.*Dxx(i,k)+Dxx(ip1,k))*Auxxx
  Gamma1(i)= - Ax(ip1,k)*Auxx - (Dxx(i,k)+Dxx(ip1,k))*Auxxx
  IF (k==2 .OR. k==nuz-1) THEN
   Psi1(i)= Fe(i,k) + (Az(i,kp1)*Fe(i,kp1)-Az(i,km1)*Fe(i,km1))*Auxz &
	+( Dxz(ip1,k)*(Fe(ip1,kp1)-Fe(ip1,km1)) & 
	- Dxz(im1,k)*(Fe(im1,kp1)-Fe(im1,km1)) )*Auxxz &
	+( Dzx(i,kp1)*(Fe(ip1,kp1)-Fe(im1,kp1)) &
	- Dzx(i,km1)*(Fe(ip1,km1)-Fe(im1,km1)) )*Auxxz &
	+( (Dzz(i,k)+Dzz(i,kp1))*(Fe(i,kp1)-Fe(i,k)) &
	- (Dzz(i,km1)+Dzz(i,k))*(Fe(i,k)-Fe(i,km1)) )*Auxzz
  ELSE
   Psi1(i)= Fe(i,k) + (2.*(Az(i,kp1)*Fe(i,kp1)-Az(i,km1)*Fe(i,km1)) &
	- (Az(i,kp2)*Fe(i,kp2)-Az(i,km2)*Fe(i,km2))/4. ) * (2.*Auxz/3.) &
	+ (2.*(Dzx(i,kp1)*Fe(ip1,kp1)-Dzx(i,km1)*Fe(ip1,km1)) &
	- (Dzx(i,kp2)*Fe(ip1,kp2)-Dzx(i,km2)*Fe(ip1,km2))/4. ) * (2.*Auxxz/3.) &
	- (2.*(Dzx(i,kp1)*Fe(im1,kp1)-Dzx(i,km1)*Fe(im1,km1)) &
	- (Dzx(i,kp2)*Fe(im1,kp2)-Dzx(i,km2)*Fe(im1,km2))/4. ) * (2.*Auxxz/3.) &
	+ Dxz(ip1,k)*(2.*(Fe(ip1,kp1)-Fe(ip1,km1)) &
	- (Fe(ip1,kp2)-Fe(ip1,km2))/4.) * (2.*Auxxz/3.) &
	+ Dxz(im1,k)*(2.*(Fe(im1,kp1)-Fe(im1,km1)) &
	- (Fe(im1,kp2)-Fe(im1,km2))/4.) * (2.*Auxxz/3.) &
	+ ((Dzz(i,kp1)+Dzz(i,k))*Fe(i,kp1)-(Dzz(i,km1)+Dzz(i,k))*Fe(i,k) &
	- (Dzz(i,kp1)*(Fe(i,kp1)+Fe(i,kp2))-Dzz(i,km1)*(Fe(i,k) &
	+ Fe(i,km1)))/8.) * (16./9.)*Auxzz &
	- ((Dzz(i,kp1)+Dzz(i,k))*Fe(i,k)-(Dzz(i,km1)+Dzz(i,k))*Fe(i,km1) &
	- (Dzz(i,kp1)*(Fe(i,k)+Fe(i,kp1))-Dzz(i,km1)*(Fe(i,k-1) &
	+ Fe(i,km2)))/8.) * (16./9.)*Auxzz &
	- ((Dzz(i,kp1)+Dzz(i,k))*(Fe(i,kp1)+Fe(i,kp2))-(Dzz(i,km1) &
	+ Dzz(i,k))*(Fe(i,k)+Fe(i,kp1))-(Dzz(i,kp1)*Fe(i,kp2) &
	- Dzz(i,km1)*Fe(i,k))/2.) * (Auxzz/9.) &
	+ ((Dzz(i,kp1)+Dzz(i,k))*(Fe(i,k)+Fe(i,km1))-(Dzz(i,km1) &
	+ Dzz(i,k))*(Fe(i,km1)+Fe(i,km2))-(Dzz(i,kp1)*Fe(i,k) &
	- Dzz(i,km1)*Fe(i,km2))/2.) * (Auxzz/9.)
  END IF
 END DO

 CALL Tridag(2,nux-1,Alpha1,Beta1,Gamma1,Psi1,Faux1,Baux1,Gaux1)

 SELECT CASE(BoundaryCondition)
  CASE("ConstantF ")
   ! Boundary conditions: (dF/dx)=0 at Ux=0, F=constant at Ux=Ulim.
   Faux1(1)= Faux1(2)    ! approximated derivative (forward)
   Faux1(nux)= Fe(nux,k)
  CASE("ConstantDF")
   ! Boundary conditions: (dF/dx)=0 at Ux=0, (dF/dx)=constant at Ux=Ulim.
   Faux1(1)= Faux1(2)    ! approximated derivative (forward)
   Aux= Faux1(nux-1)+Fepxp(k)*Dux
   IF (Aux<0.) THEN
    Faux1(nux)= 0.
   ELSE
    Faux1(nux)= Aux
   END IF
  CASE DEFAULT
   OPEN(98,FILE='Warning_Adi5p.wt')
   WRITE(98,*) ' BoundaryCondition= ',BoundaryCondition
   WRITE(98,*) ' Boundary Condition must be (ConstantF ) or (ConstantDF) !!'
   CLOSE(98)
   STOP
 END SELECT

 ! Correcting instabilities:
 CALL Cor_Ampli(Faux1,nux)

 DO i= 1,nux
  Fe(i,k)= Faux1(i)
 END DO
END DO

SELECT CASE(BoundaryCondition)
 CASE("ConstantF ")
  ! The boundary conditions in 'Uz' are automatically taken into account.
 CASE("ConstantDF")
  ! Boundary conditions in 'Uz', (dF/dz)= constant, at Uz=-Ulim and Uz=Ulim.
  DO i= 1,nux
   Aux= Fe(i,2)-Fepzm(i)*Duz
   IF (Aux<0.) THEN
    Fe(i,1)= 0.
   ELSE
    Fe(i,1)= Aux
   END IF
   Aux= Fe(i,nuz-1)+Fepzp(i)*Duz
   IF (Aux<0.) THEN
    Fe(i,nuz)= 0.
   ELSE
    Fe(i,nuz)= Aux
   END IF
  END DO
 CASE DEFAULT
   OPEN(98,FILE='Warning_Adi5p.wt')
   WRITE(98,*) ' BoundaryCondition= ',BoundaryCondition
   WRITE(98,*) ' Boundary Condition must be (ConstantF ) or (ConstantDF) !!'
   CLOSE(98)
   STOP
END SELECT

DO k= 1,nuz
 Uz= VUz(k)
 DO i= 1,nux
  Ux= VUx(i)
  U2= Ux**2+Uz**2
  IF (U2 >= Ucrit**2) THEN
   Fe_old(i,k)= Fe(i,k)
  ELSE
   Fe(i,k)= Fe_old(i,k)
  END IF
 END DO
END DO

! Second sub-interval:

CALL Coef_D

Fe_old= Fe
Auxx= DTau/4./Dux
Auxz= DTau/4./Duz
Auxxx= DTau/4./Dux/Dux
Auxxz= DTau/8./Dux/Duz
Auxzz= DTau/4./Duz/Duz

DO i= 2,nux-1
 ip1= i+1
 im1= i-1
 ip2= i+2
 im2= i-2
 DO k= 2,nuz-1
  kp1= k+1
  km1= k-1
  Alpha2(k)= Az(i,km1)*Auxz - (Dzz(i,km1)+Dzz(i,k))*Auxzz
  Beta2(k)= 1 + (Dzz(i,km1)+2.*Dzz(i,k)+Dzz(i,kp1))*Auxzz
  Gamma2(k)= - Az(i,kp1)*Auxz - (Dzz(i,k)+Dzz(i,kp1))*Auxzz
  IF (i==2 .OR. i==nux-1) THEN
   Psi2(k)= Fe(i,k) + (Ax(ip1,k)*Fe(ip1,k)-Ax(im1,k)*Fe(im1,k))*Auxx &
	+( Dxz(ip1,k)*(Fe(ip1,kp1)-Fe(ip1,km1)) & 
	- Dxz(im1,k)*(Fe(im1,kp1)-Fe(im1,km1)) )*Auxxz &
	+( Dzx(i,kp1)*(Fe(ip1,kp1)-Fe(im1,kp1)) &
	- Dzx(i,km1)*(Fe(ip1,km1)-Fe(im1,km1)) )*Auxxz &
	+( (Dxx(i,k)+Dxx(ip1,k))*(Fe(ip1,k)-Fe(i,k)) &
	- (Dxx(im1,k)+Dxx(i,k))*(Fe(i,k)-Fe(im1,k)) )*Auxxx
  ELSE
   Psi2(k)= Fe(i,k) + ((2./3.)*(Ax(ip1,k)*Fe(ip1,k)-Ax(im1,k)*Fe(im1,k)) &
	- (Ax(ip2,k)*Fe(ip2,k)-Ax(im2,k)*Fe(im2,k))/12. ) * 2.*Auxx &
	+ (2.*(Dxz(ip1,k)*Fe(ip1,kp1)-Dxz(im1,k)*Fe(im1,kp1)) &
	- (Dxz(ip2,k)*Fe(ip2,kp1)-Dxz(im2,k)*Fe(im2,kp1))/4. ) * (2.*Auxxz/3.) &
	- (2.*(Dxz(ip1,k)*Fe(ip1,km1)-Dxz(im1,k)*Fe(im1,km1)) &
	- (Dxz(ip2,k)*Fe(ip2,km1)-Dxz(im2,k)*Fe(im2,km1))/4. ) * (2.*Auxxz/3.) &
	+ Dzx(i,kp1)*(2.*(Fe(ip1,kp1)-Fe(im1,kp1)) &
	- (Fe(ip2,kp1)-Fe(im2,kp1))/4.) * (2.*Auxxz/3.) &
	+ Dxz(i,km1)*(2.*(Fe(ip1,km1)-Fe(im1,km1)) &
	- (Fe(ip2,km1)-Fe(im2,km1))/4.) * (2.*Auxxz/3.) &
	+ ((Dxx(ip1,k)+Dxx(i,k))*Fe(ip1,k)-(Dxx(im1,k)+Dxx(i,k))*Fe(i,k) &
	- (Dxx(ip1,k)*(Fe(ip1,k)+Fe(ip2,k))-Dxx(im1,k)*(Fe(i,k) &
	+ Fe(im1,k)))/8.) * (16./9.)*Auxxx &
	- ((Dxx(ip1,k)+Dxx(i,k))*Fe(i,k)-(Dxx(im1,k)+Dxx(i,k))*Fe(im1,k) &
	- (Dxx(ip1,k)*(Fe(i,k)+Fe(ip1,k))-Dxx(im1,k)*(Fe(im1,k) &
	+ Fe(im2,k)))/8.) * (16./9.)*Auxxx &
	- ((Dxx(ip1,k)+Dxx(i,k))*(Fe(ip1,k)+Fe(ip2,k))-(Dxx(im1,k) &
	+ Dxx(i,k))*(Fe(i,k)+Fe(ip1,k))-(Dxx(ip1,k)*Fe(ip2,k) &
	- Dxx(im1,k)*Fe(i,k))/2.) * (Auxxx/9.) &
	+ ((Dxx(ip1,k)+Dxx(i,k))*(Fe(i,k)+Fe(im1,k))-(Dzz(im1,k) &
	+ Dzz(i,k))*(Fe(im1,k)+Fe(im2,k))-(Dxx(ip1,k)*Fe(i,k) &
	- Dzz(im1,k)*Fe(im2,k))/2.) * (Auxxx/9.)
  END IF
 END DO

 CALL Tridag(2,nuz-1,Alpha2,Beta2,Gamma2,Psi2,Faux2,Baux2,Gaux2)

 SELECT CASE(BoundaryCondition)
  CASE("ConstantF ")
   ! Boundary conditions: F=constant at Uz=-Ulim and Uz=Ulim.
   Faux2(1)= Fe(i,1)
   Faux2(nuz)= Fe(i,nuz)
  CASE("ConstantDF")
   ! Boundary conditions in 'Uz', (dF/dz)= constant, at Uz=-Ulim and Uz=Ulim.
   Aux= Faux2(2)-Fepzm(i)*Duz
   IF (Aux<0.) THEN
    Faux2(1)= 0.
   ELSE
    Faux2(1)= Aux
   END IF
   Aux= Faux2(nuz-1)+Fepzp(i)*Duz
   IF (Aux<0.) THEN
    Faux2(nuz)= 0.
   ELSE
    Faux2(nuz)= Aux
   END IF
  CASE DEFAULT
   OPEN(98,FILE='Warning_Adi5p.wt')
   WRITE(98,*) ' BoundaryCondition= ',BoundaryCondition
   WRITE(98,*) ' Boundary Condition must be (ConstantF ) or (ConstantDF) !!'
   CLOSE(98)
   STOP
 END SELECT

 ! Correcting instabilities:
 CALL Cor_Ampli(Faux2,nuz)

 DO k= 1,nuz
  Fe(i,k)= Faux2(k)
 END DO
END DO

SELECT CASE(BoundaryCondition)
 CASE("ConstantF ")
  ! The boundary conditions at 'Uz=Ulim' are automatically taken into account.
  ! For 'Ux=0' we assume zero derivative (approximated forward derivative):
  DO k= 1,nuz
   Fe(1,k)= Fe(2,k)
  END DO
 CASE("ConstantDF")
  ! Boundary conditions: (dF/dx)=0 at Ux=0, (dF/dx)=constant at Ux=Ulim.
  DO k= 1,nuz
   Fe(1,k)= Fe(2,k)
   Aux= Fe(nux-1,k)+Fepxp(k)*Dux
   IF (Aux < 0.) THEN
    Fe(nux,k)= 0.
   ELSE
    Fe(nux,k)= Aux 
   END IF
  END DO
 CASE DEFAULT
  OPEN(98,FILE='Warning_Adi5p.wt')
  WRITE(98,*) ' BoundaryCondition= ',BoundaryCondition
  WRITE(98,*) ' Boundary Condition must be (ConstantF ) or (ConstantDF) !!'
  CLOSE(98)
  STOP
END SELECT


DO k= 1,nuz
 Uz= VUz(k)
 DO i= 1,nux
  Ux= VUx(i)
  U2= Ux**2+Uz**2
  IF (U2 >= Ucrit**2) THEN
   Fe_old(i,k)= Fe(i,k)
  ELSE
   Fe(i,k)= Fe_old(i,k)
  END IF
 END DO
END DO

RETURN
END SUBROUTINE Adi5p

SUBROUTINE Adi3p(Ucrit,Dux,Duz,DTau,BoundaryCondition)

! See Notes G-05.
! Fe_old(nux,nuz): Keeps the values of the function at the start of each
! 		 sub-interval.

USE Common_Arrays
USE Common_Params
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL :: Ux,Uz
REAL :: Auxx,Auxz,Auxxx,Auxzz,Auxxz
REAL :: DTau,Dux,Duz,U2,Ucrit
REAL, DIMENSION(nux) :: Alpha1,Beta1,Gamma1,Psi1
REAL, DIMENSION(nuz) :: Alpha2,Beta2,Gamma2,Psi2
REAL, DIMENSION(nux) :: Faux1
REAL, DIMENSION(nuz) :: Faux2
REAL, DIMENSION(nux) :: Baux1,Gaux1
REAL, DIMENSION(nuz) :: Baux2,Gaux2
REAL :: Aux
INTEGER :: i,k
INTEGER :: ip1,im1,kp1,km1
CHARACTER(LEN=10) :: BoundaryCondition

! First sub-interval:

CALL Coef_D

Fe_old= Fe
Auxx= DTau/4./Dux
Auxz= DTau/4./Duz
Auxxx= DTau/4./Dux/Dux
Auxxz= DTau/8./Dux/Duz
Auxzz= DTau/4./Duz/Duz

DO k= 2,nuz-1
 kp1= k+1
 km1= k-1
 DO i= 2,nux-1
  ip1= i+1
  im1= i-1
  Alpha1(i)= Ax(im1,k)*Auxx - (Dxx(im1,k)+Dxx(i,k))*Auxxx
  Beta1(i)= 1. + (Dxx(im1,k)+2.*Dxx(i,k)+Dxx(ip1,k))*Auxxx
  Gamma1(i)= - Ax(ip1,k)*Auxx - (Dxx(i,k)+Dxx(ip1,k))*Auxxx
  Psi1(i)= Fe(i,k) + (Az(i,kp1)*Fe(i,kp1)-Az(i,km1)*Fe(i,km1))*Auxz &
	+( Dxz(ip1,k)*(Fe(ip1,kp1)-Fe(ip1,km1)) & 
	- Dxz(im1,k)*(Fe(im1,kp1)-Fe(im1,km1)) )*Auxxz &
	+( Dzx(i,kp1)*(Fe(ip1,kp1)-Fe(im1,kp1)) &
	- Dzx(i,km1)*(Fe(ip1,km1)-Fe(im1,km1)) )*Auxxz &
	+( (Dzz(i,k)+Dzz(i,kp1))*(Fe(i,kp1)-Fe(i,k)) &
	- (Dzz(i,km1)+Dzz(i,k))*(Fe(i,k)-Fe(i,km1)) )*Auxzz
 END DO

 CALL Tridag(2,nux-1,Alpha1,Beta1,Gamma1,Psi1,Faux1,Baux1,Gaux1)

 SELECT CASE(BoundaryCondition)
  CASE("ConstantF ")
   ! Boundary conditions: (dF/dx)=0 at Ux=0, F=constant at Ux=Ulim.
   Faux1(1)= Faux1(2)    ! approximated derivative (forward)
   Faux1(nux)= Fe(nux,k)
  CASE("ConstantDF")
   ! Boundary conditions: (dF/dx)=0 at Ux=0, (dF/dx)=constant at Ux=Ulim.
   Faux1(1)= Faux1(2)    ! approximated derivative (forward)
   Aux= Faux1(nux-1)+Fepxp(k)*Dux
   IF (Aux < 0.) THEN
    Faux1(nux)= 0.
   ELSE
    Faux1(nux)= Aux
   END IF
  CASE DEFAULT
   OPEN(98,FILE='Warning_Adi3p.wt')
   WRITE(98,*) ' BoundaryCondition= ',BoundaryCondition
   WRITE(98,*) ' Boundary Condition must be (ConstantF ) or (ConstantDF) !!'
   CLOSE(98)
   STOP
 END SELECT

 ! Correcting instabilities:
 CALL Cor_Ampli(Faux1,nux)

 DO i= 1,nux
  Fe(i,k)= Faux1(i)
 END DO
END DO

SELECT CASE(BoundaryCondition)
 CASE("ConstantF ")
  ! The boundary conditions in 'Uz' are automatically taken into account.
 CASE("ConstantDF")
  ! Boundary conditions in 'Uz', (dF/dz)= constant, at Uz=-Ulim and Uz=Ulim.
  DO i= 1,nux
   Aux= Fe(i,2)-Fepzm(i)*Duz
   IF (Aux<0.) THEN
    Fe(i,1)= 0.
   ELSE
    Fe(i,1)= Aux
   END IF
   Aux= Fe(i,nuz-1)+Fepzp(i)*Duz
   IF (Aux<0.) THEN
    Fe(i,nuz)= 0.
   ELSE
    Fe(i,nuz)= Aux
   END IF
  END DO
 CASE DEFAULT
   OPEN(98,FILE='Warning_Adi5p.wt')
   WRITE(98,*) ' BoundaryCondition= ',BoundaryCondition
   WRITE(98,*) ' Boundary Condition must be (ConstantF ) or (ConstantDF) !!'
   CLOSE(98)
   STOP
END SELECT

DO k= 1,nuz
 Uz= VUz(k)
 DO i= 1,nux
  Ux= VUx(i)
  U2= Ux**2+Uz**2
  IF (U2 >= Ucrit**2) THEN
   Fe_old(i,k)= Fe(i,k)
  ELSE
   Fe(i,k)= Fe_old(i,k)
  END IF
 END DO
END DO

! Second sub-interval:

CALL Coef_D

Fe_old= Fe
Auxx= DTau/4./Dux
Auxz= DTau/4./Duz
Auxxx= DTau/4./Dux/Dux
Auxxz= DTau/8./Dux/Duz
Auxzz= DTau/4./Duz/Duz

DO i= 2,nux-1
 ip1= i+1
 im1= i-1
 DO k= 2,nuz-1
  kp1= k+1
  km1= k-1
  Alpha2(k)= Az(i,km1)*Auxz - (Dzz(i,km1)+Dzz(i,k))*Auxzz
  Beta2(k)= 1 + (Dzz(i,km1)+2.*Dzz(i,k)+Dzz(i,kp1))*Auxzz
  Gamma2(k)= - Az(i,kp1)*Auxz - (Dzz(i,k)+Dzz(i,kp1))*Auxzz
  Psi2(k)= Fe(i,k) + (Ax(ip1,k)*Fe(ip1,k)-Ax(im1,k)*Fe(im1,k))*Auxx &
	+( Dxz(ip1,k)*(Fe(ip1,kp1)-Fe(ip1,km1)) & 
	- Dxz(im1,k)*(Fe(im1,kp1)-Fe(im1,km1)) )*Auxxz &
	+( Dzx(i,kp1)*(Fe(ip1,kp1)-Fe(im1,kp1)) &
	- Dzx(i,km1)*(Fe(ip1,km1)-Fe(im1,km1)) )*Auxxz &
	+( (Dxx(i,k)+Dxx(ip1,k))*(Fe(ip1,k)-Fe(i,k)) &
	- (Dxx(im1,k)+Dxx(i,k))*(Fe(i,k)-Fe(im1,k)) )*Auxxx
 END DO

 CALL Tridag(2,nuz-1,Alpha2,Beta2,Gamma2,Psi2,Faux2,Baux2,Gaux2)

 SELECT CASE(BoundaryCondition)
  CASE("ConstantF ")
   ! Boundary conditions: F=constant at Uz=-Ulim and Uz=Ulim.
   Faux2(1)= Fe(i,1)
   Faux2(nuz)= Fe(i,nuz)
  CASE("ConstantDF")
   ! Boundary conditions in 'Uz', (dF/dz)= constant, at Uz=-Ulim and Uz=Ulim.
   Aux= Faux2(2)-Fepzm(i)*Duz
   IF (Aux<0.) THEN
    Faux2(1)= 0.
   ELSE
    Faux2(1)= Aux
   END IF
   Aux= Faux2(nuz-1)+Fepzp(i)*Duz
   IF (Aux<0.) THEN
    Faux2(nuz)= 0.
   ELSE
    Faux2(nuz)= Aux
   END IF
  CASE DEFAULT
   OPEN(98,FILE='Warning_Adi3p.wt')
   WRITE(98,*) ' BoundaryCondition= ',BoundaryCondition
   WRITE(98,*) ' Boundary Condition must be (ConstantF ) or (ConstantDF) !!'
   CLOSE(98)
   STOP
 END SELECT

 ! Correcting instabilities:
 CALL Cor_Ampli(Faux2,nuz)

 DO k= 1,nuz
  Fe(i,k)= Faux2(k)
 END DO
END DO

SELECT CASE(BoundaryCondition)
 CASE("ConstantF ")
  ! The boundary conditions at 'Uz=Ulim' are automatically taken into account.
  ! For 'Ux=0' we assume zero derivative (approximated forward derivative):
  DO k= 1,nuz
   Fe(1,k)= Fe(2,k)
  END DO
 CASE("ConstantDF")
  ! Boundary conditions: (dF/dx)=0 at Ux=0, (dF/dx)=constant at Ux=Ulim.
  DO k= 1,nuz
   Fe(1,k)= Fe(2,k)
   Aux= Fe(nux-1,k)+Fepxp(k)*Dux
   IF (Aux < 0.) THEN
    Fe(nux,k)= 0.
   ELSE
    Fe(nux,k)= Aux 
   END IF
  END DO
 CASE DEFAULT
  OPEN(98,FILE='Warning_Adi3p.wt')
  WRITE(98,*) ' BoundaryCondition= ',BoundaryCondition
  WRITE(98,*) ' Boundary Condition must be (ConstantF ) or (ConstantDF) !!'
  CLOSE(98)
  STOP
END SELECT

DO k= 1,nuz
 Uz= VUz(k)
 DO i= 1,nux
  Ux= VUx(i)
  U2= Ux**2+Uz**2
  IF (U2 >= Ucrit**2) THEN
   Fe_old(i,k)= Fe(i,k)
  ELSE
   Fe(i,k)= Fe_old(i,k)
  END IF
 END DO
END DO

RETURN
END SUBROUTINE Adi3p

!
! Subroutine for solving a system of linear simultaneous
! equations having a tridiagonal coefficient matrix.
! The equations are numbered from If through Lim, and their
! sub-diagonal, diagonal, and super-diagonal coefficients 
! are stored in the arrays A, B, and C. The computed solution
! vector V(If)...V(Lim) is stored in the array V.
!* Carnahan, Luther and Wilkes, Applied Numerical Methods,
!  John Wiley, 1969, pag. 446.
!
!* Mod. Oct/90: Beta and Gamma appear in the variable list (LFZ).
!* Mod. Aug/06: Adapted to Fortran 95.
!
SUBROUTINE Tridag(If,k,A,B,C,D,V,Beta,Gamma)
! 	PARAMETER(K=100)
IMPLICIT NONE
INTEGER, INTENT(IN) :: If,k 
INTEGER :: Lim,Ifp1,Last,i,j
REAL, INTENT(in) :: A(k),B(k),C(k),D(k)
REAL, INTENT(out) :: V(k)
REAL :: Beta(k),Gamma(k)
!
Lim= k
!... Compute intermediate arrays Beta and Gamma ...
Beta(If)=B(If)
Gamma(If)=D(If)/Beta(If)
Ifp1=If+1
DO i=Ifp1,Lim
 Beta(i)=B(i)-A(i)*C(i-1)/Beta(i-1)
 Gamma(i)=(D(i)-A(i)*Gamma(i-1))/Beta(i)
END DO
!
!... Compute final solution vector V .....
V(Lim)=Gamma(Lim)
Last=Lim-If
DO j=1,Last
 i=Lim-j
 V(i)=Gamma(i)-C(i)*V(i+1)/Beta(i)
END DO
RETURN
END SUBROUTINE Tridag

SUBROUTINE Cor_Ampli(Vf,n)
USE Math_Constants
IMPLICIT NONE
INTEGER :: n,i
REAL, DIMENSION (n) :: Vf(n)

! Corrects negative amplitudes when they are not acceptable.
! For a given real vector, replaces negative values by zero.
! Version Fortran 95, Aug 2006.

DO i= 1,n
 IF ( Vf(i) < 0. ) THEN
  Vf(i)= 0. 
 ELSE
 END IF
END DO

RETURN
END SUBROUTINE Cor_Ampli

SUBROUTINE Aitp1d(nx,Vx,Fx,Xp,Fp)
! Uses linear interpolation in order to obtain the value of a function
! F, at the point Xp. The function F is given as a set of points Fx(x).
! Uses subroutine Locate (Numerical Recipes, P. 96)
! Version Fortran 95, Aug, 2006.

IMPLICIT NONE
INTEGER :: nx,i
REAL, DIMENSION(nx) :: Vx,Fx
REAL :: Xp,Fp,Aux

CALL Locate(Vx,nx,Xp,i)

IF ( i==0 .OR. i==nx ) THEN
 FP= 0.
ELSE
 Aux= ( Xp-Vx(i) )/( Vx(i+1)-Vx(i) )
 Fp= Fx(i) + ( Fx(i+1)-Fx(i) )*Aux
END IF 
RETURN
END SUBROUTINE Aitp1d

SUBROUTINE Aitp2d(nx,ny,Vx,Vy,Fxy,Xp,Yp,Fp)
! Uses linear interpolation in order to obtain the value of a function
! F, at the point (Xp,Yp). 
! The function F is given as a set of points Fxy(x,y).
! Uses subroutine Locate (Numerical Recipes, P. 96)
! Version Fortran 95, Sep., 2006.

IMPLICIT NONE
REAL, DIMENSION(nx) :: Vx
REAL, DIMENSION(ny) :: Vy
REAL, DIMENSION(nx,ny) :: Fxy
REAL :: Xp,Yp,Fp
REAL :: F1,F2,F3,F4,T,U
INTEGER :: nx,ny,i,j

CALL Locate(Vx,nx,Xp,i)
CALL Locate(Vy,nx,Yp,j)

IF ( (i==0 .OR. i==nx) .OR. (j==0 .OR. j==ny) ) THEN
 FP= 0.
ELSE
 F1= Fxy(i,j)
 F2= Fxy(i+1,j)
 F3= Fxy(i+1,j+1)
 F4= Fxy(i,j+1)
 T= ( Xp-Vx(i) )/( Vx(i+1)-Vx(i) )
 U= ( Yp-Vy(j) )/( Vy(j+1)-Vy(j) )
 FP= (1.-T)*(1.-U)*F1 + T*(1.-U)*F2 + T*U*F3 + (1.-T)*U*F4
END IF 
RETURN
END SUBROUTINE Aitp2d

SUBROUTINE Locate(Xx,N,X,J)
! Given an array Xx of lenght N, and given a value X, returns a value 
! J such that X is between Xx(J) and Xx(J+1).
! Xx must be monotonic, either increasing or decreasing.
! J=0 or J=N is returned to indicate that X is out of range.
! See NUMERICAL RECIPES.
! Version Fortran 95, Aug 2006.

IMPLICIT NONE
INTEGER :: N,J,JL,JU,JM
REAL, DIMENSION (N) :: Xx
REAL :: X

JL= 0
JU= N+1
10 IF ( JU-JL .GT. 1 ) THEN
    JM= ( JU+JL ) / 2
    IF ( (Xx(N).GT.Xx(1)).EQV.(X.GT.Xx(JM)) ) THEN
     JL= JM
    ELSE
     JU= JM
    END IF
    GO TO 10
   END IF
   J= JL
RETURN
END SUBROUTINE Locate
 
SUBROUTINE Simpson(Vx,F,N,Res)
! Version Fortran 95, Aug 2006.
IMPLICIT NONE
INTEGER :: N,Nm1,Nm2,i
REAL, DIMENSION(N) :: F,Vx
REAL :: Res,H

! N must be odd.
! The intervals must be equal.

H= ( Vx(N) - Vx(1) )/(N-1)
Res=0.
IF ( N < 5 ) THEN
 OPEN(1,FILE='Warning_Simpson.wt')
 WRITE(1,*) ' N must be larger or equal 5, and odd! '
 CLOSE(1)
 STOP
ELSE
END IF
Nm1= N-1
Nm2= N-2
DO i= 2,Nm1,2
 Res= Res + 4.*F(i)
END DO
DO i= 3,Nm2,2
 Res= Res + 2.*F(i)
END DO
Res= Res + F(1) + F(N)
Res= Res*H/3.
RETURN

END SUBROUTINE Simpson

!
SUBROUTINE Derivxy5p(nx,ny,Vx,Vy,Fxy,Dfdx,Dfdy)
! Version Fortran 95, Sept. 2006.
! Uses 5-point derivative for the internal points.
! See Abramowitz & Stegun, 1970, Eq. 25.3.6, for 'p=0'
! Based on "deriv5p.f", College Park, 16/Feb/2001 (L. F. Ziebell)
! Needs equally spaced points!

! Evaluates the x-derivative and the y derivative of two-dimensional
! function Fxy(x,y), given by an array on (nx,ny) elements

IMPLICIT NONE
INTEGER :: i,j,nx,ny,ip1,im1,ip2,im2,jp1,jm1,jp2,jm2,nxm1,nym1,nxm2,nym2
REAL, DIMENSION(nx) :: Vx
REAL, DIMENSION(ny) :: Vy
REAL, DIMENSION(nx,ny) :: Fxy(nx,ny),Dfdx(nx,ny),Dfdy(nx,ny)
REAL :: Dx,Dy

Dx= Vx(2)-Vx(1)
Dy= Vy(2)-Vy(1)

DO j= 1,ny
 nxm1= nx-1
 nxm2= nx-2
 Dfdx(1,j)= (Fxy(2,j)-Fxy(1,j))/Dx
 Dfdx(2,j)= (Fxy(3,j)-Fxy(1,j))/(2.*Dx)
 Dfdx(nx,j)= (Fxy(nx,j)-Fxy(nxm1,j))/Dx
 Dfdx(nxm1,j)= (Fxy(nx,j)-Fxy(nxm2,j))/(2.*Dx)
 DO i = 3,nxm2
  ip2= i+2
  ip1= i+1
  im1= i-1
  im2= i-2
  Dfdx(i,j)= ( (Fxy(ip1,j)-Fxy(im1,j))*2./3. &
	- (Fxy(ip2,j)-Fxy(im2,j))/12. ) / Dx
 END DO
END DO

Do i= 1,nx
 nym1= ny-1
 nym2= ny-2
 Dfdy(i,1)= (Fxy(i,2)-Fxy(i,1))/Dy
 Dfdy(i,2)= (Fxy(i,3)-Fxy(i,1))/(2.*Dy)
 Dfdy(i,ny)= (Fxy(i,ny)-Fxy(i,nym1))/Dy
 Dfdy(i,nym1)= (Fxy(i,ny)-Fxy(i,nym2))/(2.*Dy)
 Do j = 3,nym2
  jp2= j+2
  jp1= j+1
  jm1= j-1
  jm2= j-2
  Dfdy(i,j)= ( (Fxy(i,jp1)-Fxy(i,jm1))*2./3. &
	- (Fxy(i,jp2)-Fxy(i,jm2))/12. ) / Dy
 END DO
END DO

RETURN
END SUBROUTINE Derivxy5p
!

SUBROUTINE RKQC(Y,DYDX,N,X,HTRY,EPS,YSCAL,HDID,HNEXT)
! See NUMERICAL RECIPES.
! Version Fortran 95, Feb. 2007.
IMPLICIT NONE
! EXTERNAL DERIVS
INTEGER, PARAMETER :: NMAX=15000
REAL, PARAMETER :: FCOR=.0666666667,SAFETY=0.9,ERRCON=6.E-4
INTEGER, INTENT(in) :: N
INTEGER :: I
REAL :: X,HTRY,EPS,HDID,HNEXT
REAL :: PGROW,PSHRNK,XSAV,H,HH,ERRMAX
REAL, DIMENSION(N) :: Y,DYDX,YSCAL
REAL, DIMENSION(NMAX) :: YTEMP,YSAV,DYSAV

IF (NMAX<N) THEN
 OPEN(98,FILE='Warning_Rkqc.wt')
 WRITE(98,*) ' Nmax must be smaller than N, in RKQC and RK4.'
 CLOSE(98)
 STOP
ELSE
END IF
PGROW=-0.20
PSHRNK=-0.25
XSAV=X
DO I=1,N
 YSAV(I)=Y(I)
 DYSAV(I)=DYDX(I)
END DO
H=HTRY
1 HH=0.5*H
CALL RK4(YSAV,DYSAV,N,XSAV,HH,YTEMP)
X=XSAV+HH
CALL DERIVS(X,YTEMP,DYDX)
CALL RK4(YTEMP,DYDX,N,X,HH,Y)
X=XSAV+H
IF (X == XSAV) THEN 
  OPEN(98,FILE='Warning_Rkqc.wt')
  WRITE(98,*) ' Stepsize not significant in RKQC.'
  CLOSE(98)
  STOP
ELSE
END IF 
CALL RK4(YSAV,DYSAV,N,XSAV,H,YTEMP)
ERRMAX=0.
DO I=1,N
 YTEMP(I)=Y(I)-YTEMP(I)
 ERRMAX= MAX1(ERRMAX,ABS(YTEMP(I)/YSCAL(I)))
END DO
ERRMAX=ERRMAX/EPS
IF(ERRMAX > 1.) THEN
 H=SAFETY*H*(ERRMAX**PSHRNK)
 GOTO 1
ELSE
 HDID=H
 IF(ERRMAX > ERRCON)THEN
  HNEXT=SAFETY*H*(ERRMAX**PGROW)
 ELSE
  HNEXT=4.*H
 ENDIF
ENDIF
DO I=1,N
 Y(I)=Y(I)+YTEMP(I)*FCOR
END DO
RETURN
END SUBROUTINE RKQC

SUBROUTINE RK4(Y,DYDX,N,X,H,YOUT)
! Given values for N variables Y and their derivatives DYDX
! known at X, use the fourth-order Runge-Kutta method to advan-
! ce the solution over an interval H and return the incremen-
! ted variables as YOUT, which need not be a distinct array from
! Y. The user supplies the subroutine DERIVS(X,Y,DYDX) which
! returns derivatives DYDX at X.
! See NUMERICAL RECIPES.
! Version Fortran 95, Feb. 2007.

IMPLICIT NONE
INTEGER, PARAMETER :: NMAX=15000
INTEGER, INTENT(in) :: N
INTEGER :: I
REAL, DIMENSION(N) :: Y,DYDX,YOUT
REAL, DIMENSION(NMAX) :: YT,DYT,DYM
REAL :: X,H,HH,H6,XH
        
HH=H*0.5
H6=H/6.
XH=X+HH
! CALL  DERIVS(X,Y,DYDX)
DO I=1,N
 YT(I)=Y(I)+HH*DYDX(I)
END DO
CALL DERIVS(XH,YT,DYT)
DO I=1,N
 YT(I)=Y(I)+HH*DYT(I)
END DO
CALL DERIVS(XH,YT,DYM)
DO I=1,N
 YT(I)=Y(I)+H*DYM(I)
 DYM(I)=DYT(I)+DYM(I)
END DO
CALL DERIVS(X+H,YT,DYT)
DO I=1,N
 YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.*DYM(I))
END DO
RETURN
END

SUBROUTINE Derivs(Tau,Ft,Dfdt)

USE Common_Arrays
USE Common_Params
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL, DIMENSION(4*nqx*nqz) :: Ft,Dfdt
REAL, INTENT(in) :: Tau
REAL :: Qx,Qz,Q2
REAL :: Aux
REAL, DIMENSION(nux,nuz) :: Dfdux,Dfduz
INTEGER :: ind
INTEGER :: i,k,sigma

CALL Cor_Ampli(Ft,4*nqx*nqz)
CALL Derivxy5p(nux,nuz,VUx,VUz,Fe,Dfdux,Dfduz)

ind= 0
DO i= 1,nqx
 DO k= 1,nqz
  ind= ind+1
  ILp(i,k)= Ft(ind)
 END DO
END DO
DO i= 1,nqx
 DO k= 1,nqz
  ind= ind+1
  ILm(i,k)= Ft(ind)
 END DO
END DO
DO i= 1,nqx
 DO k= 1,nqz
  ind= ind+1
  ISp(i,k)= Ft(ind)
 END DO
END DO
DO i= 1,nqx
 DO k= 1,nqz
  ind= ind+1
  ISm(i,k)= Ft(ind)
 END DO
END DO

ind= 0
DO i= 1,nqx
 Qx= VQx(i)
 DO k= 1,nqz
  Qz= VQz(k)
  Q2= Qx**2+Qz**2
  ind= ind+1
  IF (Q2<EpsMin) THEN       ! Avoid growth of waves with Q near zero
   Dfdt(ind)= 0.
  ELSE
   sigma= 1
   CALL Coef_Lwave(Qx,Qz,i,k,sigma,Dfdux,Dfduz,ILp,ILm,Aux)
   Dfdt(ind)= Aux
  END IF
 END DO
END DO
DO i= 1,nqx
 Qx= VQx(i)
 DO k= 1,nqz
  Qz= VQz(k)
  Q2= Qx**2+Qz**2
  ind= ind+1
  IF (Q2<EpsMin) THEN       ! Avoid growth of waves with Q near zero
   Dfdt(ind)= 0.
  ELSE
   sigma= -1
   CALL Coef_Lwave(Qx,Qz,i,k,sigma,Dfdux,Dfduz,ILm,ILp,Aux)
   Dfdt(ind)= Aux
  END IF
 END DO
END DO
DO i= 1,nqx
 Qx= VQx(i)
 DO k= 1,nqz
  Qz= VQz(k)
  Q2= Qx**2+Qz**2
  ind= ind+1
  IF (Q2<EpsMin) THEN       ! Avoid growth of waves with Q near zero
   Dfdt(ind)= 0.
  ELSE
   sigma= 1
   CALL Coef_Swave(Qx,Qz,i,k,sigma,Dfdux,Dfduz,ISp,Aux)
   Dfdt(ind)= Aux
  END IF
 END DO
END DO
DO i= 1,nqx
 Qx= VQx(i)
 DO k= 1,nqz
  Qz= VQz(k)
  Q2= Qx**2+Qz**2
  ind= ind+1
  IF (Q2<EpsMin) THEN       ! Avoid growth of waves with Q near zero
   Dfdt(ind)= 0.
  ELSE
   sigma= -1
   CALL Coef_Swave(Qx,Qz,i,k,sigma,Dfdux,Dfduz,ISm,Aux)
   Dfdt(ind)= Aux
  END IF
 END DO
END DO

RETURN
END SUBROUTINE Derivs

