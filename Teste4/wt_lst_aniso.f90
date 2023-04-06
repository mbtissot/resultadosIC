! Comments:
! 110205: Program to solve the time evolution of the equations for
!       electrons, Langmuir waves, ion-sound waves, and transverse waves, in
!       homogeneous medium, as formulated in notes "B11", dated Feb. 05, 2011.
!       Based on "wt_2d1d.f90", version 100826.
!       Uses subroutine "Split", with splitting method to solve the
!       equation for "Fe". 
!       The initial version does not yet include transverse waves, and should
!       be equivalent to "wt_2d1d.f90".
! 110219: Includes terms related to transverse waves, and the equation for
!       T waves.
! 110326: Corrects some mistakes found in notes B11:
!       1) L waves: in the scattering term related to L and T waves, eq. (43)
!       2) S waves: resonance conditions of decay involving L and T waves
!               (definitions of q1 ... q6).
!       3) T waves: in the scattering term related to L and T waves, eq. (50)
! 110403: Modifies the formulation according to notes E11.
! 110412: Modifies the formulation according to notes F11. These notes modify
!       the evaluation of the scattering terms, and correct a mistake in the
!       decay term involving L and S waves, in the equation for S waves. This
!       mistake was present in all previous notes.
! 110413: Modifies the dimensions of many arrays; instead of "-1:1" we use
!       dimension "2", creating vector "Iaux(2)" containing values "1" and "-1".
!       Warning: this version never worked properly. Although it compiles
!       correctly, it gives "segmentation errors", not identified, when it
!       is utilized. The version has been abandoned, and version 110425
!       was developed.
! 110425: Based on notes F11, as version 110412. The difference regarding that
!	version is that some quantities which were evaluated in subroutine
!	"Res_Cond" in version 110412 are evaluated in subroutines
!	"Coef_Lwave", "Coef_Swave" and "Coef_Twave". This is less demanding
!	in terms of allocated memory. Some tests have shown that the running
!	time of the code is not very different.
!	Introduced the option to include or not the scattering by electrons
!	(only spontaneous, see notes F-11).
!	This version corrects a mistake in the spontaneous terms for S waves
!	(which caused a small increase of the spectrum near the line Qx=0).
! 110630: Introduces a small modification in the use of boundary conditions.
!       The subroutine 'Boundary_Condition' is called after each of the
!       operators, Lx and Lz, in subroutine 'Split', instead of only at the
!       end, as in the previous versions. The subroutine is also modified
!       in the case of 'ConstantDF': When the value at the edge is negative,
!       keeps the previous value instead of attributing zero value.
!       The use of allocatable arrays is introduced, so that the dimensions
!       nux, nuz, nqx, nqz and nrz are input parameters of the code.
! 110723: Introduces the Runge-Kutta method for the wave equations, using
!       subroutine "RK4" (as in "wt_2d"). The present version considers
!       only the case of fixed time step. 
! 111115: Based on notex J11. Some mistakes found in previous versions of
! 	the notes have been corrected, and some features have been improved.
!	Particularly, some resonance conditions have been improved.
! 111211: Introduces some improvements regarding version 111115. One of them
!       is that the resonance values of "Qxp" or "Qzp" and of "Uxp" or "Uzp"
!       are evaluated only in subroutine "Res_Cond" and preserved for use
!       along the time evolution.
! 120302: Introduces a new approach for the scattering terms, following
!	notes C12 and "12b_comm", by Rudi Gaelzer. 
! 120325: Introduces an improvement in the new approach for the scattering 
!	terms, adding the contribution of induced scattering due to ions, 
!       following notes C12 dated March 11, 2012. 
! 120514: Introduces improvements on the resonance conditions for decay terms,
!       following notes D12 dated May 16, 2012.
! 120810: Further improvements on the resonance conditions for decay terms 
!       (particularly TdecayLL) and on scattering terms, following notes E12 
!       dated August 08, 2012.
! 120818: Introduces finite difference equations for the extreme points in
!       velocity space (except Ux=0, where the condition of zero derivative
!       is utilized). The subroutines for boundary conditions are not
!       necessary anymore, as well as the variable "BoundaryCondition" 
!       (see notes E12, dated August 18, 2012).
! 121229: Introduces improvements on the resonance conditions for decay terms,
!       following notes K12 dated Oct 22, 2012.
! 130121: Introduces improvements on the resonance conditions for decay terms,
!       following notes A13 dated Jan 21, 2013. Some resonance conditions
!       are squared, which can introduce spurious roots. A procedure to verify
!       which are the real roots is introduced. Preliminary results are not
!       encouraging, and the version has been abandoned, at least temporarily.
! 130126: Introduces different procedures for verification of the real roots.
!       Preliminary results are not very encouragin, and the version has been
!       abandoned, at least temporarily.
! 130128: The roots of the resonance conditions are found by a numerical 
!       procedure, which does not require squaring and generation of spurious
!       roots. Subroutine 'RTSAFE', from Numerical Recipes, is utilized to
!       locate the roots. 
! 130301: Some modifications have been made to the interpolation procedures
!       in subroutines "Funcx_RcdXXX" and "Funcz_RcdXXX". 
!	Another modification has been in the indexes of resonant values of
!	Qxdif and Qzdif, in the decay terms. The index "iresdif2" has been
!	introduced.
!       A mistake in the coefficients of the term related to decay involving
!       S and T waves ("CoefLSTdA" and "CoefLSTdB") in the subroutine
!       "Coef_Lwave": the coefficients are now multiplied by "AA", as they 
!       should be.
! 130306: Part of the coefficients of decay terms are evaluated in the
!	subroutine "Aux_Res_Cond" and stored as arrays "CxpLLSd", "CzpLLSd",
!	etc., instead of evaluated at each iteration.
!	Another modification is that the initial spectra are evaluated by
!	separated subroutines, "ILwave_Init", "ISwave_Init" and "ITwave_Init",
!	instead of a single subroutine as in previous versions.
!	This version was intended to be faster than version 130301. However,
!	the tests have shown that version 130301 is faster (about 75% of the
!	running time). Therefore, this version will not be used. We will
!	keep using version 130301.
! 130320: This version starts from version 130301; however, the roots of the
!       resonance conditions in the decay terms are now searched using
!       subroutine "RTBIS" (NUmerical Recipes), which uses the method of
!       bisection instead of the Newton-Raphson and bisection as in "RTSAFE".
!       The auxiliary routines ZBRAC2 and ZBRAC1 have been introduced. The
!       subroutines "Funcx_RcdXXX" and "Funcz_RcdXXX" have been transformed
!       into "functions". 
! 130322: Version based on version 130320. It eliminates all files related to
!       the decision based on "Qz>Qx", for decay coefficients. For decay terms
!       we look for the resonant value of "Qxp", using routines ZBRAC2 and
!       RTBIS. This reduces the amount of memory allocated to the resonance
!       conditions.
! 130920: Version based on version 130322. It gives options for the starting
!       condition of the T waves. In addition to the options existing before,
!       the new option assumes that the initial spectrum of T waves is given
!       by a condition of "turbulent equilibrium" (see notes F-12).
! 140404: The terms associated to spontaneous and induced emission for L waves
!       are evaluated using a grid with (nqx2,nqz2) points. The values for the
!       grid with (nqx,nqz) points, used for all other coefficients, are
!       obtained using interpolation. Experimentation has shown that the
!       scattering terms must be evaluated with large number of points, like
!       (71,71), while "Lemis" does not behave well for grids larger than
!       (51,51).
! 140423: The term associated to decay of L waves, involving LS waves, is also
!       evaluated using the (nqx2,nqz2) grid, and then added to the QL term and 
!       interpolated into the (nqx,nqz) grid. The results obtained are good,
!       but the code runs in a time nearly twice the time of the previous
!       version, probably due to repeated interpolation.
! 140424: All decay terms, as well as the QL terms, are evaluated using the
!       (nqx2,nqz2) grid and then interpolated into the (nqx,nqz) grid. The
!       locator indexes for the interpolation are evaluated as the code
!       starts, in subroutine "Res_Cond", stored in the arrays "IQresx2" and
!       "IQresz2", and used along the time evolution. Uses subroutine
!       "Aitp2d2", who utilizes the locator indexes obtained using "Locate".
!       The points of the (nqx,nqz) grid are interpolated into the (nqx2,nqz2)
!       grid, with locators stored into "IQresx", and "IQresz".
!       Another modification:
!       In order to increase the speed of the code, interpolations appearing
!       in the scattering terms are now made with "Aitp2d2", with locating
!       indexes evaluated in subroutine "Res_Cond" and stored to be used along
!       time evolution.
!       Correction of mistakes:
!       1) In the scattering terms, the quantity "Qs" was used before being 
!       defined. The mistake is present in all the most recent versions of the
!       code, in the segments corresponding to induced scattering by the ions. 
!       2) In the scattering for T waves, induced contribution of the ions,
!       the quantity Qstar was defined in previous versions without the SQRT.
!       3) In subroutine "Coef_Twave", the quantity "Phi" was used in the term
!       of scattering for T waves, but was not defined.
! 140505: Incorporates the correction of mistakes in the scattering terms, 
!       which has been made in version 140424. Also incorporates the evaluation
!       of locator indexes for the scattering terms and subsequent use of
!       subroutine "Aitp2d2", as in version 140424.
!       Quasilinear and decay terms evaluated using the (nqx,nqz) grid, also
!       used for the wave spectra appearing in the particle equation (as in
!       version 130920. Scattering terms, which need more resolution, are
!       evaluated using another grid (nqx2,nqz2), and their effect is then
!       interpolated for the points of grid (nqx,nqz).
! 140731: Quasilinear and decay terms for L and S waves continue to be
!       evaluated using the (nqx,nqz) grid. Scattering terms for all types of
!       waves and also the decay terms for T waves are evaluated using the grid
!       with (nqx2,nqz2), and their effect is then interpolated for the points
!       of the grid (nqx,nqz). 
! 140824: Quasilinear and decay terms for L, S, and T waves are evaluated using
!       the (nqx,nqz) grid.
!       Scattering terms for L waves are evaluated using the (nqx,nqz) grid.
!       Scattering terms for T waves are evaluated using the (nqx2,nqz2) grid,
!       and their effect is then interpolated for the points of grid (nqx,nqz).
!       Scattering terms for S waves are not taken into account.
!       This version is in many regards similar to version 140505, but it
!       incorporates some improvements adopted in version 140731, particularly
!       regarding the interpolation procedure near "Q=0".
!       Subroutine "Rebuild" is used only for "Fe", "IL", and "IS". The
!       spectrum of T waves, which can be very peaked near "Q=0" in the case
!       with a beam, is not smoothed out in this version.
!       Introduces parameter "nph", for the angular variable "Phip". It is
!       evaluates as the maximum between "nqz" and "nqz2".
! 140828: This version is basically the same as version 140731, with some of
!       the changes incorporated to version 140824:
!       Quasilinear and decay terms for L and S waves are evaluated using
!       the (nqx,nqz) grid.
!       Scattering terms for L waves are also evaluated using the (nqx,nqz) 
!       grid. Scattering terms for S waves are not taken into account.
!       Scattering and decay terms for T waves are evaluated using the
!       (nqx2,nqz2) grid, and their effect is then interpolated for the 
!       points of grid (nqx,nqz).
! 141015: This version introduces subroutine "Energy2", which evaluates the
!       energy density of waves L, S, and T, along the time evolution. It
!       also evaluates the energy associated to the fundamental and the
!       harmonics 2 and 3 of the T wave emission. These results appear in
!       file "Ewave.wt".
!       In subroutine "Output", it is now possible to generate file "ITQ1",
!       which contains the integral of the T wave intensity over angle dPhi,
!       and file "ITQ2", which contains the integral over "Q*dPhi".
! 150225: Corrects a sign mistake in components "xz", in "Coef_Coll" (this
!       subroutine introduces collisional effects, although it has not yet
!       been sufficiently tested).
!       For the collision coefficients, see page 7 of notes B15, dated 
!       Apr. 28, 2015 (equation (16), section 4). See also
!       page 27 of notes D08, dated Nov. 15, 2013.
!       Introduces character parameters "TimeEvol", "OneDfiles", "TwoDfiles",
!       "Onecolumnfiles", and "ADcoefficients", which allow to decide which
!       sets of output files to be generated, and if the time evolution will
!       be performed or if only files will be generated, from previous output
!       files. These parameters must be provided in file "Start.wt"
! 151126: This is a special version of program "wt_lst.f90", adapted to 
!       produce results for the collisional damping of L and S waves, for the 
!       case of a Maxwellian distribution, using the formalism presented in 
!       "Weak turbulence theory for collisional plasmas", by P. H. Yoon, L. F.
!       Ziebell, and E. P. Kontar, submitted for publication in September 2015.
!       Uses subroutine "Coll_Damping" to generate the collisional damping and
!       the Landau damping as well, for Maxwellian distribution
!       After printing the damping coefficients, there is a "STOP", to finish
!       the execution of the code, before the time evolution of the weak 
!       turbulence equations.
! 151208: Introduces the use of the collisional damping (approximated assuming
!       Maxwellian distribution) into the equations for the time evolution of
!       L and S waves.
!       Another modification: 
!       The subroutine "Rebuild" has been divided into three subroutines,
!       "Rebuild_L", "Rebuild_S" and "Rebuild_Fe". There are three "CALL"
!       commands, one for each. It is therefore possible to choose if the
!       spectra for L and S waves, and also Fe, will be smoothed out by
!       "Rebuild" prior creating the output files, or not. The choice is made
!       by commenting or not the CALL commands, before compiling the code. 
!       There is no logical parameter for the decision.
!       This modification is the same is introduced in "wt_lstn.f90", version
!       151205.
!       Corrects a mistake in subroutine "Aitp2d" (which only affected the
!       results in the case of use with a non-square grid, with ny diff nx).
! 170328: Program to solve the time evolution of the equations for
!       electrons, Langmuir waves, ion-sound waves, and transverse waves, in
!       homogeneous medium, with ions described by Maxwellian distribution 
!       and electrons described by a Maxwellian background and an halo of
!       type "kappa" distributions. 
!       See notes "C17", dated March 28, 2017.
!       See notes "D17", dated March 30, 2017.
!       Based on "wt_lst.f90", version 151208, and in "wt_lstk.f90",
!       version 170319.
!       The present versions includes:
!       L waves: spontaneous and induced emission, decay involving L and S
!                waves, and scattering involving L waves.
!       S waves: spontaneous and induced emission, decay involving L and S
!                waves, and decay involving L and T waves.
!       T waves: decay involving only L waves, decay involving L and S
!                waves, decay involving L and T waves, and scattering 
!                involving L and T waves. 
!       Uses subroutine "Split", with splitting method to solve the
!       equation for "Fe". 
!       Decay terms are the same as in "wt_lst.f90", version 151208.
!       Subroutine "Iwave_Init" has been modified. It now provides L, S, and
!       T spectra, depending on the choice of a parameter "WaveT".
!       The distribution functions and the wave spectra are obtained 
!       considering a 2D geometry. 
!       Includes the evaluation of collisional damping and electrostatic
!       bremsstrahlung, and introduces the use of these quantities into the
!       equations for the time evolution of L and S waves.
!       The collisional damping and the electrostatic bremsstrahlung can be
!       taken into accout into the evaluation of the initial
!       spectra of L and S waves (see notes C17, dated 28 March, 2017).
! 170430: Modify the boundary conditions in subroutine "Split": For "U<Ucrit",
!       Fnew(1,k)=Fe(i,k), and for "U>Ucrit", Fe(1,k)=Fe(2,k). This is to
!       avoid artificial flattening of the top of the electron distribution,
!       and at the same time to guarantee zero derivative at the boundary,
!       in the tail region, more affected by evolution.
!       Introduces parameter "Lsaturation", which allows to decide if the
!       initial spectrum of L waves can diverge at "Q=0", or if it will be
!       considered saturated for omega/k>c. It is applied also for evaluation
!       of the bremsstrahlung and for the collisional damping (for L waves).
!       Introduces parameters "RebuildL", "RebuildS", "RebuildFe", which allows
!       to choose if we apply "Rebuild" before saving the data, to smooth out
!       irregularities.
!       Introduces parameters "NewEffects1" (to include Bremsstrahlung and
!       collisional damping in the evaluation of the initial spectra), and
!       "NewEffects2" (to include these new effects in the equations for the
!       time evolution of the waves).
! 170819: Some changes in the evaluation of the Ai and the Dij:
!       In "Coef_A", uses the resonant values evaluated in "Res_Cond", instead
!       of calculating these values again.
!       In "Coef_D", uses the initial spectrum of L waves for resonant points
!       outside of the (Qx,Qz) grid. In previous versions, the Dij were 
!       assumed to vanish outside of the grid.
!       Correction of mistake in Mar 09, 2020: In four calls to Simpson, at
!       subroutine "Coef_Twave", in the term with "TscatLT", it 
!       was written VQx and VQz, instead of VQx3 and VQz3. Corrected, 
!       but not yet tested. (this correction was based in a correction made
!       in Jan 11, 2018, in version 171115)
! 200309: Version based on version 170819. Uses also some modifications 
!       introduced in version 171115, in the decay terms: the subroutine
!       Aux_Coef_Decay_z is as in version 171115.
!       In comparison with previous versions, uses a different expression for 
!       evaluation of the electrostatic bremsstrahlung for L waves, based on 
!       the "new" expression appearing at page 19 of notes C-16, dated April 
!       04, 2019, right before equation (31).
!       The code is transformed to DOUBLE PRECISION, by changing REAL into
!       REAL*8, COMPLEX into COMPLEX*16, and by changing 0. to 0.D0 and 1. to
!       1.D0, in the argument of functions.
! 200427: Some changes in the evaluation of resonance conditions. 
!       See file "dif_200427_200309"
! 210105: Introduces the possibility of use of ring-beam distributions as
!       initial state, for electrons. See notes B13, dated Feb. 08, 2016, 
!	equation (30).
!       This modification is based at a version developed along with Larissa
!       Petruzzellis, found in 
!       ~/zb/work/larissa/tst_2beams
!       The version was named "160210"
 

!---------------------------------------------------------------------
! List of modules, subroutines, and functions:

! MODULE Common_Params
! MODULE Common_Arrays
! MODULE Allocate_Arrays
! MODULE Math_Constants
! MODULE Phys_Constants
! MODULE Sub_Prog
! MODULE hyp_2f1_module

! SUBROUTINE Definitions
! SUBROUTINE Space_Profiles
! SUBROUTINE Init_Wave(Tau,Dqx,Dqz,Dux,Duz,Anorm,Epart,Ewave,EppEw,&
!     EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3)
! SUBROUTINE Iwave_Init(Qx,Qz,Iinit,WaveT,HBrL,HBrS,HGcL,HGcS)
! SUBROUTINE Fe_Init(Ux,Uz,Fesum,Femax,Febkgr,Fef,Feb)
! SUBROUTINE Save_Results(Tau,Dqx,Dqz,Dux,Duz,Anorm0,Epart0,Ewave0,EppEw0,&
!     Rn,Rp,Rw,Rs,EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3)
! SUBROUTINE Read_Results(Tau,Dqx,Dqz,Dux,Duz,Anorm0,Epart0,Ewave0,EppEw0,&
!     Rn,Rp,Rw,Rs,EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3)
! SUBROUTINE Res_Cond
! SUBROUTINE Coef_A
! SUBROUTINE Coef_D
! SUBROUTINE Coef_Lwave(sigma,Dfdux,Dfduz,CoefA,CoefB)
! SUBROUTINE Coef_Swave(sigma,Dfdux,Dfduz,CoefA,CoefB)
! SUBROUTINE Coef_Twave(sigma,CoefA,CoefB)
! SUBROUTINE Aux_Coef_Decay_z(nqxa,nqza,ires,iresdif,&
!  Qxp,Qzp,Qxdif,Qzdif,Wavep,Wavedif,VQxa,VQza,Vauxqxp,Vauxqxdif,&
!  Iqp,Iqdif)
! SUBROUTINE Fnorm(Anorm)
! SUBROUTINE Energy(Epart,Ewave,EppEw)
! SUBROUTINE Energy2(EwaveL,EWaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3)
! SUBROUTINE Output_Coef(CoefChoice,CoefA,CoefB)
! SUBROUTINE Output(WriteChoice)
! SUBROUTINE Output2(WriteChoice)
! SUBROUTINE Split(Dux,Duz,DTau)
! SUBROUTINE Evol_Iwave(DTau,WaveType)
! SUBROUTINE Tridag(If,k,A,B,C,D,V,Beta,Gamma)
! SUBROUTINE Cor_Ampli(Vf,n)
! SUBROUTINE Cor_Ampli_2(Vf,Vfnew,n)
! SUBROUTINE Aitp1d2(nx,Vx,Fx,Xp,Fp,i)
! SUBROUTINE Locate(Xx,N,X,J)
! SUBROUTINE Aitp2d(nx,ny,Vx,Vy,Fxy,Xp,Yp,Fp)
! SUBROUTINE Aitp2d2(nx,ny,Vx,Vy,Fxy,Xp,Yp,Fp,i,j)
! SUBROUTINE Aitp2d2b(nx,ny,Vx,Vy,Fxy,Xp,Yp,Fp,i,j)
! SUBROUTINE Simpson(Vx,F,N,Res)
! SUBROUTINE Derivxy5p2d(nx,ny,Vx,Vy,Fxy,Dfdx,Dfdy)
! SUBROUTINE Derivxy5pln2d(nx,ny,Vx,Vy,Fxy,Dfdx,Dfdy)
! SUBROUTINE Derivx5p2d(nx,ny,Vx,Fxy,Dfdx)
! SUBROUTINE Derivy5p2d(nx,ny,Vy,Fxy,Dfdy)
! SUBROUTINE Coef_Coll
! SUBROUTINE Coll_Damp
! SUBROUTINE Bremsstrahlung
! SUBROUTINE Rebuild_L
! SUBROUTINE Rebuild_S
! SUBROUTINE Rebuild_Fe
! SUBROUTINE BezierDDzzduz
! SUBROUTINE RK4(Y,DYDX,N,X,H,YOUT)
! SUBROUTINE DERIVS(Tau,Ft,Dfdt,N)
! SUBROUTINE ZBRAC2(FUNC,X1,X2,SUCCES)

! REAL*8 FUNCTION Funcx_RcdLLS(Qxp)
! REAL*8 FUNCTION Funcx_RcdLLT(Qxp)
! REAL*8 FUNCTION Funcx_RcdLST(Qxp)
! REAL*8 FUNCTION Funcx_RcdLTT(Qxp)
! REAL*8 FUNCTION Funcx_RcdSLL(Qxp)
! REAL*8 FUNCTION Funcx_RcdSLT(Qxp)
! REAL*8 FUNCTION Funcx_RcdTLL(Qxp)
! REAL*8 FUNCTION Funcx_RcdTLS(Qxp)
! REAL*8 FUNCTION Funcx_RcdTTL(Qxp)

! REAL*8 FUNCTION ZL(Qx,Qz)
! REAL*8 FUNCTION ZS(Qx,Qz)
! REAL*8 FUNCTION ZT(Qx,Qz)
! REAL*8 FUNCTION RTBIS(FUNC,X1,X2,XACC) 
! REAL*8 FUNCTION UgL(Qx,Qz)
! REAL*8 FUNCTION UgS(Qx,Qz)
! REAL*8 FUNCTION PERF(X)
! REAL*8 FUNCTION PERFC(X)
! REAL*8 FUNCTION PERFCE(X)
! REAL*8 FUNCTION Aux_GcollL(Qp)
! REAL*8 FUNCTION Aux_GcollS(Qp)
! REAL*8 FUNCTION Aux_Brem(Qp)
! REAL*8 FUNCTION Aux_BremL(Qp)

!---------------------------------------------------------------------

!---------------------------------------------------------------------
! Modules with definitions:
!---------------------------------------------------------------------

MODULE Common_Params
IMPLICIT NONE

! Read parameters:
INTEGER :: nqx,nqz,nux,nuz,nrz  ! nqx,nqz,nux,nuz,nrz must be odd numbers.
INTEGER :: nqx2,nqz2  ! nqx2,nqz2 must be odd numbers.
INTEGER :: nqmod
INTEGER :: AsympT
INTEGER :: Alphai,Alphae 
REAL*8 :: Ulim,Ucrit
REAL*8 :: Qxi,Qxf,Qzi,Qzf 
REAL*8 :: Qi,Qf
REAL*8 :: Iw0 
REAL*8 :: Kappai,Kappae
REAL*8 :: RatioNki,RatioNke
REAL*8 :: RatioNf,Uf,RTfparTs, RTfperpTs,Thf,Ufx,Ufz
REAL*8 :: RatioNb,Ub,RTbparTs, RTbperpTs,Thb,Ubx,Ubz
REAL*8 :: RTeTi,G
REAL*8 :: Ve2C2
REAL*8 :: ScaleDensInhom
CHARACTER(LEN=6) :: InitialLevel
CHARACTER(LEN=3) :: Lemis,LdecayLS,LdecayLT,LdecayST,LdecayTT,LscatLL,LscatLT
CHARACTER(LEN=3) :: Semis,SdecayLL,SdecayLT,Sscat
CHARACTER(LEN=3) :: TdecayLL,TdecayLS,TdecayTL,TscatLT
CHARACTER(LEN=3) :: CollTerm
CHARACTER(LEN=8) :: CollTermForm
CHARACTER(LEN=3) :: SpontEmis
CHARACTER(LEN=3) :: ScatElSpo
!CHARACTER(LEN=3) :: ScatElInd
CHARACTER(LEN=3) :: Gcoll
CHARACTER(LEN=3) :: Bremss
CHARACTER(LEN=3) :: NewEffects1,NewEffects2
CHARACTER(LEN=3) :: RenormFe
CHARACTER(LEN=3) :: DerivLn
CHARACTER(LEN=3) :: RebuildL,RebuildS,RebuildFe
CHARACTER(LEN=3) :: Lsaturation

! Parameters generated in the code:
REAL*8 :: RMiMe
REAL*8 :: AA,Geff 
REAL*8 :: U0x,U0z
REAL*8 :: Ue,Ui,Ue2,Ui2
REAL*8 :: Uek,Uik,Uek2,Uik2
REAL*8 :: Uf2
REAL*8 :: AuxInitialLevel
REAL*8 :: AuxAlpha
INTEGER :: Auxinterp
INTEGER :: nph

REAL*8 :: RatioGammaep051,RatioGammaem051
REAL*8 :: RatioGammae015,RatioGammaem115
REAL*8 :: RatioGammai015,RatioGammaim115

!JPi
REAL*8, PARAMETER :: QxSqr=0.4,QzSqr=0.12
!JPf

END MODULE Common_Params

MODULE Common_Arrays
IMPLICIT NONE

REAL*8, DIMENSION(:), ALLOCATABLE :: VRz
REAL*8, DIMENSION(:), ALLOCATABLE :: VQx
REAL*8, DIMENSION(:), ALLOCATABLE :: VQz
REAL*8, DIMENSION(:), ALLOCATABLE :: VQx2
REAL*8, DIMENSION(:), ALLOCATABLE :: VQz2
REAL*8, DIMENSION(:), ALLOCATABLE :: VQx3
REAL*8, DIMENSION(:), ALLOCATABLE :: VQz3
REAL*8, DIMENSION(:), ALLOCATABLE :: VQQ
REAL*8, DIMENSION(:,:), ALLOCATABLE:: ILp,ILm,ISp,ISm,ITp,ITm
REAL*8, DIMENSION(:,:), ALLOCATABLE:: ILpOld, DeltaILp, ILmOld, DeltaILm
REAL*8, DIMENSION(:,:), ALLOCATABLE:: GcollLp,GqlLp
REAL*8, DIMENSION(:,:), ALLOCATABLE:: GcollLm,GqlLm
REAL*8, DIMENSION(:), ALLOCATABLE:: GcollL1D,GqlL1D
REAL*8, DIMENSION(:,:), ALLOCATABLE:: GcollSp,GqlSp
REAL*8, DIMENSION(:,:), ALLOCATABLE:: GcollSm,GqlSm
REAL*8, DIMENSION(:), ALLOCATABLE:: GcollS1D,GqlS1D
REAL*8, DIMENSION(:,:), ALLOCATABLE:: HGcollLp,HGcollSp
REAL*8, DIMENSION(:), ALLOCATABLE:: HGcollL1D,HGcollS1D
REAL*8, DIMENSION(:,:), ALLOCATABLE:: BremL,BremS
REAL*8, DIMENSION(:), ALLOCATABLE:: BremL1D,BremS1D
REAL*8, DIMENSION(:,:), ALLOCATABLE:: HBremL,HBremS
REAL*8, DIMENSION(:), ALLOCATABLE:: HBremL1D,HBremS1D
REAL*8, DIMENSION(:), ALLOCATABLE :: VPhip
REAL*8, DIMENSION(:), ALLOCATABLE :: VUx
REAL*8, DIMENSION(:), ALLOCATABLE :: VUz
REAL*8, DIMENSION(:,:), ALLOCATABLE :: Fe0,Fe,Fi,FeOld,DeltaFe
REAL*8, DIMENSION(:,:), ALLOCATABLE :: Ax,Az,Dxx,Dxz,Dzx,Dzz
REAL*8, DIMENSION(:,:), ALLOCATABLE :: ColAx,ColAz,ColDxx,ColDxz,ColDzz
REAL*8, DIMENSION(:,:), ALLOCATABLE :: CoefARK,CoefBRK
REAL*8, DIMENSION(:), ALLOCATABLE :: Fepzm,Fepzp
REAL*8, DIMENSION(:), ALLOCATABLE :: Fepxm,Fepxp
REAL*8, DIMENSION(:,:), ALLOCATABLE :: Fepsm,Fepsp
REAL*8, DIMENSION(:), ALLOCATABLE :: VRNeNs,VRTeTs,VRTfparTs, VRTfperpTs,VRTbparTs, VRTbperpTs,VRTiTs,VDRNeNs
REAL*8, DIMENSION(:,:), ALLOCATABLE :: DDzzduz
! Resonance conditions, and interpolation indexes:
! For the values of Qx,Qz in the second grid of the wave spectra:
INTEGER, DIMENSION(:), ALLOCATABLE :: IQresx2,IQresz2
INTEGER, DIMENSION(:), ALLOCATABLE :: IQresx,IQresz
! For the values of Qx,Qz in the third grid of the wave spectra (T wave scat.):
INTEGER, DIMENSION(:), ALLOCATABLE :: IQresx3,IQresz3
REAL*8, DIMENSION(:,:), ALLOCATABLE :: ILgrid3
! For the values of Qx,Qz in the wave spectra, scattering terms:
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: IQxLLL1,IQzLLL1
INTEGER, DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: IQxLLL2,IQzLLL2
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: IQxLLT1,IQzLLT1
INTEGER, DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: IQxLLT2i,IQzLLT2i
INTEGER, DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: IQxLLT2e,IQzLLT2e
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: IQxTLT1,IQzTLT1
INTEGER, DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: IQxTLT2i,IQzTLT2i
INTEGER, DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: IQxTLT2e,IQzTLT2e
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: IQQLLLs
INTEGER, DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: IQQLLLsi
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: IQQTLTs
INTEGER, DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: IQQTLTsi
INTEGER, DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: IQQTLTse
! For the coefficients Dij:
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: Qzr1D,Qzr2D
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: Qxr1D,Qxr2D
INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE :: IQzr1D,IQzr2D
INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE :: IQxr1D,IQxr2D
! L waves, spontaneous and induced emission:
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: UzrpLql,UzrmLql
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: UxrpLql,UxrmLql
INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE :: IUzrpLql,IUzrmLql
INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE :: IUxrpLql,IUxrmLql
! L waves, decay, L and S waves:
REAL*8, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: QxpLLSd
INTEGER, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: IQxpLLSd
INTEGER, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: IQxpdifLLSd
! L waves, decay, L and T waves:
REAL*8, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: QxpLLTd
INTEGER, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: IQxpLLTd
INTEGER, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: IQxpdifLLTd
! L waves, decay, S and T waves:
REAL*8, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: QxpLSTd
INTEGER, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: IQxpLSTd
INTEGER, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: IQxpdifLSTd
! L waves, decay, T and T waves:
REAL*8, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: QxpLTTd
INTEGER, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: IQxpLTTd
INTEGER, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: IQxpdifLTTd
! S waves, spontaneous and induced emission:
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: UzrpSql,UzrmSql
REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: UxrpSql,UxrmSql
INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE :: IUzrpSql,IUzrmSql
INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE :: IUxrpSql,IUxrmSql
! S waves, decay, L and L waves:
REAL*8, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: QxpSLLd
INTEGER, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: IQxpSLLd
INTEGER, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: IQxpdifSLLd
! S waves, decay, L and T waves:
REAL*8, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: QxpSLTd
INTEGER, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: IQxpSLTd
INTEGER, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: IQxpdifSLTd
! T waves, decay, L and L waves:
REAL*8, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: QxpTLLd
INTEGER, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: IQxpTLLd
INTEGER, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: IQxpdifTLLd
! T waves, decay, L and S waves:
REAL*8, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: QxpTLSd
INTEGER, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: IQxpTLSd
INTEGER, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: IQxpdifTLSd
! T waves, decay, T and T waves:
REAL*8, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: QxpTTLd
INTEGER, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: IQxpTTLd
INTEGER, DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: IQxpdifTTLd
! Auxiliary for resonance conditions:
REAL*8, DIMENSION(:), ALLOCATABLE :: Aux1_Rcd
INTEGER, DIMENSION(:), ALLOCATABLE :: Aux2_Rcd
! Auxiliary for collisional damping:
REAL*8, DIMENSION(:), ALLOCATABLE :: Aux1_Gcoll
INTEGER, DIMENSION(:), ALLOCATABLE :: Aux2_Gcoll
REAL*8, DIMENSION(:), ALLOCATABLE :: Aux1_Brem

END MODULE Common_Arrays

MODULE Math_Constants
IMPLICIT none

REAL*8 :: Pi= 3.1415926535897932384626433832795
REAL*8 :: Sqtwo= 1.414213562
REAL*8 :: Infinity= 1.0E+30
REAL*8 :: Degree= 0.017453293
REAL*8 :: EpsMin= 1.0E-16
REAL*8 :: Xacc= 1.0E-6
REAL*8 :: Qmin= 1.0E-4

COMPLEX*16:: Zi= (0.,1.)
COMPLEX*16:: Zzero= (0.,0.)

END MODULE Math_Constants

MODULE Phys_Constants
IMPLICIT none

REAL*8 :: MeC2= 511.01E+3          ! electron mass (eV)
REAL*8 :: MpC2= 938.50E+6          ! proton mass (eV)
REAL*8 :: C_SI= 2.997925E+8        ! speed of light (m/s)
REAL*8 :: C_cgs= 2.997925E+10      ! speed of light (cm/s)

END MODULE Phys_Constants

!---------------------------------------------------------------------
! Module with Functions and Subroutines:
!---------------------------------------------------------------------

MODULE Sub_Prog
CONTAINS

SUBROUTINE Allocate_Arrays
USE Common_Params
USE Common_Arrays
IMPLICIT NONE

ALLOCATE (VRz(nrz))
ALLOCATE (VQx(nqx))
ALLOCATE (VQz(nqz))
ALLOCATE (VQx2(nqx2))
ALLOCATE (VQz2(nqz2))
ALLOCATE (VQx3(nqx2))
ALLOCATE (VQz3(nqz2))
ALLOCATE (VQQ(nqmod))
ALLOCATE (ILp(nqx,nqz),ILm(nqx,nqz))
ALLOCATE (ILpOld(nqx, nqz), DeltaILp(nqx, nqz))
ALLOCATE (ILmOld(nqx, nqz), DeltaILm(nqx, nqz))
ALLOCATE (ISp(nqx,nqz),ISm(nqx,nqz))
ALLOCATE (ITp(nqx,nqz),ITm(nqx,nqz))
ALLOCATE (GcollLp(nqx,nqz),GqlLp(nqx,nqz))
ALLOCATE (GcollLm(nqx,nqz),GqlLm(nqx,nqz))
ALLOCATE (GcollL1D(nqmod),GqlL1D(nqmod))
ALLOCATE (GcollSp(nqx,nqz),GqlSp(nqx,nqz))
ALLOCATE (GcollSm(nqx,nqz),GqlSm(nqx,nqz))
ALLOCATE (GcollS1D(nqmod),GqlS1D(nqmod))
ALLOCATE (HGcollLp(nqx,nqz),HGcollSp(nqx,nqz))
ALLOCATE (HGcollL1D(nqmod),HGcollS1D(nqmod))
ALLOCATE (BremL(nqx,nqz),BremS(nqx,nqz))
ALLOCATE (BremL1D(nqmod),BremS1D(nqmod))
ALLOCATE (HBremL(nqx,nqz),HBremS(nqx,nqz))
ALLOCATE (HBremL1D(nqmod),HBremS1D(nqmod))
ALLOCATE (VPhip(nph))
ALLOCATE (VUx(nux))
ALLOCATE (VUz(nuz))
ALLOCATE (Fe0(nux,nuz),Fe(nux,nuz),Fi(nux,nuz),FeOld(nux,nuz),DeltaFe(nux,nuz))
ALLOCATE (Ax(nux,nuz),Az(nux,nuz))
ALLOCATE (Dxx(nux,nuz),Dxz(nux,nuz),Dzx(nux,nuz),Dzz(nux,nuz))
ALLOCATE (ColAx(nux,nuz),ColAz(nux,nuz))
ALLOCATE (CoefARK(nqx,nqz),CoefBRK(nqx,nqz))
ALLOCATE (ColDxx(nux,nuz),ColDxz(nux,nuz),ColDzz(nux,nuz))
ALLOCATE (Fepzm(nux),Fepzp(nux))
ALLOCATE (Fepxm(nuz),Fepxp(nuz))
ALLOCATE (Fepsm(nux,nuz),Fepsp(nux,nuz))
ALLOCATE (VRNeNs(nrz),VRTeTs(nrz),VRTfparTs(nrz), VRTfperpTs(nrz),VRTbparTs(nrz), VRTbperpTs(nrz))
ALLOCATE (VRTiTs(nrz),VDRNeNs(nrz))

ALLOCATE (DDzzduz(nux, nuz))
! Resonance conditions:
! For the coefficients Dij:
ALLOCATE (Qzr1D(nux,nuz,nqx,-1:1),Qzr2D(nux,nuz,nqx,-1:1))
ALLOCATE (Qxr1D(nux,nuz,nqz,-1:1),Qxr2D(nux,nuz,nqz,-1:1))
ALLOCATE (IQzr1D(nux,nuz,nqx,-1:1),IQzr2D(nux,nuz,nqx,-1:1))
ALLOCATE (IQxr1D(nux,nuz,nqz,-1:1),IQxr2D(nux,nuz,nqz,-1:1))
! Auxiliary for resonance conditions:
ALLOCATE (Aux1_Rcd(3))
ALLOCATE (Aux2_Rcd(4))
! Auxiliary for collisional damping:
ALLOCATE (Aux1_Gcoll(4))
ALLOCATE (Aux2_Gcoll(1))
ALLOCATE (Aux1_Brem(5))

END SUBROUTINE Allocate_Arrays

SUBROUTINE Definitions
USE Common_Params
USE Common_Arrays
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL*8 :: RMpMe,RMeMp
INTEGER :: m

m= 1   ! Auxiliary quantity, for the spatial profile.

RMpMe= MpC2/MeC2
RMeMp= MeC2/MpC2
RMiMe= RMpMe
AA= SQRT(1.+3./RTeTi)/SQRT(RMiMe)/SQRT(2.)
Geff= G/(2.*SQRT(2.)*(4.*Pi)**2)
SELECT CASE(Alphai)
 CASE(1)
  Uik2= (Kappai-1.)/Kappai /RTeTi/RMiMe * VRTeTs(m)
 CASE(0)
  Uik2= 1./RTeTi/RMiMe * VRTeTs(m)
 CASE DEFAULT
  OPEN(98,FILE='Warning_Definitions.wt')
  WRITE(98,*) ' Alphai= ',Alphai
  WRITE(98,*) ' Alphai must be 1 or 0 !!'
  CLOSE(98)
  STOP
END SELECT
SELECT CASE(Alphae)
 CASE(1)
  Uek2= (Kappae-1.)/Kappae * VRTeTs(m)
 CASE(0)
  Uek2= VRTeTs(m)
 CASE DEFAULT
  OPEN(98,FILE='Warning_Definitions.wt')
  WRITE(98,*) ' Alphae= ',Alphae
  WRITE(98,*) ' Alphae must be 1 or 0 !!'
  CLOSE(98)
  STOP
END SELECT
Uek= SQRT(Uek2)
Uik= SQRT(Uik2)
Ue2= VRTeTs(m)
Ui2= 1./RTeTi/RMiMe * VRTeTs(m)
Ue= SQRT(Ue2)
Ui= SQRT(Ui2)
Uf2= 1./Ve2C2
!Uf2= 1.E+6

IF (Kappae<=30.) THEN
 RatioGammaem051= GAMMA(Kappae+Alphae-0.5)/GAMMA(Kappae+Alphae-1.)
 RatioGammaem115= GAMMA(Kappae+Alphae-1.)/GAMMA(Kappae+Alphae-1.5)
ELSE
 RatioGammaem051= SQRT((Kappae+Alphae-1.5)/(Kappae+Alphae-2.)) &
  * (Kappae+Alphae)**0.5*EXP(-0.5) &
  * (1.-1.5/(Kappae+Alphae))**(Kappae+Alphae-1.5) &
  / (1.-2./(Kappae+Alphae))**(Kappae+Alphae-2.)
 RatioGammaem115= SQRT((Kappae+Alphae-2.)/(Kappae+Alphae-2.5)) &
  * (Kappae+Alphae)**(0.5)*EXP(-0.5) &
  * (1.-2./(Kappae+Alphae))**(Kappae+Alphae-2.) &
  /((1.-2.5/(Kappae+Alphae))**(Kappae+Alphae-2.5))
END IF
RatioGammaep051= RatioGammaem051*(Kappae+Alphae-0.5)
RatioGammae015= RatioGammaem115*(Kappae+Alphae-1.)
IF (Kappai<=30.) THEN
 RatioGammaim115= GAMMA(Kappai+Alphai-1.)/GAMMA(Kappai+Alphai-1.5)
ELSE
 RatioGammaim115= SQRT((Kappai+Alphai-2.)/(Kappai+Alphai-2.5)) &
  * (Kappai+Alphai)**(0.5)*EXP(-0.5) &
  * (1.-2./(Kappai+Alphai))**(Kappai+Alphai-2.) &
  /((1.-2.5/(Kappai+Alphai))**(Kappai+Alphai-2.5))
END IF
RatioGammai015= RatioGammaim115*(Kappai+Alphai-1.)

END SUBROUTINE Definitions

SUBROUTINE Space_Profiles
USE Common_Params
USE Common_Arrays
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL*8 :: Rz,Drz,R
REAL*8 :: Rzf,Rzi
INTEGER :: m

Rzf= 101.
Rzi= 100.

IF (nrz==1) THEN
 VRz(1)= Rzi
ELSE
 Drz= (Rzf-Rzi)/(nrz-1)
 DO m= 1,nrz
  VRz(m)= Rzi+(m-1)*Drz
 END DO
 VRz(nrz)= Rzf
END IF

DO m= 1,nrz
 Rz= VRz(m)
 R= Rz

! TODO
! Não entendi direito o que essa parte faz. Isso seta os vetores abaixo com os mesmos valores sempre?

 VRNeNs(m)= 1.
 VDRNeNs(m)= 0.
 VRTeTs(m)= 1.!RTsTe   ! Ratio of Tstar/Tparallel of background plasma
 !VRTeTs(m)= RTsTe ! Ratio of Tstar/Tperp     of background plasma

 VRTfparTs(m)= RTfparTs   ! Ratio of Tstar/Tparallel of forward beam
 VRTfperpTs(m)= RTfperpTs ! Ratio of Tstar/Tperp     of forward beam

 VRTbparTs(m)= RTbparTs  ! Ratio of Tstar/Tparallel of backward beam
 VRTbperpTs(m)= RTbperpTs ! Ratio of Tstar/Tperp     of backward beam

 VRTiTs(m)= 1./RTeTi
END DO

RETURN
END SUBROUTINE Space_Profiles


!---------------------------------------------------------------------
! Functions and Subroutines:
!---------------------------------------------------------------------

SUBROUTINE Init_Wave(Tau,Dqx,Dqz,Dux,Duz,Anorm,Epart,Ewave,EppEw,&
   EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3)
USE Common_Params
USE Common_Arrays
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL*8 :: Fesum,Femax,Febkgr,Fef,Feb
REAL*8 :: Iinit
REAL*8, INTENT(out) :: Tau
REAL*8, INTENT(out) :: Dqx,Dqz,Dux,Duz
REAL*8 :: Dqx2,Dqz2,Dq
REAL*8 :: DPhip,Phip
REAL*8 :: Ux,Uz,U,U2,Qx,Qz 
REAL*8 :: HBrL,HBrS,HGcL,HGcS
REAL*8, INTENT(out) :: Anorm,Epart,Ewave,EppEw
REAL*8, INTENT(out) :: EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3
REAL*8 :: Aux
REAL*8, DIMENSION(nux,nuz) :: Dfdux,Dfduz
INTEGER :: i,j,k,nu2
CHARACTER(LEN=1) :: WaveT

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
Dqx2= (Qxf-0.1*Qmin-Qxi-0.1*Qmin)/(nqx2-1)
DO i= 1,nqx2
 VQx2(i)= Qxi+0.1*Qmin+(i-1)*Dqx2
END DO
Dqz2= (Qzf-0.1*Qmin-Qzi-0.1*Qmin)/(nqz2-1)
DO i= 1,nqz2
 VQz2(i)= Qzi+0.1*Qmin+(i-1)*Dqz2
END DO
Aux= SQRT(2./3./Ve2C2)
Dqx2= (Aux*Qxf-0.1*Qmin-Qxi-0.1*Qmin)/(nqx2-1)
DO i= 1,nqx2
 VQx3(i)= Qxi+0.1*Qmin+(i-1)*Dqx2
END DO
VQx3(nqx2)= Aux*Qxf
Dqz2= (Aux*Qzf-0.1*Qmin-Qzi-0.1*Qmin)/(nqz2-1)
DO i= 1,nqz2
 VQz3(i)= Qzi+0.1*Qmin+(i-1)*Dqz2
END DO
VQz3(nqz2)= Aux*Qzf
! One-dimensional array (absolute value of Q):
Dq= (Qf-Qi)/(nqmod-1)
DO i= 1,nqmod
 VQQ(i)= Qi+(i-1)*Dq
END DO
VQQ(nqmod)= Qf

! Initializes an auxiliary array, to be used in the scattering terms:
DPhip= (Pi/2.-0.)/(nph-1)
DO i= 1,nph
 Phip= 0.+(i-1)*DPhip
 VPhip(i)= Phip
END DO
VPhip(nph)= Pi/2.

! Initializes the collisional damping and bremsstrahlung coefficients, assumed
! to be constant along time evolution:
CALL Coll_Damping
CALL Bremsstrahlung

! Initializes the spectrum of L waves and S waves:
SELECT CASE(InitialLevel)

 CASE("Choose")
  ILp= Iw0
  ILm= Iw0
  ISp= Iw0
  ISm= Iw0
  ITp= Iw0
  ITm= Iw0
  AuxInitialLevel= 0.E0
  DO i= 1,nqx
   DO j= 1,nqz
    ILp(i,j)= Iw0
    ILm(i,j)= Iw0
    ISp(i,j)= Iw0
    ISm(i,j)= Iw0
    ITp(i,j)= Iw0
    ITm(i,j)= Iw0
   END DO
  END DO

 CASE("Auto  ")
  AuxInitialLevel= 1.E0
  DO i= 1,nqx
   Qx= VQx(i)
   DO j= 1,nqz
    Qz= VQz(j)
    HBrL= HBremL(i,j)
    HBrS= HBremS(i,j)
    HGcL= HGcollLp(i,j)
    HGcS= HGcollSp(i,j)
    WaveT= "L"
    CALL Iwave_Init(Qx,Qz,Iinit,WaveT,HBrL,HBrS,HGcL,HGcS)
    ILp(i,j)= Iinit
    ILm(i,j)= Iinit
    WaveT= "S"
    CALL Iwave_Init(Qx,Qz,Iinit,WaveT,HBrL,HBrS,HGcL,HGcS)
    ISp(i,j)= Iinit
    ISm(i,j)= Iinit
    WaveT= "T"
    CALL Iwave_Init(Qx,Qz,Iinit,WaveT,HBrL,HBrS,HGcL,HGcS)
    ITp(i,j)= Iinit
    ITm(i,j)= Iinit
   END DO
  END DO

 CASE("Null  ")
  Iw0= 0.E0
  ILp= Iw0
  ILm= Iw0
  ISp= Iw0
  ISm= Iw0
  ITp= Iw0
  ITm= Iw0
  AuxInitialLevel= 0.E0
  DO i= 1,nqx
   DO j= 1,nqz
    ILp(i,j)= Iw0
    ILm(i,j)= Iw0
    ISp(i,j)= Iw0
    ISm(i,j)= Iw0
    ITp(i,j)= Iw0
    ITm(i,j)= Iw0
   END DO
  END DO

 CASE DEFAULT
  OPEN(98,FILE='Warning_Init_Wave.wt')
  WRITE(98,*) ' InitialLevel= ',InitialLevel
  WRITE(98,*) ' InitialLevel must be (Choose) or (Auto  ) or (Null  ) !!'
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
DO k= 1,nuz
 Uz= VUz(k)
 DO i= 1,nux
  Ux= VUx(i)
  U2= Ux**2+Uz**2
  U= SQRT(U2)
  CALL Fe_Init(Ux,Uz,Fesum,Femax,Febkgr,Fef,Feb)
  Fe(i,k)= Fesum
  Fi(i,k)= (RMiMe*RTeTi)*EXP(-RMiMe*RTeTi*U2)/(Pi)
 END DO
END DO
CALL Fnorm(Anorm)
CALL Energy(Epart,Ewave,EppEw)
CALL Energy2(EwaveL,EWaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3)
Fe0= Fe

! Initializes the boundary conditions for the eletron distribution: 
SELECT CASE(DerivLn)
 CASE("Yes")
  CALL Derivxy5pln2d(nux,nuz,VUx,VUz,Fe0,Dfdux,Dfduz)
 CASE("No ")
  CALL Derivxy5p2d(nux,nuz,VUx,VUz,Fe0,Dfdux,Dfduz)
END SELECT
DO k= 1,nuz
 Fepxm(k)= Dfdux(2,k)
 Fepxp(k)= Dfdux(nux-1,k)
END DO
DO i= 1,nux
 Fepzp(i)= Dfduz(i,nuz-1)
 Fepzm(i)= Dfduz(i,2)
END DO

RETURN
END SUBROUTINE Init_Wave

SUBROUTINE Iwave_Init(Qx,Qz,Iinit,WaveT,HBrL,HBrS,HGcL,HGcS)
USE Common_Params
USE Common_Arrays
USE Math_Constants
!USE hyp_2f1_module
IMPLICIT NONE
REAL*8, INTENT(in) :: Qx,Qz
REAL*8 :: Q,Qm
REAL*8 :: Zlq,Zsq,Zsoq,Ztq,Aux!,Muq
REAL*8, INTENT(out) :: Iinit
REAL*8 :: ILinit,ISinit,ITinit
REAL*8, INTENT(in) :: HBrL,HBrS,HGcL,HGcS
REAL*8 Aux1,Aux2,Aux3,Aux4
REAL*8 :: Aux1i,Aux1e,Aux2i,Aux2e
REAL*8 :: Auxnum,Auxden
REAL*8 :: Auxnume,Auxdene
REAL*8 :: Auxnumi,Auxdeni
INTEGER :: m
CHARACTER(LEN=1), INTENT(in) :: WaveT

m= 1   ! Auxiliary quantity, for the spatial profile.

Q= SQRT(Qx**2+Qz**2)
IF (Q>=Qmin) THEN
 Q= Q
ELSE
 Q= Qmin
END IF

SELECT CASE(WaveT)
 CASE("L")
  !IF (AuxInitialLevel==0.) THEN
  IF (ABS(AuxInitialLevel).lt.EpsMin) THEN
   ILinit= Iw0
  ELSE
   IF (Lsaturation=="Yes" .AND. Q<1./SQRT(1./Ve2C2-1.5)) THEN
    Q= 1./SQRT(1./Ve2C2-1.5)
   ELSE
    Q= Q
   END IF
   Zlq= ZL(Q,0.D0)
   ! Contribution from the Maxwellian background:
   Aux= Geff*VRNeNs(m)*VRTeTs(m)/2./(Zlq**2)
   Auxnum= (1.-RatioNke)*EXP(-Zlq**2/Q**2/VRTeTs(m))
   Auxden= Auxnum
   ! Contribution from the halo kappa distribution:
   Aux1= Kappae*Uek2
   Aux2= 1.+(Zlq/Q)**2/Aux1
   Aux3= RatioNke*Ue/SQRT(Aux1)/Aux2**(Kappae+Alphae-1.)
   Auxnum= Auxnum + Aux3 * RatioGammaem115
   Auxden= Auxden + Aux3 * Ue**2/Aux1/Aux2 * RatioGammae015
  SELECT CASE(NewEffects1)
   CASE("Yes")
    ILinit= Aux*(Auxnum + HBrL)/(Auxden - 2.*HGcL)
   CASE("No ")
    IF (RatioNke<=1.E-16) THEN
     ILinit= Aux
    ELSE
     IF (Zlq**2/Q**2/VRTeTs(m)>100.) THEN ! Maxwellian contribution vanishes
      ILinit= Aux*Uek2/Ue2*Kappae/(Kappae+Alphae-1.)*Aux2
     ELSE
      ILinit= Aux*(Auxnum)/(Auxden)
     END IF
    END IF
   CASE DEFAULT
  END SELECT 
  END IF
  Iinit= ILinit
 CASE("S")
  !IF (AuxInitialLevel==0.) THEN
  IF (ABS(AuxInitialLevel).lt.EpsMin) THEN
   ISinit= Iw0
  ELSE
   Zlq= ZL(Q,0.D0)
   Zsq= ZS(Q,0.D0)
   Zsoq= AA*SQRT(VRTeTs(m))/SQRT(1.+Q**2/2.*VRTeTs(m)/VRNeNs(m))
   !Muq= Q**3*AA/2.
   ! Contribution from the Maxwellian background:
   Auxnume= (1.-RatioNke)*EXP(-Zsoq**2/VRTeTs(m))
   Auxnumi= (1.-RatioNki)*Ue/Ui*EXP(-Zsoq**2/Ui2)
   Auxdene= Auxnume
   Auxdeni= Auxnumi * RTeTi
   ! Contribution from the halo kappa distribution:
   Aux1e= Kappae*Uek2
   Aux2e= 1.+Zsoq**2/Aux1e
   Aux1i= Kappai*Uik2
   Aux2i= 1.+Zsoq**2/Aux1i
   Aux= VRNeNs(m)*VRTeTs(m)*Geff/2./Zlq/Zsq
   Aux3= RatioNke*Ue/SQRT(Aux1e)/Aux2e**(Kappae+Alphae-1.)
   Aux4= RatioNki*Ui/SQRT(Aux1i)/Aux2i**(Kappai+Alphai-1.)
   Auxnume= Auxnume+ Aux3 * RatioGammaem115
   Auxnumi= Auxnumi+ Aux4 * RatioGammaim115
   Auxdene= Auxdene+ Aux3 * Ue**2/Aux1e/Aux2e * RatioGammae015
   Auxdeni= Auxdeni+ Aux4 * Ui**2/Aux1i/Aux2i * RatioGammai015
  SELECT CASE(NewEffects1)
   CASE("Yes")
    ISinit= Aux*(Auxnume+Auxnumi+HBrS)/(Auxdene+Auxdeni-2.*HGcS)
   CASE("No ")
    ISinit= Aux*(Auxnume+Auxnumi)/(Auxdene+Auxdeni)
   CASE DEFAULT
  END SELECT 
  END IF
  Iinit= ISinit
 CASE("T")
  !IF (AuxInitialLevel==0.) THEN
  IF (ABS(AuxInitialLevel).lt.EpsMin) THEN
   ITinit= EpsMin
  ELSE
   IF (AsympT==1) THEN
    Qm= Q*SQRT(2./3./Ve2C2)
    Ztq= ZL(Qm,0.D0)
    ! Contribution from the Maxwellian background:
    Aux= Geff*VRNeNs(m)/VRTeTs(m)/(Ztq**2)
    Auxnum= (1.-RatioNke)*EXP(-Ztq**2/Qm**2/VRTeTs(m))
    Auxden= Auxnum
    ! Contribution from the halo kappa distribution:
    Aux1= Kappae*Uek2
    Aux2= 1.+(Ztq/Qm)**2/Aux1
    Aux3= RatioNke*Ue/SQRT(Aux1)/Aux2**(Kappae+Alphae-1.)
    Auxnum= Auxnum + Aux3 * RatioGammaem115
    Auxden= Auxden + Aux3 * Ue**2/Aux1/Aux2 * RatioGammae015
    SELECT CASE(NewEffects1)
     CASE("Yes")
      OPEN(98,FILE='Warning_Bremsstrahlung.wt')
      WRITE(98,*) ' NewEffects1= ',NewEffects1,'  AsympT=',AsympT
      WRITE(98,*) 'New effects are not implemented for initial T wave spectrum'
      CLOSE(98)
      STOP
     CASE("No ")
      IF (RatioNke<=1.E-16) THEN
       ITinit= Aux
      ELSE
       IF (Ztq**2/Qm**2/VRTeTs(m)>100.) THEN ! Maxwellian contr. vanishes
        ITinit= Aux*Uek2/Ue2*Kappae/(Kappae+Alphae-1.)*Aux2
       ELSE
        ITinit= Aux*(Auxnum)/(Auxden)
       END IF
      END IF
     CASE DEFAULT
    END SELECT 
   ELSE
    ITinit= 0.
   END IF
  END IF
  Iinit= ITinit
 CASE DEFAULT
END SELECT
RETURN
END SUBROUTINE Iwave_Init

SUBROUTINE Fe_Init(Ux,Uz,Fesum,Femax,Febkgr,Fef,Feb)
USE Common_Params
USE Common_Arrays
USE Math_Constants
IMPLICIT NONE
REAL*8, INTENT(in) :: Ux,Uz
REAL*8, INTENT(out) :: Fesum,Femax,Febkgr,Fef,Feb
!REAL*8 :: Kappa
INTEGER :: m

m= 1   ! Auxiliary quantity, for the spatial profile.

Ufx= Uf*SIN(Thf*Pi/180.) ! ufperp
Ufz= Uf*COS(Thf*Pi/180.) ! ufpar
Ubx= Ub*SIN(Thb*Pi/180.) ! ubperp
Ubz= Ub*COS(Thb*Pi/180.) ! ubpar

U0x= 0. ! u0perp

U0z= -(RatioNf*Ufz+RatioNb*Ubz)/(1.-RatioNf-RatioNb) ! uopar

Femax= (1.-RatioNf-RatioNb-RatioNke)/(Pi*VRTeTs(m)) & ! background plasma maxwellian
  *EXP(-((Ux-U0x)**2+(Uz-U0z)**2)/VRTeTs(m)) ! background plasma maxwellian

Febkgr= RatioNke*((Kappae+Alphae-1.)/Kappae)/(Pi*Uek2) & ! Dist. Kappa: n mudar essa
  *(1.+((Ux-U0x)**2+(Uz-U0z)**2)/Kappae/Uek2)**(-(Kappae+Alphae)) ! Dist. Kappa: n mudar essa

Fef= (RatioNf/(Pi)) / (SQRT(VRTfperpTs(m)*VRTfparTs(m))) * (1/(1.+PERF(Ufx/SQRT(VRTfperpTs(m))))) & ! forward beam
   * EXP(-(((Ux-Ufx)**2))/VRTfperpTs(m) -(((Uz-Ufz)**2))/VRTfparTs(m)) ! forward beam

Feb= (RatioNb/(Pi)) / (SQRT(VRTbperpTs(m)*VRTbparTs(m))) * (1/(1.+PERF(Ubx/SQRT(VRTbperpTs(m))))) & ! backward beam
   * EXP(-(((Ux-Ubx)**2))/VRTbperpTs(m) - (((Uz-Ubz)**2))/VRTbparTs(m)) ! backward beam

Fesum= Femax+Febkgr+Fef+Feb

RETURN

END SUBROUTINE Fe_Init

SUBROUTINE Save_Results(Tau,Dqx,Dqz,Dux,Duz,Anorm0,Epart0,Ewave0,EppEw0,&
        Rn,Rp,Rw,Rs,EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3)
USE Common_Params
USE Common_Arrays
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL*8, INTENT(in) :: Tau 
REAL*8, INTENT(in) :: Dqx,Dqz,Dux,Duz
REAL*8, INTENT(in) :: Anorm0,Ewave0,Epart0,EppEw0
REAL*8, INTENT(in) :: Rn,Rp,Rw,Rs
REAL*8, INTENT(in) :: EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3
INTEGER :: i,k,iflag

Iflag= 1
OPEN(1,FILE='Out.wt')
WRITE(1,*) Iflag
WRITE(1,*) nux,nuz
WRITE(1,*) nqx,nqz
WRITE(1,*) nqx2,nqz2
WRITE(1,*) nqmod
WRITE(1,*) nrz
WRITE(1,*) Ulim
WRITE(1,*) Ucrit
WRITE(1,*) Qxi 
WRITE(1,*) Qxf 
WRITE(1,*) Qzi 
WRITE(1,*) Qzf 
WRITE(1,*) Qi
WRITE(1,*) Qf 
WRITE(1,*) InitialLevel
WRITE(1,*) Iw0 
WRITE(1,*) AsympT 
WRITE(1,*) RatioNke 
WRITE(1,*) RatioNki 
WRITE(1,*) Kappae,Alphae
WRITE(1,*) Kappai,Alphai
WRITE(1,*) RatioNf 
WRITE(1,*) Uf 
WRITE(1,*) Thf 
WRITE(1,*) RTfparTs 
WRITE(1,*) RTfperpTs 
WRITE(1,*) RatioNb 
WRITE(1,*) Ub 
WRITE(1,*) Thb 
WRITE(1,*) RTbparTs 
WRITE(1,*) RTbperpTs 
WRITE(1,*) RTeTi 
WRITE(1,*) G
WRITE(1,*) Ve2C2 
WRITE(1,*) Lemis
WRITE(1,*) LdecayLS
WRITE(1,*) LdecayLT
WRITE(1,*) LdecayST
WRITE(1,*) LdecayTT
WRITE(1,*) LscatLL
WRITE(1,*) LscatLT
WRITE(1,*) Semis
WRITE(1,*) SdecayLL
WRITE(1,*) SdecayLT
WRITE(1,*) Sscat
WRITE(1,*) TdecayLL
WRITE(1,*) TdecayLS
WRITE(1,*) TdecayTL
WRITE(1,*) TscatLT
WRITE(1,*) CollTerm
WRITE(1,*) CollTermForm
WRITE(1,*) SpontEmis
WRITE(1,*) ScatElSpo
!WRITE(1,*) ScatElInd
WRITE(1,*) NewEffects1
WRITE(1,*) NewEffects2
WRITE(1,*) Gcoll
WRITE(1,*) Bremss
WRITE(1,*) RenormFe
WRITE(1,*) DerivLn
WRITE(1,*) RebuildL
WRITE(1,*) RebuildS
WRITE(1,*) RebuildFe
WRITE(1,*) Lsaturation
! Parameters generated in the code:
WRITE(1,*) RMiMe
WRITE(1,*) AA,Geff 
WRITE(1,*) U0x,U0z
WRITE(1,*) AuxInitialLevel
WRITE(1,*) Auxinterp

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
DO i= 1,nqx2
 WRITE(1,*) VQx2(i)
END DO
DO k= 1,nqz2
 WRITE(1,*) VQz2(k)
END DO
DO i= 1,nqx2
 WRITE(1,*) VQx3(i)
END DO
DO k= 1,nqz2
 WRITE(1,*) VQz3(k)
END DO
DO i= 1,nph
 WRITE(1,*) VPhip(i)
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
DO i= 1,nqx
 DO k= 1,nqz
  WRITE(1,*) ITp(i,k)
 END DO
END DO
DO i= 1,nqx
 DO k= 1,nqz
  WRITE(1,*) ITm(i,k)
 END DO
END DO
DO i= 1,nux
 WRITE(1,*) Fepzp(i),Fepzm(i)
END DO
DO k= 1,nuz
 WRITE(1,*) Fepxp(k),Fepxm(k)
END DO
WRITE(1,*) nux
WRITE(1,*) nuz
WRITE(1,*) nqx
WRITE(1,*) nqz
WRITE(1,*) nqx2
WRITE(1,*) nqz2
WRITE(1,*) nqmod
WRITE(1,*) nrz
WRITE(1,*) Ulim
WRITE(1,*) Qxi
WRITE(1,*) Qxf
WRITE(1,*) Qzi
WRITE(1,*) Qzf
WRITE(1,*) Qi
WRITE(1,*) Qf
WRITE(1,*) Dqx
WRITE(1,*) Dqz
WRITE(1,*) Dux
WRITE(1,*) Duz
WRITE(1,*) U0x
WRITE(1,*) U0z
WRITE(1,*) Anorm0
WRITE(1,*) Epart0
WRITE(1,*) Ewave0
WRITE(1,*) EppEw0
WRITE(1,*) Rn
WRITE(1,*) Rp
WRITE(1,*) Rw
WRITE(1,*) Rs
WRITE(1,*) EwaveL
WRITE(1,*) EwaveS
WRITE(1,*) EwaveT
WRITE(1,*) EwaveTF
WRITE(1,*) EwaveTH
WRITE(1,*) EwaveT3
WRITE(1,*) Tau
CLOSE(1)

RETURN
END SUBROUTINE Save_Results 

SUBROUTINE Read_Results(Tau,Dqx,Dqz,Dux,Duz,Anorm0,Epart0,Ewave0,EppEw0,&
     Rn,Rp,Rw,Rs,EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3)
USE Common_Params
USE Common_Arrays
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL*8, INTENT(out) :: Tau
REAL*8, INTENT(out) :: Dqx,Dqz,Dux,Duz
REAL*8, INTENT(out) :: Anorm0,Ewave0,Epart0,EppEw0
REAL*8, INTENT(out) :: Rn,Rp,Rw,Rs
REAL*8, INTENT(out) :: EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3
REAL*8 :: Ulimback,Qxiback,Qxfback,Qziback,Qzfback,Qiback,Qfback
INTEGER :: nuxback,nuzback,nqxback,nqzback,nqx2back,nqz2back,nrzback,nqmodback
INTEGER :: i,k

! Parameters generated in the code:
READ(1,*) RMiMe
READ(1,*) AA,Geff 
READ(1,*) U0x,U0z
READ(1,*) AuxInitialLevel
READ(1,*) Auxinterp
!
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
DO i= 1,nqx2
 READ(1,*) VQx2(i)
END DO
DO k= 1,nqz2
 READ(1,*) VQz2(k)
END DO
DO i= 1,nqx2
 READ(1,*) VQx3(i)
END DO
DO k= 1,nqz2
 READ(1,*) VQz3(k)
END DO
DO i= 1,nph
 READ(1,*) VPhip(i)
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
DO i= 1,nqx
 DO k= 1,nqz
  READ(1,*) ITp(i,k)
 END DO
END DO
DO i= 1,nqx
 DO k= 1,nqz
  READ(1,*) ITm(i,k)
 END DO
END DO
DO i= 1,nux
 READ(1,*) Fepzp(i),Fepzm(i)
END DO
DO k= 1,nuz
 READ(1,*) Fepxp(k),Fepxm(k)
END DO
READ(1,*) nuxback
READ(1,*) nuzback
READ(1,*) nqxback
READ(1,*) nqzback
READ(1,*) nqx2back
READ(1,*) nqz2back
READ(1,*) nqmodback
READ(1,*) nrzback
READ(1,*) Ulimback
READ(1,*) Qxiback
READ(1,*) Qxfback
READ(1,*) Qziback
READ(1,*) Qzfback
READ(1,*) Qiback
READ(1,*) Qfback
READ(1,*) Dqx
READ(1,*) Dqz
READ(1,*) Dux
READ(1,*) Duz
READ(1,*) U0x
READ(1,*) U0z
READ(1,*) Anorm0
READ(1,*) Epart0
READ(1,*) Ewave0
READ(1,*) EppEw0
READ(1,*) Rn
READ(1,*) Rp
READ(1,*) Rw
READ(1,*) Rs
READ(1,*) EwaveL
READ(1,*) EwaveS
READ(1,*) EwaveT
READ(1,*) EwaveTF
READ(1,*) EwaveTH
READ(1,*) EwaveT3
READ(1,*) Tau
CLOSE(1)

IF(nuxback/=nux .OR. nuzback/=nuz .OR. nqxback/=nqx .OR. nqzback/=nqz &
        .OR. nrzback/=nrz .OR. nqx2back/=nqx2 .OR. nqz2back/=nqz2 &
        .OR. nqmodback/=nqmod) THEN
 OPEN(98,FILE='Warning_Read_Results.wt')
 WRITE(98,*) ' nuxback= ',nuxback,' nux= ',nux
 WRITE(98,*) ' nuzback= ',nuzback,' nuz= ',nuz
 WRITE(98,*) ' nqxback= ',nqxback,' nqx= ',nqx
 WRITE(98,*) ' nqzback= ',nqzback,' nqz= ',nqz
 WRITE(98,*) ' nqx2back= ',nqx2back,' nqx2= ',nqx2
 WRITE(98,*) ' nqz2back= ',nqz2back,' nqz2= ',nqz2
 WRITE(98,*) ' nrzback= ',nrzback,' nrz= ',nrz
 WRITE(98,*) ' nqmodback= ',nqmodback,' nqmod= ',nqmod
 WRITE(98,*) ' The quantities nux, nuz, nqx, nqz, nqx2, nqz2, nrz, nqmod, &
   & can not be modified !'
 CLOSE(98)
 STOP
ELSE
END IF

IF(Ulimback/=Ulim .OR. Qxiback/=Qxi .OR. Qxfback/=Qxf .OR. Qziback/=Qzi &
        .OR. Qzfback/=Qzf .OR. Qiback/=Qi .OR. Qfback/=Qf) THEN
 OPEN(98,FILE='Warning_Read_Results.wt')
 WRITE(98,*) ' Ulimback= ',Ulimback,' Ulim= ',Ulim
 WRITE(98,*) ' Qxiback= ',Qxiback,' Qxi= ',Qxi
 WRITE(98,*) ' Qxfback= ',Qxfback,' Qxf= ',Qxf
 WRITE(98,*) ' Qziback= ',Qziback,' Qzi= ',Qzi
 WRITE(98,*) ' Qzfback= ',Qzfback,' Qzf= ',Qzf
 WRITE(98,*) ' Qiback= ',Qiback,' Qi= ',Qi
 WRITE(98,*) ' Qfback= ',Qfback,' Qf= ',Qf
 WRITE(98,*) ' The quantities Ulim, Qxi, Qxf, Qzi, Qzf, Qi, Qf, &
   & can not be modified !'
 CLOSE(98)
 STOP
ELSE
END IF

RETURN
END SUBROUTINE Read_Results 

SUBROUTINE Res_Cond
USE Common_Params
USE Common_Arrays
USE Math_Constants
USE Phys_Constants
!USE hyp_2f1_module
IMPLICIT NONE
REAL*8 :: Ux,Uz,Qx,Qz
REAL*8 :: AbsUz,AbsUx
REAL*8 :: Q1,Q2,Qxp,Qzp
REAL*8 :: Qx2,Qz2,Q,Uresp,Uresm
REAL*8 :: Qxdif
REAL*8 :: Qp,Qp2
REAL*8 :: Zlq,Zsq,Ztq
REAL*8 :: Muq,Aux
REAL*8 :: X1,X2
REAL*8 :: Qstar,Qs,Phi,Phip,Aux1,Aalpha
REAL*8 :: AuxNum,CPhip,SPhip
REAL*8 :: Iinit
REAL*8 :: HBrL,HBrS,HGcL,HGcS
INTEGER :: i,l,m,sigma,sigmap,S1,S2,Asp,Aspp
INTEGER :: iqx,kqz,kqzp
INTEGER :: ires
INTEGER :: iresq,kresq
INTEGER :: j
INTEGER :: SS
LOGICAL :: Succes 
CHARACTER(LEN=1) :: WaveT

j= 1   ! Auxiliary quantity, for the spatial profile.

Qzr1D= 0.
Qzr2D= 0.
Qxr1D= 0.
Qxr2D= 0.
IQzr1D= 0.
IQzr2D= 0.
IQxr1D= 0.
IQxr2D= 0.

! For the values of Qx,Qz in the second grid of the wave spectra:
ALLOCATE (IQresx2(nqx2),IQresz2(nqz2))
! The points of the grid (nqx2,nqz2) are "inside" the grid (nqx,nqz), and
! therefore are correctly interpolated.
DO iqx= 1,nqx2
 Qx= VQx2(iqx)
 CALL Locate(VQx,nqx,Qx,iresq)
 IQresx2(iqx)= iresq
END DO
DO kqz= 1,nqz2
 Qz= VQz2(kqz)
 CALL Locate(VQz,nqz,Qz,kresq)
 IQresz2(kqz)= kresq
END DO
! For the values of Qx,Qz in the first grid of the wave spectra:
ALLOCATE (IQresx(nqx),IQresz(nqz))
! The extreme points of the grid (nqx,nqz) are "outside" of the grid
! (nqx2,nqz2), and therefore need a correction for correct interpolation.
! This is made with routine "Aitp2d2b"
DO iqx= 1,nqx
 Qx= VQx(iqx)
 CALL Locate(VQx2,nqx2,Qx,iresq)
 IQresx(iqx)= iresq
END DO
DO kqz= 1,nqz
 Qz= VQz(kqz)
 CALL Locate(VQz2,nqz2,Qz,kresq)
 IQresz(kqz)= kresq
END DO

! For the values of Qx,Qz in the third grid of the wave spectra (to be used
! in the scattering term in the equation for T waves):
ALLOCATE (IQresx3(nqx2),IQresz3(nqz2))
ALLOCATE (ILgrid3(nqx2,nqz2))
! The first points in the grid (small values of Qx and Qz) are "inside" the 
! grid (nqx,nqz), and therefore are correctly interpolated.
DO iqx= 1,nqx2
 Qx= VQx3(iqx)
 CALL Locate(VQx,nqx,Qx,iresq)
 IQresx3(iqx)= iresq
END DO
DO kqz= 1,nqz2
 Qz= VQz3(kqz)
 CALL Locate(VQz,nqz,Qz,kresq)
 IQresz3(kqz)= kresq
END DO
WaveT= "L"
DO iqx= 1,nqx2
 Qx= VQx3(iqx)
 DO kqz= 1,nqz2
  Qz= VQz3(kqz)
  Q2= Qx**2+Qz**2
  Q= SQRT(Q2)
  !IF(Q<=5.E-3) Q=5.E-3
  HBrS= 0.
  HGcS= 0.
  SELECT CASE(NewEffects1)
   CASE("Yes")
    CALL Locate(VQQ,nqmod,Q,ires)
    CALL Aitp1d2(nqmod,VQQ,HBremL1D,Q,Aux,ires)
    HBrL= Aux
    CALL Aitp1d2(nqmod,VQQ,HGcollL1D,Q,Aux,ires)
    HGcL= Aux
   CASE("No ")
    HBrL= 0.
    HGcL= 0.
   CASE DEFAULT
  END SELECT
  CALL Iwave_Init(Qx,Qz,Iinit,WaveT,HBrL,HBrS,HGcL,HGcS)
  ILgrid3(iqx,kqz)= Iinit
 END DO
END DO

! For the values of Qx,Qz in the wave spectra, scattering terms:
! L waves, scattering involving L and L waves
ALLOCATE (IQxLLL1(nqx,nqz,nph),IQzLLL1(nqx,nqz,nph))
ALLOCATE (IQxLLL2(nqx,nqz,nph,-1:1,-1:1,-1:1),&
          IQzLLL2(nqx,nqz,nph,-1:1,-1:1,-1:1))
DO iqx= 1,nqx
 Qx= VQx(iqx)
 DO kqz= 1,nqz
  Qz= VQz(kqz)
  Qx2= Qx**2
  Qz2= Qz**2
  Q2= Qx2+Qz2
  Q= SQRT(Q2)
  Phi= ACOS(Qz/Q)
  Zlq= ZL(Qx,Qz)
  Qstar= Q
  DO i= 1,nph
   Phip= VPhip(i)
   Qxp= Qstar*SIN(Phip)
   Qzp= Qstar*COS(Phip)
   CALL Locate(VQx,nqx,Qxp,ires)
   IQxLLL1(iqx,kqz,i)= ires
   CALL Locate(VQz,nqz,Qzp,ires)
   IQzLLL1(iqx,kqz,i)= ires
   DO S2= -1,1,2
    DO S1= -1,1,2
     Aux1= (Q2+Qstar**2-2.*S1*Q*Qstar*COS(Phi-S1*S2*Phip))
     IF (Aux1<1.E-6) THEN
      Aalpha= RMiMe/VRTiTs(j)*(9./4.)*(Qstar/Zlq)**2/1.E-6 
     ELSE
      Aalpha= RMiMe/VRTiTs(j)*(9./4.)*(Qstar/Zlq)**2/Aux1 
     END IF
     DO SS= -1,1,2
      Qs= Qstar+SS/SQRT(2.*Aalpha)
      Qxp= Qs*SIN(Phip)
      Qzp= Qs*COS(Phip)
      CALL Locate(VQx,nqx,Qxp,ires)
      IQxLLL2(iqx,kqz,i,S2,S1,SS)= ires
      CALL Locate(VQz,nqz,Qzp,ires)
      IQzLLL2(iqx,kqz,i,S2,S1,SS)= ires
     END DO
    END DO
   END DO
  END DO
 END DO
END DO

! L waves, scattering involving L and T waves
ALLOCATE (IQxLLT1(nqx,nqz,nph),IQzLLT1(nqx,nqz,nph))
ALLOCATE (IQxLLT2i(nqx,nqz,nph,-1:1,-1:1,-1:1),&
          IQzLLT2i(nqx,nqz,nph,-1:1,-1:1,-1:1))
ALLOCATE (IQxLLT2e(nqx,nqz,nph,-1:1,-1:1,-1:1),&
          IQzLLT2e(nqx,nqz,nph,-1:1,-1:1,-1:1))
DO iqx= 1,nqx
 Qx= VQx(iqx)
 DO kqz= 1,nqz
  Qz= VQz(kqz)
  Qx2= Qx**2
  Qz2= Qz**2
  Q2= Qx2+Qz2
  Q= SQRT(Q2)
  Phi= ACOS(Qz/Q)
  Zlq= ZL(Qx,Qz)
  Qstar= SQRT(3./2.*Ve2C2)*Q
  DO i= 1,nph
   Phip= VPhip(i)
   Qxp= Qstar*SIN(Phip)
   Qzp= Qstar*COS(Phip)
   CALL Locate(VQx,nqx,Qxp,ires)
   IQxLLT1(iqx,kqz,i)= ires
   CALL Locate(VQz,nqz,Qzp,ires)
   IQzLLT1(iqx,kqz,i)= ires
   DO S2= -1,1,2
    DO S1= -1,1,2
     Aux1= (Q2+Qstar**2-2.*S1*Q*Qstar*COS(Phi-S1*S2*Phip))
     IF (Aux1<1.E-6) THEN
      Aalpha= RMiMe/VRTiTs(j)*(1./Ve2C2**2)*(Qstar/Zlq)**2/1.E-6
     ELSE
      Aalpha= RMiMe/VRTiTs(j)*(1./Ve2C2**2)*(Qstar/Zlq)**2/Aux1 
     END IF
     DO SS= -1,1,2
      Qs= Qstar+SS/SQRT(2.*Aalpha)
      Qxp= Qs*SIN(Phip)
      Qzp= Qs*COS(Phip)
      CALL Locate(VQx,nqx,Qxp,ires)
      IQxLLT2i(iqx,kqz,i,S2,S1,SS)= ires
      CALL Locate(VQz,nqz,Qzp,ires)
      IQzLLT2i(iqx,kqz,i,S2,S1,SS)= ires
     END DO
     IF (Aux1<1.E-6) THEN
      Aalpha= 1.D0/VRTeTs(j)*(1./Ve2C2**2)*(Qstar/Zlq)**2/1.E-6
     ELSE
      Aalpha= 1.D0/VRTeTs(j)*(1./Ve2C2**2)*(Qstar/Zlq)**2/Aux1 
     END IF
     DO SS= -1,1,2
      Qs= Qstar+SS/SQRT(2.*Aalpha)
      Qxp= Qs*SIN(Phip)
      Qzp= Qs*COS(Phip)
      CALL Locate(VQx,nqx,Qxp,ires)
      IQxLLT2e(iqx,kqz,i,S2,S1,SS)= ires
      CALL Locate(VQz,nqz,Qzp,ires)
      IQzLLT2e(iqx,kqz,i,S2,S1,SS)= ires
     END DO
    END DO
   END DO
  END DO
 END DO
END DO

! T waves, scattering involving L and T waves
ALLOCATE (IQxTLT1(nqx2,nqz2,nph),IQzTLT1(nqx2,nqz2,nph))
ALLOCATE (IQxTLT2i(nqx2,nqz2,nph,-1:1,-1:1,-1:1),&
          IQzTLT2i(nqx2,nqz2,nph,-1:1,-1:1,-1:1))
ALLOCATE (IQxTLT2e(nqx2,nqz2,nph,-1:1,-1:1,-1:1),&
          IQzTLT2e(nqx2,nqz2,nph,-1:1,-1:1,-1:1))
DO iqx= 1,nqx2
 Qx= VQx2(iqx)
 DO kqz= 1,nqz2
  Qz= VQz2(kqz)
  Qx2= Qx**2
  Qz2= Qz**2
  Q2= Qx2+Qz2
  Q= SQRT(Q2)
  Phi= ACOS(Qz/Q)
  Zlq= ZL(Qx,Qz)
  Ztq= ZT(Qx,Qz)
  Qstar= SQRT(2./3./Ve2C2/VRTeTs(j))*Q
  DO i= 1,nph
   Phip= VPhip(i)
   Qxp= Qstar*SIN(Phip)
   Qzp= Qstar*COS(Phip)
   CALL Locate(VQx,nqx,Qxp,ires)
   IQxTLT1(iqx,kqz,i)= ires
   CALL Locate(VQz,nqz,Qzp,ires)
   IQzTLT1(iqx,kqz,i)= ires
   DO S2= -1,1,2
    DO S1= -1,1,2
     Aux1= (Q2+Qstar**2-2.*S1*Q*Qstar*COS(Phi-S1*S2*Phip))
     IF (Aux1<1.E-6) THEN
      Aalpha= RMiMe/VRTiTs(j)*(9./4.)*(Qstar/Ztq)**2/1.E-6
     ELSE
      Aalpha= RMiMe/VRTiTs(j)*(9./4)*(Qstar/Ztq)**2/Aux1 
     END IF
     DO SS= -1,1,2
      Qs= Qstar+SS/SQRT(2.*Aalpha)
      Qxp= Qs*SIN(Phip)
      Qzp= Qs*COS(Phip)
      CALL Locate(VQx,nqx,Qxp,ires)
      IQxTLT2i(iqx,kqz,i,S2,S1,SS)= ires
      CALL Locate(VQz,nqz,Qzp,ires)
      IQzTLT2i(iqx,kqz,i,S2,S1,SS)= ires
     END DO
     IF (Aux1<1.E-6) THEN
      Aalpha= 1.D0/VRTeTs(j)*(9./4.)*(Qstar/Ztq)**2/1.E-6
     ELSE
      Aalpha= 1.D0/VRTeTs(j)*(9./4)*(Qstar/Ztq)**2/Aux1 
     END IF
     DO SS= -1,1,2
      Qs= Qstar+SS/SQRT(2.*Aalpha)
      Qxp= Qs*SIN(Phip)
      Qzp= Qs*COS(Phip)
      CALL Locate(VQx,nqx,Qxp,ires)
      IQxTLT2e(iqx,kqz,i,S2,S1,SS)= ires
      CALL Locate(VQz,nqz,Qzp,ires)
      IQzTLT2e(iqx,kqz,i,S2,S1,SS)= ires
     END DO
    END DO
   END DO
  END DO
 END DO
END DO

! Resonant values of Q, spontaneous and induced scattering, approximated form:
! L waves, LL scattering:
ALLOCATE (IQQLLLs(nqx,nqz,nph),IQQLLLsi(nqx,nqz,nph,-1:1,-1:1,-1:1))
DO iqx= 1,nqx
 Qx= VQx(iqx)
 DO kqz= 1,nqz
  Qz= VQz(kqz)
  Q2= Qx**2+Qz**2
  Q= SQRT(Q2)
  Zlq= ZL(Qx,Qz)
  Qstar= Q
  DO i= 1,nph
   Phip= VPhip(i)
   Sphip= SIN(Phip)
   Cphip= COS(Phip)
   Qxp= Qstar*SIN(Phip)
   Qzp= Qstar*COS(Phip)
   Qp2= Qxp**2+Qzp**2
   Qp= SQRT(Qp2)
   CALL Locate(VQQ,nqmod,Qp,ires)
   IQQLLLs(iqx,kqz,i)= ires
   DO sigma= -1,1,2
    DO sigmap= -1,1,2
     !DO S1= -1,1,2
      Aux1= (Q2+Qstar**2-2.*sigma*sigmap*Q*Qstar*COS(Phip))
      IF (Aux1<1.E-6) THEN
       AuxNum= SQRT(1.E-6)
       Aalpha= RMiMe/VRTiTs(j)*(9./4.)*(Qstar/Zlq)**2/1.E-6 
      ELSE
       AuxNum= SQRT(Aux1)
       Aalpha= RMiMe/VRTiTs(j)*(9./4.)*(Qstar/Zlq)**2/Aux1 
      END IF
      Aux= (COS(Phip))**2*AuxNum
      DO SS= -1,1,2
       Qs= Qstar+SS/SQRT(2.*Aalpha)
       Qxp= Qs*Sphip
       Qzp= Qs*Cphip
       Qp2= Qxp**2+Qzp**2
       Qp= SQRT(Qp2)
       CALL Locate(VQQ,nqmod,Qp,ires)
       IQQLLLsi(iqx,kqz,i,sigma,sigmap,SS)= ires
      END DO
     !END DO
    END DO
   END DO
  END DO
 END DO
END DO
! T waves, LT scattering:
ALLOCATE (IQQTLTs(nqx2,nqz2,nph),IQQTLTsi(nqx2,nqz2,nph,-1:1,-1:1,-1:1))
ALLOCATE (IQQTLTse(nqx2,nqz2,nph,-1:1,-1:1,-1:1))
DO iqx= 1,nqx2
 Qx= VQx2(iqx)
 DO kqz= 1,nqz2
  Qz= VQz2(kqz)
  Q2= Qx**2+Qz**2
  Q= SQRT(Q2)
  Ztq= ZT(Qx,Qz)
  Qstar= SQRT(2./3./Ve2C2)*Q
  DO i= 1,nph
   Phip= VPhip(i)
   Sphip= SIN(Phip)
   Cphip= COS(Phip)
   Qxp= Qstar*SIN(Phip)
   Qzp= Qstar*COS(Phip)
   Qp2= Qxp**2+Qzp**2
   Qp= SQRT(Qp2)
   CALL Locate(VQQ,nqmod,Qp,ires)
   IQQTLTs(iqx,kqz,i)= ires
   DO sigma= -1,1,2
    DO sigmap= -1,1,2
     !DO S1= -1,1,2
      Aux1= (Q2+Qstar**2-2.*sigma*sigmap*Q*Qstar*COS(Phip))
      IF (Aux1<1.E-6) THEN
       AuxNum= SQRT(1.E-6)
       Aalpha= RMiMe/VRTiTs(j)*(9./4.)*(Qstar/Ztq)**2/1.E-6
      ELSE
       AuxNum= SQRT(Aux1)
       Aalpha= RMiMe/VRTiTs(j)*(9./4.)*(Qstar/Ztq)**2/Aux1 
      END IF
      Aux= (COS(Phip))**2*AuxNum
      DO SS= -1,1,2
       Qs= Qstar+SS/SQRT(2.*Aalpha)
       Qxp= Qs*Sphip
       Qzp= Qs*Cphip
       Qp2= Qxp**2+Qzp**2
       Qp= SQRT(Qp2)
       CALL Locate(VQQ,nqmod,Qp,ires)
       IQQTLTsi(iqx,kqz,i,sigma,sigmap,SS)= ires
      END DO
     !END DO
    END DO
   END DO
   DO sigma= -1,1,2
    DO sigmap= -1,1,2
     !DO S1= -1,1,2
      Aux1= (Q2+Qstar**2-2.*sigma*sigmap*Q*Qstar*COS(Phip))
      IF (Aux1<1.E-6) THEN
       AuxNum= SQRT(1.E-6)
       Aalpha= 1./VRTeTs(j)*(9./4.)*(Qstar/Ztq)**2/1.E-6
      ELSE
       AuxNum= SQRT(Aux1)
       Aalpha= 1./VRTeTs(j)*(9./4.)*(Qstar/Ztq)**2/Aux1 
      END IF
      Aux= (COS(Phip))**2*AuxNum
      DO SS= -1,1,2
       Qs= Qstar+SS/SQRT(2.*Aalpha)
       Qxp= Qs*Sphip
       Qzp= Qs*Cphip
       Qp2= Qxp**2+Qzp**2
       Qp= SQRT(Qp2)
       CALL Locate(VQQ,nqmod,Qp,ires)
       IQQTLTse(iqx,kqz,i,sigma,sigmap,SS)= ires
      END DO
     !END DO
    END DO
   END DO
  END DO
 END DO
END DO

! For the coefficients Ai and Dij:
DO m= 1,nuz
 Uz= VUz(m)
 !IF (Uz==0.) Uz= INT(SIGN(1.,Uz))*EpsMin
 IF (ABS(Uz).LT.1.E-16) Uz= INT(SIGN(1.D0,Uz))*EpsMin
 AbsUz= ABS(Uz)
 DO l= 1,nux
  Ux= VUx(l)
  !IF (Ux==0.) Ux= EpsMin
  IF (ABS(Ux).LT.1.E-16) Ux= EpsMin
  AbsUx= ABS(Ux)
  IF (AbsUz>=Ux) THEN
   DO sigma= -1,1,2
    DO i= 1,nqx
     Qx= VQx(i)
     Q1= (sigma*SQRT(VRNeNs(j))-Qx*Ux)/Uz
     IF (Q1>0.) THEN
      CALL Locate(VQz,nqz,Q1,ires)
      IQzr1D(l,m,i,sigma)= ires
      Qzr1D(l,m,i,sigma)= Q1
     ELSE
      IQzr1D(l,m,i,sigma)= -1
      Qzr1D(l,m,i,sigma)= 0.
     END IF
     Q2= (sigma*SQRT(VRNeNs(j))+Qx*Ux)/Uz
     IF (Q2>0.) THEN
      CALL Locate(VQz,nqz,Q2,ires)
      IQzr2D(l,m,i,sigma)= ires
      Qzr2D(l,m,i,sigma)= Q2
     ELSE
      IQzr2D(l,m,i,sigma)= -1
      Qzr2D(l,m,i,sigma)= 0.
     END IF
    END DO
   END DO
  ELSE
   DO sigma= -1,1,2
    DO i= 1,nqz
     Qz= VQz(i)
     Q1= (sigma*SQRT(VRNeNs(j))-Qz*Uz)/Ux
     IF (Q1>0.) THEN
      CALL Locate(VQx,nqx,Q1,ires)
      IQxr1D(l,m,i,sigma)= ires
      Qxr1D(l,m,i,sigma)= Q1
     ELSE
      IQxr1D(l,m,i,sigma)= -1
      Qxr1D(l,m,i,sigma)= 0.
     END IF
     Q2= (-sigma*SQRT(VRNeNs(j))+Qz*Uz)/Ux
     IF (Q2>0.) THEN
      CALL Locate(VQx,nqx,Q2,ires)
      IQxr2D(l,m,i,sigma)= ires
      Qxr2D(l,m,i,sigma)= Q2
     ELSE
      IQxr2D(l,m,i,sigma)= -1
      Qxr2D(l,m,i,sigma)= 0.
     END IF
    END DO
   END DO
  END IF
END DO
END DO

IF(Lemis=="Yes") THEN
! L waves, contribution due to spontaneous and induced emission:
! L waves, spontaneous and induced emission:
ALLOCATE (UzrpLql(nqx,nqz,nux,-1:1),UzrmLql(nqx,nqz,nux,-1:1))
ALLOCATE (UxrpLql(nqx,nqz,nuz,-1:1),UxrmLql(nqx,nqz,nuz,-1:1))
ALLOCATE (IUzrpLql(nqx,nqz,nux,-1:1),IUzrmLql(nqx,nqz,nux,-1:1))
ALLOCATE (IUxrpLql(nqx,nqz,nuz,-1:1),IUxrmLql(nqx,nqz,nuz,-1:1))
UzrpLql= 0.
UzrmLql= 0.
UxrpLql= 0.
UxrmLql= 0.
IUzrpLql= 0
IUzrmLql= 0
IUxrpLql= 0
IUxrmLql= 0
DO iqx= 1,nqx
 Qx= VQx(iqx)
 DO kqz= 1,nqz
  Qz= VQz(kqz)
  Qx2= Qx**2
  Qz2= Qz**2
  Q2= Qx2+Qz2
  Q= SQRT(Q2)
  Zlq= ZL(Qx,ABS(Qz))
  IF (Qz>Qx) THEN
   DO l= 1,nux
    Ux= VUx(l)
    DO sigma= -1,1,2
     Uresp= (sigma*Zlq+Qx*Ux)/Qz
     Uresm= (sigma*Zlq-Qx*Ux)/Qz
     UzrpLql(iqx,kqz,l,sigma)= Uresp
     UzrmLql(iqx,kqz,l,sigma)= Uresm
     CALL Locate(VUz,nuz,Uresp,ires)
     IUzrpLql(iqx,kqz,l,sigma)= ires
     CALL Locate(VUz,nuz,Uresm,ires)
     IUzrmLql(iqx,kqz,l,sigma)= ires
    END DO
   END DO
  ELSE
   DO l= 1,nuz
    Uz= VUz(l)
    DO sigma= -1,1,2
     Uresp= (sigma*Zlq-Qz*Uz)/Qx
     Uresm= -(sigma*Zlq-Qz*Uz)/Qx
     CALL Locate(VUx,nux,Uresp,ires)
     IUxrpLql(iqx,kqz,l,sigma)= ires
     UxrpLql(iqx,kqz,l,sigma)= Uresp
     CALL Locate(VUx,nux,Uresm,ires)
     IUxrmLql(iqx,kqz,l,sigma)= ires
     UxrmLql(iqx,kqz,l,sigma)= Uresm
    END DO
   END DO
  END IF
 END DO
END DO
ELSE
END IF

IF(LdecayLS=="Yes") THEN
! L waves, decay associated to L and S waves:
! L waves, decay, L and S waves:
ALLOCATE (QxpLLSd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
ALLOCATE (IQxpLLSd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
ALLOCATE (IQxpdifLLSd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
QxpLLSd= 0.
IQxpLLSd= 0
IQxpdifLLSd= 0
DO iqx= 1,nqx
 Qx= VQx(iqx)
 DO kqz= 1,nqz
  Qz= VQz(kqz)
  Qx2= Qx**2
  Qz2= Qz**2
  Q2= Qx2+Qz2
  DO Asp= -1,1,2
   DO Aspp= -1,1,2
    DO S1= -1,1,2
     DO S2= -1,1,2
      Aux1_Rcd(1)= Qx
      Aux1_Rcd(2)= Qz
      Aux2_Rcd(1)= Asp
      Aux2_Rcd(2)= Aspp
      Aux2_Rcd(3)= S1
      Aux2_Rcd(4)= S2
!      IF (Qz>Qx) THEN
!      IF (Qz<0.) THEN
!      ELSE
       DO kqzp= 1,nqz
        Qzp= VQz(kqzp)
        Aux1_Rcd(3)= Qzp
        X1= 0.
        X2= 1.0
        CALL ZBRAC2(Funcx_RcdLLS,X1,X2,Succes) 
        IF (Succes .eqv. .true.) THEN 
         Qxp= RTBIS(Funcx_RcdLLS,X1,X2,Xacc)
        ELSE
         Qxp= 0.
        END IF
        IF(Qxp>0.)THEN
         Qxdif= Qx-S1*Qxp
         CALL Locate(VQx,nqx,Qxp,ires)
         IQxpLLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
         QxpLLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= Qxp
         CALL Locate(VQx,nqx,ABS(Qxdif),ires)
         IQxpdifLLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
        ELSE
         IQxpLLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= -1
         QxpLLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= 0.D0
         IQxpdifLLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= -1
        END IF
       END DO
!      END IF
     END DO
    END DO
   END DO
  END DO
 END DO
END DO
ELSE
END IF

IF(LdecayLT=="Yes") THEN
! L waves, decay associated to L and T waves:
! L waves, decay, L and T waves:
ALLOCATE (QxpLLTd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
ALLOCATE (IQxpLLTd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
ALLOCATE (IQxpdifLLTd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
QxpLLTd= 0.
IQxpLLTd= 0
IQxpdifLLTd= 0
DO iqx= 1,nqx
 Qx= VQx(iqx)
 DO kqz= 1,nqz
  Qz= VQz(kqz)
  Qx2= Qx**2
  Qz2= Qz**2
  Q2= Qx2+Qz2
  DO Asp= -1,1,2
   DO Aspp= -1,1,2
    DO S1= -1,1,2
     DO S2= -1,1,2
      Aux1_Rcd(1)= Qx
      Aux1_Rcd(2)= Qz
      Aux2_Rcd(1)= Asp
      Aux2_Rcd(2)= Aspp
      Aux2_Rcd(3)= S1
      Aux2_Rcd(4)= S2
!      IF (Qz>Qx) THEN
!      IF (Qz<0.) THEN
!      ELSE
       DO kqzp= 1,nqz
        Qzp= VQz(kqzp)
        Aux1_Rcd(3)= Qzp
        X1= 0.
        X2= 1.0
        CALL ZBRAC2(Funcx_RcdLLT,X1,X2,Succes) 
        IF (Succes .eqv. .true.) THEN 
         Qxp= RTBIS(Funcx_RcdLLT,X1,X2,Xacc)
        ELSE
         Qxp= 0.
        END IF
        IF(Qxp>0.)THEN
         Qxdif= Qx-S1*Qxp
         CALL Locate(VQx,nqx,Qxp,ires)
         IQxpLLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
         QxpLLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= Qxp
         CALL Locate(VQx,nqx,ABS(Qxdif),ires)
         IQxpdifLLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
        ELSE
         IQxpLLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= -1
         QxpLLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= 0.D0
         IQxpdifLLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= -1
        END IF
       END DO
!      END IF
     END DO
    END DO
   END DO
  END DO
 END DO
END DO
ELSE
END IF

IF(LdecayST=="Yes") THEN
! L waves, decay associated to S and T waves:
! L waves, decay, S and T waves:
ALLOCATE (QxpLSTd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
ALLOCATE (IQxpLSTd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
ALLOCATE (IQxpdifLSTd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
QxpLSTd= 0.
IQxpLSTd= 0
IQxpdifLSTd= 0
DO iqx= 1,nqx
 Qx= VQx(iqx)
 DO kqz= 1,nqz
  Qz= VQz(kqz)
  Qx2= Qx**2
  Qz2= Qz**2
  Q2= Qx2+Qz2
  DO Asp= -1,1,2
   DO Aspp= -1,1,2
    DO S1= -1,1,2
     DO S2= -1,1,2
      Aux1_Rcd(1)= Qx
      Aux1_Rcd(2)= Qz
      Aux2_Rcd(1)= Asp
      Aux2_Rcd(2)= Aspp
      Aux2_Rcd(3)= S1
      Aux2_Rcd(4)= S2
!      IF (Qz>Qx) THEN
!      IF (Qz<0.) THEN
!      ELSE
       DO kqzp= 1,nqz
        Qzp= VQz(kqzp)
        Aux1_Rcd(3)= Qzp
        X1= 0.
        X2= 1.0
        CALL ZBRAC2(Funcx_RcdLST,X1,X2,Succes) 
        IF (Succes .eqv. .true.) THEN 
         Qxp= RTBIS(Funcx_RcdLST,X1,X2,Xacc)
        ELSE
         Qxp= 0.
        END IF
        IF(Qxp>0.)THEN
         Qxdif= Qx-S1*Qxp
         CALL Locate(VQx,nqx,Qxp,ires)
         IQxpLSTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
         QxpLSTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= Qxp
         CALL Locate(VQx,nqx,ABS(Qxdif),ires)
         IQxpdifLSTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
        ELSE
         IQxpLSTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= -1
         QxpLLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= 0.D0
         IQxpdifLSTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= -1
        END IF
       END DO
!      END IF
     END DO
    END DO
   END DO
  END DO
 END DO
END DO
ELSE
END IF

IF(LdecayTT=="Yes") THEN
! L waves, decay associated to T and T waves:
! L waves, decay, T and T waves:
ALLOCATE (QxpLTTd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
ALLOCATE (IQxpLTTd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
ALLOCATE (IQxpdifLTTd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
QxpLTTd= 0.
IQxpLTTd= 0
IQxpdifLTTd= 0
DO iqx= 1,nqx
 Qx= VQx(iqx)
 DO kqz= 1,nqz
  Qz= VQz(kqz)
  Qx2= Qx**2
  Qz2= Qz**2
  Q2= Qx2+Qz2
  DO Asp= -1,1,2
   DO Aspp= -1,1,2
    DO S1= -1,1,2
     DO S2= -1,1,2
      Aux1_Rcd(1)= Qx
      Aux1_Rcd(2)= Qz
      Aux2_Rcd(1)= Asp
      Aux2_Rcd(2)= Aspp
      Aux2_Rcd(3)= S1
      Aux2_Rcd(4)= S2
!      IF (Qz>Qx) THEN
!      IF (Qz<0.) THEN
!      ELSE
       DO kqzp= 1,nqz
        Qzp= VQz(kqzp)
        Aux1_Rcd(3)= Qzp
        X1= 0.
        X2= 1.0
        CALL ZBRAC2(Funcx_RcdLTT,X1,X2,Succes) 
        IF (Succes .eqv. .true.) THEN 
         Qxp= RTBIS(Funcx_RcdLTT,X1,X2,Xacc)
        ELSE
         Qxp= 0.
        END IF
        IF(Qxp>0.)THEN
         Qxdif= Qx-S1*Qxp
         CALL Locate(VQx,nqx,Qxp,ires)
         IQxpLTTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
         QxpLTTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= Qxp
         CALL Locate(VQx,nqx,ABS(Qxdif),ires)
         IQxpdifLTTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
        ELSE
         IQxpLTTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= -1
         QxpLTTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= 0.D0
         IQxpdifLTTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= -1
        END IF
       END DO
!      END IF
     END DO
    END DO
   END DO
  END DO
 END DO
END DO
ELSE
END IF

IF(Semis=="Yes") THEN
! S waves, contribution due to spontaneous and induced emission:
! S waves, spontaneous and induced emission:
ALLOCATE (UzrpSql(nqx,nqz,nux,-1:1),UzrmSql(nqx,nqz,nux,-1:1))
ALLOCATE (UxrpSql(nqx,nqz,nuz,-1:1),UxrmSql(nqx,nqz,nuz,-1:1))
ALLOCATE (IUzrpSql(nqx,nqz,nux,-1:1),IUzrmSql(nqx,nqz,nux,-1:1))
ALLOCATE (IUxrpSql(nqx,nqz,nuz,-1:1),IUxrmSql(nqx,nqz,nuz,-1:1))
UzrpSql= 0.
UzrmSql= 0.
UxrpSql= 0.
UxrmSql= 0.
IUzrpSql= 0
IUzrmSql= 0
IUxrpSql= 0
IUxrmSql= 0
DO iqx= 1,nqx
 Qx= VQx(iqx)
 DO kqz= 1,nqz
  Qz= VQz(kqz)
  Q2= Qx**2+Qz**2
  Q= SQRT(Q2)
  Muq= Q**3*AA/2.
  Zlq= ZL(Qx,Qz)
  Zsq= ZS(Qx,Qz)
  IF (Qz>Qx) THEN 
   DO l= 1,nux
    Ux= VUx(l)
    DO sigma= -1,1,2
     Uresp= (sigma*Zsq+Qx*Ux)/Qz
     Uresm= (sigma*Zsq-Qx*Ux)/Qz
     CALL Locate(VUz,nuz,Uresp,ires)
     IUzrpSql(iqx,kqz,l,sigma)= ires
     UzrpSql(iqx,kqz,l,sigma)= Uresp
     CALL Locate(VUz,nuz,Uresm,ires)
     IUzrmSql(iqx,kqz,l,sigma)= ires
     UzrmSql(iqx,kqz,l,sigma)= Uresm
    END DO
   END DO
  ELSE
   DO l= 1,nuz
    Uz= VUz(l)
    DO sigma= -1,1,2
     Uresp= (sigma*Zsq-Qz*Uz)/Qx
     Uresm= -(sigma*Zsq-Qz*Uz)/Qx
     CALL Locate(VUx,nux,Uresp,ires)
     IUxrpSql(iqx,kqz,l,sigma)= ires
     UxrpSql(iqx,kqz,l,sigma)= Uresp
     CALL Locate(VUx,nux,Uresm,ires)
     IUxrmSql(iqx,kqz,l,sigma)= ires
     UxrmSql(iqx,kqz,l,sigma)= Uresm
    END DO
   END DO
  END IF
 END DO
END DO
ELSE
END IF

IF(SdecayLL=="Yes") THEN
! S waves, decay associated to L and L waves:
! S waves, decay, L and L waves:
ALLOCATE (QxpSLLd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
ALLOCATE (IQxpSLLd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
ALLOCATE (IQxpdifSLLd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
QxpSLLd= 0.
IQxpSLLd= 0
IQxpdifSLLd= 0
DO iqx= 1,nqx
 Qx= VQx(iqx)
 DO kqz= 1,nqz
  Qz= VQz(kqz)
  Qx2= Qx**2
  Qz2= Qz**2
  Q2= Qx2+Qz2
  DO Asp= -1,1,2
   DO Aspp= -1,1,2
    DO S1= -1,1,2
     DO S2= -1,1,2
      Aux1_Rcd(1)= Qx
      Aux1_Rcd(2)= Qz
      Aux2_Rcd(1)= Asp
      Aux2_Rcd(2)= Aspp
      Aux2_Rcd(3)= S1
      Aux2_Rcd(4)= S2
!      IF (Qz>Qx) THEN
!      IF (Qz<0.) THEN
!      ELSE
       DO kqzp= 1,nqz
        Qzp= VQz(kqzp)
        Aux1_Rcd(3)= Qzp
        X1= 0.
        X2= 1.0
        CALL ZBRAC2(Funcx_RcdSLL,X1,X2,Succes) 
        IF (Succes .eqv. .true.) THEN 
         Qxp= RTBIS(Funcx_RcdSLL,X1,X2,Xacc)
        ELSE
         Qxp= 0.
        END IF
        IF(Qxp>0.)THEN
         Qxdif= Qx-S1*Qxp
         CALL Locate(VQx,nqx,Qxp,ires)
         IQxpSLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
         QxpSLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= Qxp
         CALL Locate(VQx,nqx,ABS(Qxdif),ires)
         IQxpdifSLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
        ELSE
         IQxpSLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= -1
         QxpSLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= 0.D0
         IQxpdifSLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= -1
        END IF
       END DO
!      END IF
     END DO
    END DO
   END DO
  END DO
 END DO
END DO
ELSE
END IF

IF(SdecayLT=="Yes") THEN
! S waves, decay associated to L and T waves:
! S waves, decay, L and T waves:
ALLOCATE (QxpSLTd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
ALLOCATE (IQxpSLTd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
ALLOCATE (IQxpdifSLTd(nqx,nqz,nqz,-1:1,-1:1,-1:1,-1:1))
QxpSLTd= 0.
IQxpSLTd= 0
IQxpdifSLTd= 0
DO iqx= 1,nqx
 Qx= VQx(iqx)
 DO kqz= 1,nqz
  Qz= VQz(kqz)
  Qx2= Qx**2
  Qz2= Qz**2
  Q2= Qx2+Qz2
  DO Asp= -1,1,2
   DO Aspp= -1,1,2
    DO S1= -1,1,2
     DO S2= -1,1,2
      Aux1_Rcd(1)= Qx
      Aux1_Rcd(2)= Qz
      Aux2_Rcd(1)= Asp
      Aux2_Rcd(2)= Aspp
      Aux2_Rcd(3)= S1
      Aux2_Rcd(4)= S2
!      IF (Qz>Qx) THEN
!      IF (Qz<0.) THEN
!      ELSE
       DO kqzp= 1,nqz
        Qzp= VQz(kqzp)
        Aux1_Rcd(3)= Qzp
        X1= 0.
        X2= 1.0
        CALL ZBRAC2(Funcx_RcdSLT,X1,X2,Succes) 
        IF (Succes .eqv. .true.) THEN 
         Qxp= RTBIS(Funcx_RcdSLT,X1,X2,Xacc)
        ELSE
         Qxp= 0.
        END IF
        IF(Qxp>0.)THEN
         Qxdif= Qx-S1*Qxp
         CALL Locate(VQx,nqx,Qxp,ires)
         IQxpSLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
         QxpSLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= Qxp
         CALL Locate(VQx,nqx,ABS(Qxdif),ires)
         IQxpdifSLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
        ELSE
         IQxpSLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= -1
         QxpSLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= 0.D0
         IQxpdifSLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= -1
        END IF
       END DO
!      END IF
     END DO
    END DO
   END DO
  END DO
 END DO
END DO
ELSE
END IF

IF(TdecayLL=="Yes") THEN
! T waves, decay associated to L and L waves:
! T waves, decay, L and L waves:
ALLOCATE (QxpTLLd(nqx2,nqz2,nqz2,-1:1,-1:1,-1:1,-1:1))
ALLOCATE (IQxpTLLd(nqx2,nqz2,nqz2,-1:1,-1:1,-1:1,-1:1))
ALLOCATE (IQxpdifTLLd(nqx2,nqz2,nqz2,-1:1,-1:1,-1:1,-1:1))
QxpTLLd= 0.
IQxpTLLd= 0
IQxpdifTLLd= 0
DO iqx= 1,nqx2
 Qx= VQx2(iqx)
 DO kqz= 1,nqz2
  Qz= VQz2(kqz)
  Qx2= Qx**2
  Qz2= Qz**2
  Q2= Qx2+Qz2
  DO Asp= -1,1,2
   DO Aspp= -1,1,2
    DO S1= -1,1,2
     DO S2= -1,1,2
      Aux1_Rcd(1)= Qx
      Aux1_Rcd(2)= Qz
      Aux2_Rcd(1)= Asp
      Aux2_Rcd(2)= Aspp
      Aux2_Rcd(3)= S1
      Aux2_Rcd(4)= S2
!      IF (Qz>Qx) THEN
!      IF (Qz<0.) THEN
!      ELSE
       DO kqzp= 1,nqz2
        Qzp= VQz2(kqzp)
        Aux1_Rcd(3)= Qzp
        X1= 0.
        X2= 1.0
        CALL ZBRAC2(Funcx_RcdTLL,X1,X2,Succes) 
        IF (Succes .eqv. .true.) THEN 
         Qxp= RTBIS(Funcx_RcdTLL,X1,X2,Xacc)
        ELSE
         Qxp= 0.
        END IF
        IF(Qxp>0.)THEN
         Qxdif= Qx-S1*Qxp
         CALL Locate(VQx2,nqx2,Qxp,ires)
         IQxpTLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
         QxpTLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= Qxp
         CALL Locate(VQx2,nqx2,ABS(Qxdif),ires)
         IQxpdifTLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
        ELSE
         IQxpTLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= -1
         QxpTLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= 0.D0
         IQxpdifTLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= -1
        END IF
       END DO
!      END IF
     END DO
    END DO
   END DO
  END DO
 END DO
END DO
ELSE
END IF

IF(TdecayLS=="Yes") THEN
! T waves, decay associated to L and S waves:
! T waves, decay, L and S waves:
ALLOCATE (QxpTLSd(nqx2,nqz2,nqz2,-1:1,-1:1,-1:1,-1:1))
ALLOCATE (IQxpTLSd(nqx2,nqz2,nqz2,-1:1,-1:1,-1:1,-1:1))
ALLOCATE (IQxpdifTLSd(nqx2,nqz2,nqz2,-1:1,-1:1,-1:1,-1:1))
QxpTLSd= 0.
IQxpTLSd= 0
IQxpdifTLSd= 0
DO iqx= 1,nqx2
 Qx= VQx2(iqx)
 DO kqz= 1,nqz2
  Qz= VQz2(kqz)
  Qx2= Qx**2
  Qz2= Qz**2
  Q2= Qx2+Qz2
  DO Asp= -1,1,2
   DO Aspp= -1,1,2
    DO S1= -1,1,2
     DO S2= -1,1,2
      Aux1_Rcd(1)= Qx
      Aux1_Rcd(2)= Qz
      Aux2_Rcd(1)= Asp
      Aux2_Rcd(2)= Aspp
      Aux2_Rcd(3)= S1
      Aux2_Rcd(4)= S2
!      IF (Qz>Qx) THEN
!      IF (Qz<0.) THEN
!      ELSE
       DO kqzp= 1,nqz2
        Qzp= VQz2(kqzp)
        Aux1_Rcd(3)= Qzp
        X1= 0.
        X2= 1.0
        CALL ZBRAC2(Funcx_RcdTLS,X1,X2,Succes) 
        IF (Succes .eqv. .true.) THEN 
         Qxp= RTBIS(Funcx_RcdTLS,X1,X2,Xacc)
        ELSE
         Qxp= 0.
        END IF
        IF(Qxp>0.)THEN
         Qxdif= Qx-S1*Qxp
         CALL Locate(VQx2,nqx2,Qxp,ires)
         IQxpTLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
         QxpTLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= Qxp
         CALL Locate(VQx2,nqx2,ABS(Qxdif),ires)
         IQxpdifTLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
        ELSE
         IQxpTLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= -1
         QxpTLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= 0.D0
         IQxpdifTLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= -1
        END IF
       END DO
!      END IF
     END DO
    END DO
   END DO
  END DO
 END DO
END DO
ELSE
END IF

IF(TdecayTL=="Yes") THEN
! T waves, decay associated to T and L waves:
! T waves, decay, T and L waves:
ALLOCATE (QxpTTLd(nqx2,nqz2,nqz2,-1:1,-1:1,-1:1,-1:1))
ALLOCATE (IQxpTTLd(nqx2,nqz2,nqz2,-1:1,-1:1,-1:1,-1:1))
ALLOCATE (IQxpdifTTLd(nqx2,nqz2,nqz2,-1:1,-1:1,-1:1,-1:1))
QxpTTLd= 0.
IQxpTTLd= 0
IQxpdifTTLd= 0
DO iqx= 1,nqx2
 Qx= VQx2(iqx)
 DO kqz= 1,nqz2
  Qz= VQz2(kqz)
  Qx2= Qx**2
  Qz2= Qz**2
  Q2= Qx2+Qz2
  DO Asp= -1,1,2
   DO Aspp= -1,1,2
    DO S1= -1,1,2
     DO S2= -1,1,2
      Aux1_Rcd(1)= Qx
      Aux1_Rcd(2)= Qz
      Aux2_Rcd(1)= Asp
      Aux2_Rcd(2)= Aspp
      Aux2_Rcd(3)= S1
      Aux2_Rcd(4)= S2
!      IF (Qz>Qx) THEN
!      IF (Qz<0.) THEN
!      ELSE
       DO kqzp= 1,nqz2
        Qzp= VQz2(kqzp)
        Aux1_Rcd(3)= Qzp
        X1= 0.
        X2= 1.0
        CALL ZBRAC2(Funcx_RcdTTL,X1,X2,Succes) 
        IF (Succes .eqv. .true.) THEN 
         Qxp= RTBIS(Funcx_RcdTTL,X1,X2,Xacc)
        ELSE
         Qxp= 0.
        END IF
        IF(Qxp>0.)THEN
         Qxdif= Qx-S1*Qxp
         CALL Locate(VQx2,nqx2,Qxp,ires)
         IQxpTTLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
         QxpTTLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= Qxp
         CALL Locate(VQx2,nqx2,ABS(Qxdif),ires)
         IQxpdifTTLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= ires
        ELSE
         IQxpTTLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= -1
         QxpTTLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= 0.D0
         IQxpdifTTLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)= -1
        END IF
       END DO
!      END IF
     END DO
    END DO
   END DO
  END DO
 END DO
END DO
ELSE
END IF

RETURN
END SUBROUTINE Res_Cond

SUBROUTINE Coef_A
USE Common_Params
USE Common_Arrays
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL*8 :: Ux,Uz,Qx,Qx2,Qz,Qz2
REAL*8 :: AbsUz,AbsUx
REAL*8 :: Q1,Q12
REAL*8 :: Q2,Q22
REAL*8 :: Auxx,Auxz,AuxSE
REAL*8, DIMENSION(nqx) :: VxAx,VxAz
REAL*8, DIMENSION(nqz) :: VzAx,VzAz
INTEGER :: i,j,l,m,sigma

!   Initialization of the Ai:
Ax= 0.
Az= 0.

IF(SpontEmis=="Yes") THEN
 AuxSE= 1.E0
ELSE
 AuxSE= 0.E0
END IF

j= 1   ! Auxiliary quantity, for the spatial profile.
DO m= 1,nuz
 Uz= VUz(m)
 !IF (Uz==0.) Uz= INT(SIGN(1.,Uz))*EpsMin
 IF (ABS(Uz).LT.EpsMin) Uz= INT(SIGN(1.D0,Uz))*EpsMin
 AbsUz= ABS(Uz)
 DO l= 1,nux
  Ux= VUx(l)
  !IF (Ux==0.) Ux= EpsMin
  IF (ABS(Ux).LT.EpsMin) Ux= EpsMin
  AbsUx= ABS(Ux)
  IF (AbsUz>=Ux) THEN
   DO sigma= -1,1,2
    VxAx= 0.
    VxAz= 0.
    DO i= 1,nqx
     Qx= VQx(i)
     Qx2= Qx**2
     Q1= Qzr1D(l,m,i,sigma)
     IF (Q1>0.) THEN
      Q12= Q1**2
      VxAx(i)= VxAx(i)+(Qx*Ux+Q1*Uz)/(Qx2+Q12) * Qx
      VxAz(i)= VxAz(i)+(Qx*Ux+Q1*Uz)/(Qx2+Q12) * Q1
     ELSE
      VxAx(i)= VxAx(i)+0.
      VxAz(i)= VxAz(i)+0.
     END IF
     Q2= Qzr2D(l,m,i,sigma)
     IF (Q2>0.) THEN
      Q22= Q2**2
      VxAx(i)= VxAx(i)+(-Qx*Ux+Q2*Uz)/(Qx2+Q22) * (-Qx)
      VxAz(i)= VxAz(i)+(-Qx*Ux+Q2*Uz)/(Qx2+Q22) * Q2
     ELSE
      VxAx(i)= VxAx(i)+0.
      VxAz(i)= VxAz(i)+0.
     END IF
    END DO
    CALL Simpson(VQx,VxAx,nqx,Auxx)
    CALL Simpson(VQx,VxAz,nqx,Auxz)
    Ax(l,m)= Ax(l,m)+(2.*Geff/AbsUz)*Auxx * AuxSE
    Az(l,m)= Az(l,m)+(2.*Geff/AbsUz)*Auxz * AuxSE
   END DO
  ELSE
   DO sigma= -1,1,2
    VzAx= 0.
    VzAz= 0.
    DO i= 1,nqz
     Qz= VQz(i)
     Qz2= Qz**2
     Q1= Qxr1D(l,m,i,sigma)
     IF (Q1>0.) THEN
      Q12= Q1**2
      VzAx(i)= VzAx(i)+(Q1*Ux+Qz*Uz)/(Q12+Qz2) * Q1
      VzAz(i)= VzAz(i)+(Q1*Ux+Qz*Uz)/(Q12+Qz2) * Qz
     ELSE
      VzAx(i)= VzAx(i)+0.
      VzAz(i)= VzAz(i)+0.
     END IF
     Q2= Qxr2D(l,m,i,sigma)
     IF (Q2>0.) THEN
      Q22= Q2**2
      VzAx(i)= VzAx(i)+(-Q2*Ux+Qz*Uz)/(Q22+Qz2) * (-Q2)
      VzAz(i)= VzAz(i)+(-Q2*Ux+Qz*Uz)/(Q22+Qz2) * Qz
     ELSE
      VzAx(i)= VzAx(i)+0.
      VzAz(i)= VzAz(i)+0.
     END IF
    END DO
    CALL Simpson(VQz,VzAx,nqz,Auxx)
    CALL Simpson(VQz,VzAz,nqz,Auxz)
    Ax(l,m)= Ax(l,m)+(2.*Geff/AbsUx)*Auxx * AuxSE
    Az(l,m)= Az(l,m)+(2.*Geff/AbsUx)*Auxz * AuxSE
   END DO
  END IF
 END DO
END DO

RETURN
END SUBROUTINE Coef_A

SUBROUTINE Coef_D
USE Common_Params
USE Common_Arrays
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL*8 :: Ux,Uz,Qx,Qx2,Qz,Qz2
REAL*8 :: AbsUz,AbsUx
REAL*8 :: EsLq1,EsLq2
REAL*8 :: Q1,Q12
REAL*8 :: Q2,Q22
REAL*8 :: Auxxx,Auxxz,Auxzz
REAL*8 :: HBrL,HBrS,HGcL,HGcS
REAL*8, DIMENSION(nqx,nqz) :: Iwave
REAL*8, DIMENSION(nqz) :: Vxq
REAL*8, DIMENSION(nqx) :: Vzq
REAL*8, DIMENSION(nqx) :: VxDxx,VxDxz,VxDzz
REAL*8, DIMENSION(nqz) :: VzDxx,VzDxz,VzDzz
INTEGER :: i,l,m,sigma
INTEGER :: ires1,ires2
INTEGER :: j

j= 1   ! Auxiliary quantity, for the spatial profile.
!   Initialization of the Dij:
Dxx= 0.
Dxz= 0.
Dzx= 0.
Dzz= 0.

DO m= 1,nuz
 Uz= VUz(m)
 !IF (Uz==0.) Uz= INT(SIGN(1.,Uz))*EpsMin
 IF (ABS(Uz).lt.EpsMin) Uz= INT(SIGN(1.D0,Uz))*EpsMin
 AbsUz= ABS(Uz)
 DO l= 1,nux
  Ux= VUx(l)
  !IF (Ux==0.) Ux= EpsMin
  IF (ABS(Ux).lt.EpsMin) Ux= EpsMin
  AbsUx= ABS(Ux)
  IF (AbsUz>=Ux) THEN
   DO sigma= -1,1,2
    VxDxx= 0.
    VxDxz= 0.
    VxDzz= 0.
    IF (sigma==1) THEN
     Iwave= ILp
    ELSE
     Iwave= ILm
    END IF
    DO i= 1,nqx
     Qx= VQx(i)
     Qx2= Qx**2
     Q1= Qzr1D(l,m,i,sigma)
     ires1= IQzr1D(l,m,i,sigma)
     SELECT CASE(ires1)
      CASE(-1)   ! Qz<0
       EsLq1= 0.D0
      CASE(0)
       EsLq1= Iwave(i,1)
      CASE(1:)   ! ires1 >= 1
       IF (ires1<nqz) THEN
        Vxq(:)= Iwave(i,:)
        CALL Aitp1d2(nqz,VQz,Vxq,Q1,EsLq1,ires1)
       ELSE
        HBrL=0. ! The new effects are not relevant for large Q and are neglected
        HBrS=0.
        HGcL=0.
        HGcS=0.
       CALL Iwave_Init(Qx,Q1,EsLq1,"L",HBrL,HBrS,HGcL,HGcS)
       END IF
      CASE DEFAULT
       OPEN(98,FILE='Warning_Select_Case.wt')
       WRITE(98,*) ' ires1= ',ires1
       WRITE(98,*) ' ires1 must be integer greater or equal -1 !!'
       CLOSE(98)
     END SELECT
     Q12= Q1*Q1
     VxDxx(i)= VxDxx(i)+Qx2/(Qx2+Q12) * EsLq1
     VxDxz(i)= VxDxz(i)+Qx*Q1/(Qx2+Q12) * EsLq1
     VxDzz(i)= VxDzz(i)+Q12/(Qx2+Q12) * EsLq1
     Q2= Qzr2D(l,m,i,sigma)
     ires2= IQzr2D(l,m,i,sigma)
     SELECT CASE(ires2)
      CASE(-1)   ! Qz<0
       EsLq2= 0.D0
      CASE(0)
       EsLq2= Iwave(i,1)
      CASE(1:)   ! ires2 >= 1
       IF (ires2<nqz) THEN
        Vxq(:)= Iwave(i,:)
        CALL Aitp1d2(nqz,VQz,Vxq,Q2,EsLq2,ires2)
       ELSE
        HBrL=0. ! The new effects are not relevant for large Q and are neglected
        HBrS=0.
        HGcL=0.
        HGcS=0.
       CALL Iwave_Init(Qx,Q2,EsLq2,"L",HBrL,HBrS,HGcL,HGcS)
       END IF
      CASE DEFAULT
       OPEN(98,FILE='Warning_Select_Case.wt')
       WRITE(98,*) ' ires2= ',ires2
       WRITE(98,*) ' ires2 must be integer greater or equal -1 !!'
       CLOSE(98)
     END SELECT
     Q22= Q2*Q2
     VxDxx(i)= VxDxx(i)+Qx2/(Qx2+Q22) * EsLq2
     VxDxz(i)= VxDxz(i)+(-Qx*Q2)/(Qx2+Q22) * EsLq2
     VxDzz(i)= VxDzz(i)+Q22/(Qx2+Q22) * EsLq2
    END DO
    CALL Simpson(VQx,VxDxx,nqx,Auxxx)
    CALL Simpson(VQx,VxDxz,nqx,Auxxz)
    CALL Simpson(VQx,VxDzz,nqx,Auxzz)
    Dxx(l,m)= Dxx(l,m)+(2./AbsUz)*Auxxx
    Dxz(l,m)= Dxz(l,m)+(2./AbsUz)*Auxxz
    Dzx(l,m)= Dzx(l,m)+(2./AbsUz)*Auxxz
    Dzz(l,m)= Dzz(l,m)+(2./AbsUz)*Auxzz
   END DO
  ELSE
   DO sigma= -1,1,2
    VzDxx= 0.
    VzDxz= 0.
    VzDzz= 0.
    IF (sigma==1) THEN
     Iwave= ILp
    ELSE
     Iwave= ILm
    END IF
    DO i= 1,nqz
     Qz= VQz(i)
     Qz2= Qz**2
     Q1= Qxr1D(l,m,i,sigma)
     ires1= IQxr1D(l,m,i,sigma)
     SELECT CASE(ires1)
      CASE(-1)   ! Q<0
       EsLq1= 0.D0
      CASE(0)
       EsLq1= Iwave(1,i)
      CASE(1:)   ! ires1 >= 1
       IF (ires1<nqx) THEN
        Vzq(:)= Iwave(:,i)
        CALL Aitp1d2(nqx,VQx,Vzq,Q1,EsLq1,ires1)
       ELSE
        HBrL=0. ! The new effects are not relevant for large Q and are neglected
        HBrS=0.
        HGcL=0.
        HGcS=0.
       CALL Iwave_Init(Q1,Qz,EsLq1,"L",HBrL,HBrS,HGcL,HGcS)
       END IF
      CASE DEFAULT
       OPEN(98,FILE='Warning_Select_Case.wt')
       WRITE(98,*) ' ires1= ',ires1
       WRITE(98,*) ' ires1 must be integer greater or equal -1 !!'
       CLOSE(98)
     END SELECT
     Q12= Q1*Q1
     VzDxx(i)= VzDxx(i)+Q12/(Q12+Qz2) * EsLq1
     VzDxz(i)= VzDxz(i)+Q1*Qz/(Q12+Qz2) * EsLq1
     VzDzz(i)= VzDzz(i)+Qz2/(Q12+Qz2) * EsLq1
     Q2= Qxr2D(l,m,i,sigma)
     ires2= IQxr2D(l,m,i,sigma)
     SELECT CASE(ires2)
      CASE(-1)   ! Q<0
       EsLq2= 0.D0
      CASE(0)
       EsLq2= Iwave(1,i)
      CASE(1:)   ! ires1 >= 1
       IF (ires2<nqx) THEN
        Vzq(:)= Iwave(:,i)
        CALL Aitp1d2(nqx,VQx,Vzq,Q2,EsLq2,ires2)
       ELSE
        HBrL=0. ! The new effects are not relevant for large Q and are neglected
        HBrS=0.
        HGcL=0.
        HGcS=0.
       CALL Iwave_Init(Q2,Qz,EsLq2,"L",HBrL,HBrS,HGcL,HGcS)
       END IF
      CASE DEFAULT
       OPEN(98,FILE='Warning_Select_Case.wt')
       WRITE(98,*) ' ires2= ',ires2
       WRITE(98,*) ' ires2 must be integer greater or equal -2 !!'
       CLOSE(98)
     END SELECT
     Q22= Q2*Q2
     VzDxx(i)= VzDxx(i)+Q22/(Q22+Qz2) * EsLq2
     VzDxz(i)= VzDxz(i)+(-Q2*Qz)/(Q22+Qz2) * EsLq2
     VzDzz(i)= VzDzz(i)+Qz2/(Q22+Qz2) * EsLq2
    END DO
    CALL Simpson(VQz,VzDxx,nqz,Auxxx)
    CALL Simpson(VQz,VzDxz,nqz,Auxxz)
    CALL Simpson(VQz,VzDzz,nqz,Auxzz)
    Dxx(l,m)= Dxx(l,m)+(2./AbsUx)*Auxxx
    Dxz(l,m)= Dxz(l,m)+(2./AbsUx)*Auxxz
    Dzx(l,m)= Dzx(l,m)+(2./AbsUx)*Auxxz
    Dzz(l,m)= Dzz(l,m)+(2./AbsUx)*Auxzz
   END DO
  END IF
 END DO
END DO

RETURN
END SUBROUTINE Coef_D

SUBROUTINE Coef_Lwave(sigma,Dfdux,Dfduz,CoefA,CoefB)
USE Common_Params
USE Common_Arrays
USE Math_Constants
USE Phys_Constants
!USE hyp_2f1_module
IMPLICIT NONE
REAL*8 :: Ux,Qx,Qz
REAL*8 :: Qx2,Qz2,Q2,Q,Uresp,Uresm
REAL*8 :: Zlq,Zlqp,Zlqdif
REAL*8 :: Ztqp,Ztqdif
REAL*8 :: Aux,Aux0,Aux1
REAL*8 :: Aux2,Aux3,Aux4,Aux3m,Aux4m
REAL*8 :: AuxSE,D1,AuxCoef, D1Real, D1Aprox, D1Zerado
REAL*8 :: AuxScatElSpo
!REAL*8 :: AuxScatElInd
REAL*8 :: CoefEA,CoefLLSdA,CoefLLTdA,CoefLSTdA,CoefLTTdA,CoefLLLsA,CoefLLTsA
REAL*8 :: CoefEB,CoefLLSdB,CoefLLTdB,CoefLSTdB,CoefLTTdB,CoefLLLsB,CoefLLTsB
REAL*8 :: BremLA,GcollLB
REAL*8 :: Qxp,Qzp,Qxp2,Qzp2
REAL*8 :: Qxpaux,Qzpaux
REAL*8 :: Qxdif,Qzdif
REAL*8 :: HBrL,HBrS,HGcL,HGcS
REAL*8 :: Iqp,Iqdif,Iqpp,Iqpm
REAL*8 :: Fesum,Femax,Febkgr,Fef,Feb
REAL*8 :: Uz
REAL*8 :: Auxe,Auxe2
REAL*8 :: Iinit
REAL*8 :: Phip,Phi,Sphip,Cphip
REAL*8 :: Aalpha,AuxNum
REAL*8 :: Qstar,Qs
REAL*8 :: Beta
REAL*8, DIMENSION(nux) :: VintuxA,VauxuxA
REAL*8, DIMENSION(nux) :: VintuxB,VauxuxB
REAL*8, DIMENSION(nuz) :: VauxuzA,VintuzA
REAL*8, DIMENSION(nuz) :: VauxuzB,VintuzB
REAL*8, DIMENSION(nux,nuz), INTENT(in) :: Dfdux,Dfduz
REAL*8, DIMENSION(nqz) :: VauxqzA,VauxqzB
REAL*8, DIMENSION(nqx,nqz), INTENT(out) :: CoefA,CoefB
REAL*8, DIMENSION(nqx,nqz,-1:1) :: EcalL
REAL*8, DIMENSION(nqx,nqz,-1:1) :: EcalS
REAL*8, DIMENSION(nqx,nqz,-1:1) :: EcalT
REAL*8, DIMENSION(nqx) :: Vauxqxp!,Vauxqxdif
REAL*8, DIMENSION(nqx) :: VintqxA,VintqxB
REAL*8, DIMENSION(nqz) :: VintqzA,VintqzB
REAL*8, DIMENSION(nph) :: VintPhipA,VintPhipB
REAL*8, DIMENSION(nqx) :: VQxa
REAL*8, DIMENSION(nqz) :: VQza
INTEGER, INTENT(in) :: sigma
INTEGER :: i,k,l
INTEGER :: m,nqxa,nqza
INTEGER :: sigmap,sigmapp
INTEGER :: iqx,kqz,kqzp
INTEGER :: iresp,iresm
INTEGER :: ires,iresdif,iresdif2,kres
INTEGER :: S1,S2,S3,SS
INTEGER :: Asp,Aspp
CHARACTER(LEN=1) :: WaveT

m= 1   ! Auxiliary to the space profiles (just one point, for the moment)
Beta= VRTeTs(m)/VRNeNs(m)
nqxa= nqx
nqza= nqz
VQxa= VQx
VQza= VQz

IF(SpontEmis=="Yes") THEN
 AuxSE= 1.E0
ELSE
 AuxSE= 0.E0
END IF
IF(ScatElSpo=="Yes") THEN
 AuxScatElSpo= 1.E0
ELSE
 AuxScatElSpo= 0.E0
END IF
!IF(ScatElInd=="Yes") THEN
! AuxScatElInd= 1.E0
!ELSE
! AuxScatElInd= 0.E0
!END IF

CoefA= 0.
CoefB= 0.
DO iqx= 1,nqx
 Qx= VQx(iqx)
 DO kqz= 1,nqz
  Qz= VQz(kqz)
  Qx2= Qx**2
  Qz2= Qz**2
  Q2= Qx2+Qz2
  Q= SQRT(Q2)
  Zlq= ZL(Qx,Qz)
  Phi= ACOS(Qz/Q)
!TODO
  IF(Lemis=="Yes") THEN
   ! Contribution due to spontaneous and induced emission:
   IF (Qz>Qx) THEN
    DO l= 1,nux
     Ux= VUx(l)
     Uresp= UzrpLql(iqx,kqz,l,sigma)
     Uresm= UzrmLql(iqx,kqz,l,sigma)
     iresp= IUzrpLql(iqx,kqz,l,sigma)
     iresm= IUzrmLql(iqx,kqz,l,sigma)
     IF (iresp==0 .OR. iresp==nuz) THEN
      Uz= Uresp
      CALL Fe_Init(Ux,Uz,Fesum,Femax,Febkgr,Fef,Feb)
      VintuxA(l)= AuxSE * VRNeNs(m)*Geff*Fesum
      VintuxB(l)= -(sigma*Zlq) &
       * ( -Qx*2.*(Ux-U0x)*Femax/VRTeTs(m) &
       - Qx*2.*(Ux-U0x)*(Kappae+Alphae)/(Kappae*Uek2)*Febkgr &
       * (1+((Ux-U0x)**2+(Uz-U0z)**2)/(Kappae*Uek2))**(-1) &
       - Qx*2.*(Ux-Ufx)*Fef/VRTfperpTs(m) &
       - Qx*2.*(Ux-Ubx)*Feb/VRTbperpTs(m) &
       + Qz*2.*(Uz-U0z)*Femax/VRTeTs(m) &
       + Qz*2.*(Uz-U0z)*(Kappae+Alphae)/(Kappae*Uek2)*Febkgr &
       * (1+((Ux-U0x)**2+(Uz-U0z)**2)/(Kappae*Uek2))**(-1) &
       + Qz*2.*(Uz-Ufz)*Fef/VRTfparTs(m) &
       + Qz*2.*(Uz-Ubz)*Feb/VRTbparTs(m) )
     ELSE
      VauxuzA(:)= AuxSE * VRNeNs(m)*Geff*Fe(l,:)
      VauxuzB(:)= (sigma*Zlq)*(-Qx*Dfdux(l,:) &
        +Qz*Dfduz(l,:))
      CALL Aitp1d2(nuz,VUz,VauxuzA,Uresp,Aux,iresp)
      VintuxA(l)= Aux
      CALL Aitp1d2(nuz,VUz,VauxuzB,Uresp,Aux,iresp)
      VintuxB(l)= Aux
     END IF
     IF (iresm==0 .OR. iresm==nuz) THEN
      Uz= Uresm
      CALL Fe_Init(Ux,Uz,Fesum,Femax,Febkgr,Fef,Feb)
      VintuxA(l)= VintuxA(l) + AuxSE*VRNeNs(m)*Geff*Fesum
      VintuxB(l)= VintuxB(l)+ (sigma*Zlq) &
       * ( -Qx*2.*(Ux-U0x)*Femax/VRTeTs(m) &
       - Qx*2.*(Ux-U0x)*(Kappae+Alphae)/(Kappae*Uek2)*Febkgr &
       * (1+((Ux-U0x)**2+(Uz-U0z)**2)/(Kappae*Uek2))**(-1) &
       - Qx*2.*(Ux-Ufx)*Fef/VRTfperpTs(m) &
       - Qx*2.*(Ux-Ubx)*Feb/VRTbperpTs(m) &
       - Qz*2.*(Uz-U0z)*Femax/VRTeTs(m) &
       - Qz*2.*(Uz-U0z)*(Kappae+Alphae)/(Kappae*Uek2)*Febkgr &
       * (1+((Ux-U0x)**2+(Uz-U0z)**2)/(Kappae*Uek2))**(-1) &
       - Qz*2.*(Uz-Ufz)*Fef/VRTfparTs(m) &
       - Qz*2.*(Uz-Ubz)*Feb/VRTbparTs(m) )
     ELSE
      VauxuzA(:)= AuxSE * VRNeNs(m)*Geff*Fe(l,:)
      VauxuzB(:)= (sigma*Zlq)*(Qx*Dfdux(l,:) &
        +Qz*Dfduz(l,:))
      CALL Aitp1d2(nuz,VUz,VauxuzA,Uresm,Aux,iresm)
      VintuxA(l)= VintuxA(l) + Aux
      CALL Aitp1d2(nuz,VUz,VauxuzB,Uresm,Aux,iresm)
      VintuxB(l)= VintuxB(l) + Aux
     END IF
    END DO
    CALL Simpson(VUx,VintuxA,nux,Aux)
    CoefEA= (Pi/Q2/(ABS(Qz))) * VRNeNs(m)*Aux
    CALL Simpson(VUx,VintuxB,nux,Aux)
    CoefEB= (Pi/Q2/(ABS(Qz))) * VRNeNs(m)*Aux
   ELSE
    DO l= 1,nuz
     Uz= VUz(l)
     Uresp= UxrpLql(iqx,kqz,l,sigma)
     Uresm= UxrmLql(iqx,kqz,l,sigma)
     iresp= IUxrpLql(iqx,kqz,l,sigma)
     iresm= IUxrmLql(iqx,kqz,l,sigma)
     IF (iresm==0 .OR. iresm==nux) THEN
      IF (iresm==nux) THEN
       Ux= Uresm
       CALL Fe_Init(Ux,Uz,Fesum,Femax,Febkgr,Fef,Feb)
       VintuzA(l)= AuxSE * VRNeNs(m)*Geff*Fesum
       VintuzB(l)= -(sigma*Zlq) &
        * ( -Qx*2.*(Ux-U0x)*Femax/VRTeTs(m) &
        - Qx*2.*(Ux-U0x)*(Kappae+Alphae)/(Kappae*Uek2)*Febkgr &
        * (1+((Ux-U0x)**2+(Uz-U0z)**2)/(Kappae*Uek2))**(-1) &
        - Qx*2.*(Ux-Ufx)*Fef/VRTfperpTs(m) &
        - Qx*2.*(Ux-Ubx)*Feb/VRTbperpTs(m) &
        + Qz*2.*(Uz-U0z)*Femax/VRTeTs(m) &
        + Qz*2.*(Uz-U0z)*(Kappae+Alphae)/(Kappae*Uek2)*Febkgr &
        * (1+((Ux-U0x)**2+(Uz-U0z)**2)/(Kappae*Uek2))**(-1) &
        + Qz*2.*(Uz-Ufz)*Fef/VRTfparTs(m) &
        + Qz*2.*(Uz-Ubz)*Feb/VRTbparTs(m) )
      ELSE
       VintuzA(l)= 0.0
       VintuzB(l)= 0.0
      END IF
     ELSE
      VauxuxA(:)= AuxSE * VRNeNs(m)*Geff*Fe(:,l)
      VauxuxB(:)= (sigma*Zlq)*(-Qx*Dfdux(:,l) &
        +Qz*Dfduz(:,l))
      CALL Aitp1d2(nux,VUx,VauxuxA,Uresm,Aux,iresm)
      VintuzA(l)= Aux
      CALL Aitp1d2(nux,VUx,VauxuxB,Uresm,Aux,iresm)
      VintuzB(l)= Aux
     END IF
     IF (iresp==0 .OR. iresp==nux) THEN
      IF (iresp==nux) THEN
       Ux= Uresp
       CALL Fe_Init(Ux,Uz,Fesum,Femax,Febkgr,Fef,Feb)
       VintuzA(l)= VintuzA(l) + AuxSE * VRNeNs(m)*Geff*Fesum
       VintuzB(l)= VintuzB(l) + (sigma*Zlq) &
        * ( -Qx*2.*(Ux-U0x)*Femax/VRTeTs(m) &
        - Qx*2.*(Ux-U0x)*(Kappae+Alphae)/(Kappae*Uek2)*Febkgr &
        * (1+((Ux-U0x)**2+(Uz-U0z)**2)/(Kappae*Uek2))**(-1) &
        - Qx*2.*(Ux-Ufx)*Fef/VRTfperpTs(m) &
        - Qx*2.*(Ux-Ubx)*Feb/VRTbperpTs(m) &
        - Qz*2.*(Uz-U0z)*Femax/VRTeTs(m) &
        - Qz*2.*(Uz-U0z)*(Kappae+Alphae)/(Kappae*Uek2)*Febkgr &
        * (1+((Ux-U0x)**2+(Uz-U0z)**2)/(Kappae*Uek2))**(-1) &
        - Qz*2.*(Uz-Ufz)*Fef/VRTfparTs(m) &
        - Qz*2.*(Uz-Ubz)*Feb/VRTbparTs(m) )
      ELSE
       VintuzA(l)= VintuzA(l)+ 0.0
       VintuzB(l)= VintuzB(l)+ 0.0
      END IF
     ELSE
      VauxuxA(:)= AuxSE * VRNeNs(m)*Geff*Fe(:,l)
      VauxuxB(:)= (sigma*Zlq)*(Qx*Dfdux(:,l) &
        +Qz*Dfduz(:,l))
      CALL Aitp1d2(nux,VUx,VauxuxA,Uresp,Aux,iresp)
      VintuzA(l)= VintuzA(l) + Aux
      CALL Aitp1d2(nux,VUx,VauxuxB,Uresp,Aux,iresp)
      VintuzB(l)= VintuzB(l) + Aux
     END IF
    END DO
    CALL Simpson(VUz,VintuzA,nuz,Aux)
    CoefEA= (Pi/Q2/(ABS(Qx))) * VRNeNs(m)*Aux
    CALL Simpson(VUz,VintuzB,nuz,Aux)
    CoefEB= (Pi/Q2/(ABS(Qx))) * VRNeNs(m)*Aux
   END IF
  ELSE
   CoefEA= 0.
   CoefEB= 0.
  END IF

  IF(LdecayLS=="Yes") THEN
   ! Contribution due to spontaneous and induced decay involving L and S waves:
   EcalL(:,:,-1)= ILm(:,:)
   EcalL(:,:,+1)= ILp(:,:)
   EcalS(:,:,-1)= ISm(:,:)
   EcalS(:,:,+1)= ISp(:,:)
!   IF (Qz>Qx) THEN
!   IF (Qz<0.) THEN
!   ELSE
    VauxqzA= 0.   ! Initializes the integrand of the Qzp integrals
    VauxqzB= 0.
    DO kqzp= 1,nqz
     Qzp= VQz(kqzp)
     Qzp2= Qzp**2
     DO sigmap= -1,1,2
      Asp= sigmap/sigma
      Vauxqxp(:)= EcalL(:,kqzp,sigmap)
      DO sigmapp= -1,1,2
       Aspp= sigmapp/sigma
       DO S1= -1,1,2
        DO S2= -1,1,2
         Qxp= QxpLLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
         ires= IQxpLLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
         iresdif= IQxpdifLLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
         Qxp2= Qxp**2
         Qxdif= Qx-S1*Qxp
         Qzdif= Qz-S2*Qzp
         S3= INT(SIGN(1.D0,Qzdif))
         Zlqp= ZL(Qxp,Qzp)
         Zlqdif= ZL(Qxdif,ABS(Qzdif))

         D1Real= -S2*sigmap*(1.5*VRTeTs(m)/SQRT(VRNeNs(m))*Qxp &
          /SQRT(1.+1.5*Beta*(Qxp2+Qzp2))) &
          + S1*S3*sigmapp*AA*SQRT(VRTeTs(m))*Qxdif &
          /SQRT(Qxdif**2+Qzdif**2+EpsMin) &
          /(1.+Beta*(Qxdif**2+Qzdif**2)/2.)**(1.5)

         D1Aprox = -S2*sigmap*(1.5*VRTeTs(m)/SQRT(VRNeNs(m))*Qxp) &
          + S1*S3*sigmapp*AA*SQRT(VRTeTs(m))*Qxdif &
          /SQRT(Qxdif**2+Qzdif**2+EpsMin) &
          /(1.+Beta*(Qxdif**2+Qzdif**2)/2.)**(1.5)

         D1Zerado = -S2*sigmap*(1.5*VRTeTs(m)/SQRT(VRNeNs(m))*Qxp &
          /SQRT(1.+1.5*Beta*(Qxp2+Qzp2)))
         
         D1= D1Aprox

         AuxCoef= SQRT((Qxdif)**2+(Qzdif)**2)*(S1*Qx*Qxp+S2*Qz*Qzp)**2 &
          /(Qxp2+Qzp2+100*EpsMin)/(ABS(D1)+100*EpsMin)


         iresdif2= ABS(kqz-S2*kqzp)+S2
         !Vauxqxdif(:)= EcalS(:,ABS(iresdif2+S2*SIGN(1,iresdif2)),sigmapp)
         CALL Aux_Coef_Decay_z(nqxa,nqza,sigmapp,ires,iresdif,iresdif2,&
          Qxp,Qzp,Qxdif,Qzdif,"L","S",VQxa,VQza,Vauxqxp,EcalS,&
          Iqp,Iqdif)
         VauxqzA(kqzp)= VauxqzA(kqzp)+AuxCoef*sigma*Zlq*Iqp*Iqdif
         VauxqzB(kqzp)= VauxqzB(kqzp)-AuxCoef*(S2*sigmap*Zlqp*Iqdif &
           +S3*sigmapp*Zlqdif*Iqp)
        END DO
       END DO
      END DO
     END DO
    END DO
    CALL Simpson(VQz,VauxqzA,nqz,Aux)
    CoefLLSdA= Aux
    CALL Simpson(VQz,VauxqzB,nqz,Aux)
    CoefLLSdB= Aux
!   END IF
   CoefLLSdA= AA/VRNeNs(m)/SQRT(VRNeNs(m)*VRTeTs(m))*(sigma*Zlq/Q2) * CoefLLSdA
   CoefLLSdB= AA/VRNeNs(m)/SQRT(VRNeNs(m)*VRTeTs(m))*(sigma*Zlq/Q2) * CoefLLSdB
  ELSE
   CoefLLSdA= 0.
   CoefLLSdB= 0.
  END IF

  IF(LdecayLT=="Yes") THEN
   ! Contribution due to spontaneous and induced decay involving L and T waves:
   EcalL(:,:,-1)= ILm(:,:)
   EcalL(:,:,+1)= ILp(:,:)
   EcalT(:,:,-1)= ITm(:,:)
   EcalT(:,:,+1)= ITp(:,:)
!   IF (Qz>Qx) THEN
!   IF (Qz<0.) THEN
!   ELSE
    VauxqzA= 0.   ! Initializes the integrand of the Qzp integrals
    VauxqzB= 0.
    DO kqzp= 1,nqz
     Qzp= VQz(kqzp)
     Qzp2= Qzp**2
     DO sigmap= -1,1,2
      Asp= sigmap/sigma
      Vauxqxp(:)= EcalL(:,kqzp,sigmap)
      DO sigmapp= -1,1,2
       Aspp= sigmapp/sigma
       DO S1= -1,1,2
        DO S2= -1,1,2
         Qxp= QxpLLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
         ires= IQxpLLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
         iresdif= IQxpdifLLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
         Qxp2= Qxp**2
         Qxdif= Qx-S1*Qxp
         Qzdif= Qz-S2*Qzp
         S3= INT(SIGN(1.D0,Qzdif))
         Zlqp= ZL(Qxp,Qzp)
         Ztqdif= ZT(Qxdif,ABS(Qzdif))
         D1= -S2*sigmap*(1.5*Qxp/SQRT(1.+1.5*(Qxp2+Qzp2)) &
          +S1*S3*sigmapp*Qxdif/Ve2C2/SQRT(1.+(Qxdif**2+Qzdif**2)/Ve2C2))
         AuxCoef= (Q2*(Qxp2+Qzp2)-(S1*Qx*Qxp+S2*Qz*Qzp)**2) &
          /(Qxp2+Qzp2+EpsMin)/(Qxdif**2+Qzdif**2+EpsMin) &
          *((Qxp2+Qzp2)/S2/sigmap/Zlqp+Q2/sigma/Zlq)**2/(ABS(D1)+EpsMin)
         iresdif2= ABS(kqz-S2*kqzp)+S2
         !Vauxqxdif(:)= EcalT(:,ABS(iresdif2+S2*SIGN(1,iresdif2)),sigmapp)
         CALL Aux_Coef_Decay_z(nqxa,nqza,sigmapp,ires,iresdif,iresdif2,&
          Qxp,Qzp,Qxdif,Qzdif,"L","T",VQxa,VQza,Vauxqxp,EcalT,&
          Iqp,Iqdif)
         VauxqzA(kqzp)= VauxqzA(kqzp)+AuxCoef*sigma*Zlq*Iqp*Iqdif
         VauxqzB(kqzp)= VauxqzB(kqzp)-AuxCoef*(S2*sigmap*Zlqp*Iqdif &
           +S3*sigmapp*2.*Ztqdif*Iqp)
        END DO
       END DO
      END DO
     END DO
    END DO
    CALL Simpson(VQz,VauxqzA,nqz,Aux)
    CoefLLTdA= Aux
    CALL Simpson(VQz,VauxqzB,nqz,Aux)
    CoefLLTdB= Aux
!   END IF
   CoefLLTdA= (1./16.) * (sigma*Zlq/Q2) * CoefLLTdA
   CoefLLTdB= (1./16.) * (sigma*Zlq/Q2) * CoefLLTdB
  ELSE
   CoefLLTdA= 0.
   CoefLLTdB= 0.
  END IF

  IF(LdecayST=="Yes") THEN
   ! Contribution due to spontaneous and induced decay involving S and T waves:
   EcalS(:,:,-1)= ISm(:,:)
   EcalS(:,:,+1)= ISp(:,:)
   EcalT(:,:,-1)= ITm(:,:)
   EcalT(:,:,+1)= ITp(:,:)
!   IF (Qz>Qx) THEN
!   IF (Qz<0.) THEN
!   ELSE
    VauxqzA= 0.   ! Initializes the integrand of the Qzp integrals
    VauxqzB= 0.
    DO kqzp= 1,nqz
     Qzp= VQz(kqzp)
     Qzp2= Qzp**2
     DO sigmap= -1,1,2
      Asp= sigmap/sigma
      Vauxqxp(:)= EcalS(:,kqzp,sigmap)
      DO sigmapp= -1,1,2
       Aspp= sigmapp/sigma
       DO S1= -1,1,2
        DO S2= -1,1,2
         Qxp= QxpLSTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
         ires= IQxpLSTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
         iresdif= IQxpdifLSTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
         Qxp2= Qxp**2
         Qxdif= Qx-S1*Qxp
         Qzdif= Qz-S2*Qzp
         S3= INT(SIGN(1.D0,Qzdif))
         Zlqp= ZL(Qxp,Qzp)
         Ztqdif= ZT(Qxdif,ABS(Qzdif))
         D1= -S2*sigmap*AA*Qxp/SQRT(Qxp2+Qzp2+EpsMin) &
          /(1.+(Qxp2+Qzp2)/2.)**(1.5) &
          + S1*S3*sigmapp*Qxdif/Ve2C2/SQRT(1.+(Qxdif**2+Qzdif**2)/Ve2C2)
         AuxCoef= SQRT(Qxp2+Qzp2)*(Q2*(Qxp2+Qzp2)-(S1*Qx*Qxp+S2*Qz*Qzp)**2) &
          /(Qxdif**2+Qzdif**2+EpsMin)/(ABS(D1)+EpsMin)
         iresdif2= ABS(kqz-S2*kqzp)+S2
         !Vauxqxdif(:)= EcalT(:,ABS(iresdif2+S2*SIGN(1,iresdif2)),sigmapp)
         CALL Aux_Coef_Decay_z(nqxa,nqza,sigmapp,ires,iresdif,iresdif2,&
          Qxp,Qzp,Qxdif,Qzdif,"S","T",VQxa,VQza,Vauxqxp,EcalT,&
          Iqp,Iqdif)
         VauxqzA(kqzp)= VauxqzA(kqzp)+AuxCoef*sigma*Zlq*Iqp*Iqdif
         VauxqzB(kqzp)= VauxqzB(kqzp)-AuxCoef*(S2*sigmap*Zlqp*Iqdif &
           +S3*sigmapp*2.*Ztqdif*Iqp)
        END DO
       END DO
      END DO
     END DO
    END DO
    CALL Simpson(VQz,VauxqzA,nqz,Aux)
    CoefLSTdA= Aux
    CALL Simpson(VQz,VauxqzB,nqz,Aux)
    CoefLSTdB= Aux
!   END IF
   CoefLSTdA= (AA/2.) * (sigma*Zlq/Q2) * CoefLSTdA
   CoefLSTdB= (AA/2.) * (sigma*Zlq/Q2) * CoefLSTdB
  ELSE
   CoefLSTdA= 0.
   CoefLSTdB= 0.
  END IF

  IF(LdecayTT=="Yes") THEN
   ! Contribution due to spontaneous and induced decay involving two T waves:
   EcalT(:,:,-1)= ITm(:,:)
   EcalT(:,:,+1)= ITp(:,:)
!   IF (Qz>Qx) THEN
!   IF (Qz<0.) THEN
!   ELSE
    VauxqzA= 0.   ! Initializes the integrand of the Qzp integrals
    VauxqzB= 0.
    DO kqzp= 1,nqz
     Qzp= VQz(kqzp)
     Qzp2= Qzp**2
     DO sigmap= -1,1,2
      Asp= sigmap/sigma
      Vauxqxp(:)= EcalT(:,kqzp,sigmap)
      DO sigmapp= -1,1,2
       Aspp= sigmapp/sigma
       DO S1= -1,1,2
        DO S2= -1,1,2
         Qxp= QxpLTTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
         ires= IQxpLTTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
         iresdif= IQxpdifLTTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
         Qxp2= Qxp**2
         Qxdif= Qx-S1*Qxp
         Qzdif= Qz-S2*Qzp
         S3= INT(SIGN(1.D0,Qzdif))
         Ztqp= ZT(Qxp,Qzp)
         Ztqdif= ZT(Qxdif,ABS(Qzdif))
         D1= -S2*sigmap*Qxp/Ve2C2/SQRT(1.+(Qxp2+Qzp2)/Ve2C2) &
          + S1*S3*sigmapp*Qxdif/Ve2C2/SQRT(1.+(Qxdif**2+Qzdif**2)/Ve2C2)
         AuxCoef= (Q2+(S1*Qxp*Qxdif+S2*Qzp*Qzdif)**2 &
          /(Qxdif**2+Qzdif**2+EpsMin))/(ABS(D1)+EpsMin)/Ztqp**2/Ztqdif**2
         iresdif2= ABS(kqz-S2*kqzp)+S2
         !Vauxqxdif(:)= EcalT(:,ABS(iresdif2+S2*SIGN(1,iresdif2)),sigmapp)
         CALL Aux_Coef_Decay_z(nqxa,nqza,sigmapp,ires,iresdif,iresdif2,&
          Qxp,Qzp,Qxdif,Qzdif,"T","T",VQxa,VQza,Vauxqxp,EcalT,&
          Iqp,Iqdif)
         VauxqzA(kqzp)= VauxqzA(kqzp)+AuxCoef*sigma*Zlq*Iqp*Iqdif
         VauxqzB(kqzp)= VauxqzB(kqzp)-AuxCoef*(S2*sigmap*2.*Ztqp*Iqdif &
           +S3*sigmapp*2.*Ztqdif*Iqp)
        END DO
       END DO
      END DO
     END DO
    END DO
    CALL Simpson(VQz,VauxqzA,nqz,Aux)
    CoefLTTdA= Aux
    CALL Simpson(VQz,VauxqzB,nqz,Aux)
    CoefLTTdB= Aux
!   END IF
   CoefLTTdA= (sigma*Zlq/16.) * CoefLTTdA
   CoefLTTdB= (sigma*Zlq/16.) * CoefLTTdB
  ELSE
   CoefLTTdA= 0.
   CoefLTTdB= 0.
  END IF

  IF(LscatLL=="Yes") THEN
   ! Contribution due to spontaneous and induced scattering involving L waves:
   EcalL(:,:,-1)= ILm(:,:)
   EcalL(:,:,+1)= ILp(:,:)
   ! Electron contribution, Kappa:
   VintqxA= 0.
   VintqxB= 0.
   DO i= 1,nqx
    Qxp= VQx(i)
    VintqzA= 0.
    VintqzB= 0.
    DO k= 1,nqz
     Qzp= VQz(k)
     Zlqp= ZL(Qxp,Qzp)
     DO sigmap= -1,1,2
      DO S1= -1,1,2
       Qxdif= Qx-S1*Qxp
       S2= sigma*sigmap  ! The other value of S2 gives negligible contribution.
       Qzdif= Qz-S2*Qzp
       IF((Qxdif**2+Qzdif**2)<=1.E-6) THEN
        Auxe= 0.
        Auxe2= 0.
       ELSE
        Aux0= (S1*Qx*Qxp+S2*Qz*Qzp)**2/(Qxp**2+Qzp**2+EpsMin) &
         /SQRT(Qxdif**2+Qzdif**2)
        Aux= (sigma*Zlq-S2*sigmap*Zlqp)**2/(Qxdif**2+Qzdif**2)
        Aux1= Kappae*Uek2
        Aux2= 1.+Aux/Aux1
        Aux3m= Aux0 * (1.-RatioNke)*EXP(-Aux/VRTeTs(m))*SQRT(1./VRTeTs(m))
        Aux3= Aux0 * RatioNke/SQRT(Aux1)/Aux2**(Kappae+Alphae-0.5)
        Aux4m= Aux3m * 2./VRTeTs(m)
        Aux4= Aux3 * 2./Aux1/Aux2
        Auxe= Aux3m + Aux3 * RatioGammaem051 
        Auxe2= (Aux4m+Aux4*RatioGammaep051)*(sigma*Zlq-sigma*Zlqp)
       END IF
       VintqzA(k)= VintqzA(k) + (Auxe*AuxScatElSpo)*(sigma*Zlq) &
        * (AuxSE*Geff*(sigma*Zlq)*EcalL(i,k,sigmap))
       VintqzB(k)= VintqzB(k) &
        + (Auxe*AuxScatElSpo)*(AuxSE*Geff*(sigma*Zlq)*(-sigma*Zlqp)) &
        + Auxe2*(sigma*Zlq-sigma*Zlqp)*EcalL(i,k,sigmap)
      END DO
     END DO
    END DO
    CALL Simpson(VQz,VintqzA,nqz,Aux)
    VintqxA(i)= Aux
    CALL Simpson(VQz,VintqzB,nqz,Aux)
    VintqxB(i)= Aux
   END DO
   CALL Simpson(VQx,VintqxA,nqx,Aux)
   CoefLLLsA= 1./SQRT(Pi)/VRNeNs(m)/Q2 * Aux
   CALL Simpson(VQx,VintqxB,nqx,Aux)
   CoefLLLsB= 1./SQRT(Pi)/VRNeNs(m)/Q2 * Aux
   ! Ion contribution, spontaneous scattering, Maxw:
   Qstar= Q
   VintPhipA= 0.
   VintPhipB= 0.
   DO i= 1,nph
    Phip= VPhip(i)
    Qxp= Qstar*SIN(Phip)
    Qzp= Qstar*COS(Phip)
    IF( (Qxp>=VQx(nqx)) .OR. (Qzp>=VQz(nqz)) ) THEN
     HBrL=0. ! The new effects are not relevant for large Q and are neglected
     HBrS=0.
     HGcL=0.
     HGcS=0.
     WaveT= "L"
     CALL Iwave_Init(Qxp,Qzp,Iinit,WaveT,HBrL,HBrS,HGcL,HGcS)
     Iqpp= Iinit
     Iqpm= Iinit
    ELSE
     IF (Qxp<=VQx(1)) THEN
      ires= 1
      Qxpaux= VQx(1)+EpsMin
      IF (Qzp<=VQz(1)) THEN
       kres= 1
       Iqpp= ILp(1,1)
       Iqpm= ILm(1,1)
      ELSE
       kres= IQzLLL1(iqx,kqz,i)
       Qzpaux= Qzp
       CALL Aitp2d2(nqx,nqz,VQx,VQz,ILp,Qxpaux,Qzpaux,Iqpp,ires,kres)
       CALL Aitp2d2(nqx,nqz,VQx,VQz,ILm,Qxpaux,Qzpaux,Iqpm,ires,kres)
      END IF
     ELSE
      ires= IQxLLL1(iqx,kqz,i)
      Qxpaux= Qxp
      IF (Qzp<=VQz(1)) THEN
       kres= 1
       Qzpaux= VQz(1)+EpsMin
       CALL Aitp2d2(nqx,nqz,VQx,VQz,ILp,Qxpaux,Qzpaux,Iqpp,ires,kres)
       CALL Aitp2d2(nqx,nqz,VQx,VQz,ILm,Qxpaux,Qzpaux,Iqpm,ires,kres)
      ELSE
       kres= IQzLLL1(iqx,kqz,i)
       Qzpaux= Qzp
       CALL Aitp2d2(nqx,nqz,VQx,VQz,ILp,Qxpaux,Qzpaux,Iqpp,ires,kres)
       CALL Aitp2d2(nqx,nqz,VQx,VQz,ILm,Qxpaux,Qzpaux,Iqpm,ires,kres)
      END IF
     END IF
    END IF
    DO S1= -1,1,2
     sigmap= 1.
     Aux= (COS(Phi-S1*sigma*sigmap*Phip))**2
     VintPhipA(i)= VintPhipA(i)+Aux*Iqpp
     VintPhipB(i)= VintPhipB(i)+Aux
     sigmap= -1.
     Aux= (COS(Phi-S1*sigma*sigmap*Phip))**2
     VintPhipA(i)= VintPhipA(i)+Aux*Iqpm
     VintPhipB(i)= VintPhipB(i)+Aux
    END DO
   END DO
   CALL Simpson(VPhip,VintPhipA,nph,Aux)
   CoefLLLsA= CoefLLLsA + (2./3.)/VRNeNs(m) * (Zlq)**3 * AuxSE*Geff*Aux
   CALL Simpson(VPhip,VintPhipB,nph,Aux)
   CoefLLLsB= CoefLLLsB - (2./3.)/VRNeNs(m) * (Zlq)**3 * AuxSE*Geff*Aux
   ! Induced scattering due to the ions:
   VintPhipB= 0.
   DO i= 1,nph
    Phip= VPhip(i)
    Sphip= SIN(Phip)
    Cphip= COS(Phip)
    DO sigmap= -1,1,2
     DO S1= -1,1,2
      Aux1= (Q2+Qstar**2-2.*S1*Q*Qstar*COS(Phi-S1*sigma*sigmap*Phip))
      IF (Aux1<1.E-6) THEN
       AuxNum= SQRT(1.E-6)
       Aalpha= RMiMe/VRTiTs(m)*(9./4.)*(Qstar/Zlq)**2/1.E-6 
      ELSE
       AuxNum= SQRT(Aux1)
       Aalpha= RMiMe/VRTiTs(m)*(9./4.)*(Qstar/Zlq)**2/Aux1 
      END IF
      Aux= (COS(Phi-S1*sigma*sigmap*Phip))**2*AuxNum
      DO SS= -1,1,2
       Qs= Qstar+SS/SQRT(2.*Aalpha)
       Qxp= Qs*Sphip
       Qzp= Qs*Cphip
       IF( (Qxp>=VQx(nqx)) .OR. (Qzp>=VQz(nqz)) ) THEN
        HBrL=0. ! The new effects are not relevant for large Q and are neglected
        HBrS=0.
        HGcL=0.
        HGcS=0.
        WaveT= "L"
        CALL Iwave_Init(Qxp,Qzp,Iinit,WaveT,HBrL,HBrS,HGcL,HGcS)
        Iqp= Iinit
       ELSE
        S2= sigma*sigmap
        IF (Qxp<=VQx(1)) THEN
         ires= 1
         Qxpaux= VQx(1)+EpsMin
         IF (Qzp<=VQz(1)) THEN
          kres= 1
          IF (sigmap==1) THEN
           Iqp= ILp(1,1) 
          ELSE
           Iqp= ILm(1,1) 
          END IF
         ELSE
          kres= IQzLLL2(iqx,kqz,i,S2,S1,SS)
          Qzpaux= Qzp
          IF (sigmap==1) THEN
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ILp,Qxpaux,Qzpaux,Iqpp,ires,kres)
          ELSE
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ILm,Qxpaux,Qzpaux,Iqpm,ires,kres)
          END IF
         END IF
        ELSE
         ires= IQxLLL2(iqx,kqz,i,S2,S1,SS)
         Qxpaux= Qxp
         IF (Qzp<=VQz(1)) THEN
          kres= 1
          Qzpaux= VQz(1)+EpsMin
          IF (sigmap==1) THEN
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ILp,Qxpaux,Qzpaux,Iqpp,ires,kres)
          ELSE
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ILm,Qxpaux,Qzpaux,Iqpm,ires,kres)
          END IF
         ELSE
          kres= IQzLLL2(iqx,kqz,i,S2,S1,SS)
          Qzpaux= Qzp
          IF (sigmap==1) THEN
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ILp,Qxpaux,Qzpaux,Iqpp,ires,kres)
          ELSE
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ILm,Qxpaux,Qzpaux,Iqpm,ires,kres)
          END IF
         END IF
        END IF
       END IF
       VintPhipB(i)= VintPhipB(i)+Aux*Iqp * SS*Qs/Qstar
      END DO !SS
     END DO !S1
    END DO !sigmap
   END DO
   CALL Simpson(VPhip,VintPhipB,nph,Aux)
   CoefLLLsB= CoefLLLsB + (4./3.)/VRNeNs(m)/SQRT(VRTiTs(m)*RMiMe) &
      * EXP(-0.5)/SQRT(2.) * (Zlq)**2 * Aux
  ELSE
   CoefLLLsA= 0.
   CoefLLLsB= 0.
  END IF

  IF(LscatLT=="Yes") THEN
   ! Contribution due to spontaneous and induced scattering involving L and T 
   ! waves:
   EcalT(:,:,-1)= ITm(:,:)
   EcalT(:,:,+1)= ITp(:,:)
   ! Electron contribution, Kappa:
   VintqxA= 0.
   VintqxB= 0.
   DO i= 1,nqx
    Qxp= VQx(i)
    VintqzA= 0.
    VintqzB= 0.
    DO k= 1,nqz
     Qzp= VQz(k)
     Ztqp= ZT(Qxp,Qzp)
     DO sigmap= -1,1,2
      DO S1= -1,1,2
       Qxdif= Qx-S1*Qxp
       S2= sigma*sigmap  ! The other value of S2 gives negligible contribution.
       Qzdif= Qz-S2*Qzp
       IF((Qxdif**2+Qzdif**2)<=1.E-6) THEN
        Auxe= 0.
        Auxe2= 0.
       ELSE
        Aux0= (Q2*(Qxp**2+Qzp**2)-(S1*Qx*Qxp+S2*Qz*Qzp)**2) &
         /(Qxp**2+Qzp**2+EpsMin) &
         /SQRT(Qxdif**2+Qzdif**2)
        Aux= (sigma*Zlq-S2*sigmap*Ztqp)**2/(Qxdif**2+Qzdif**2)
        Aux1= Kappae*Uek2
        Aux2= 1.+Aux/Aux1
        Aux3= Aux0 * RatioNke/SQRT(Aux1)/Aux2**(Kappae+Alphae-0.5)
        Aux4= Aux3 * 2./Aux1/Aux2
        Auxe= Aux3 * RatioGammaem051 
        Auxe2= (Aux4m+Aux4*RatioGammaep051)*(sigma*Zlq-sigma*Ztqp)
       END IF
       VintqzA(k)= VintqzA(k) + (Auxe*AuxScatElSpo)*(sigma*Zlq) &
        * (AuxSE*Geff*(sigma*Zlq)*EcalT(i,k,sigmap)/2.)
       VintqzB(k)= VintqzB(k) &
        + (Auxe*AuxScatElSpo)*(AuxSE*Geff*(sigma*Zlq)*(-sigma*Ztqp)) &
        + Auxe2*(sigma*Zlq-sigma*Ztqp)*EcalT(i,k,sigmap)/2.
      END DO
     END DO
    END DO
    CALL Simpson(VQz,VintqzA,nqz,Aux)
    VintqxA(i)= Aux
    CALL Simpson(VQz,VintqzB,nqz,Aux)
    VintqxB(i)= Aux
   END DO
   CALL Simpson(VQx,VintqxA,nqx,Aux)
   CoefLLTsA= 1./SQRT(Pi)/VRNeNs(m)/Q2 * Aux
   CALL Simpson(VQx,VintqxB,nqx,Aux)
   CoefLLTsB= 1./SQRT(Pi)/VRNeNs(m)/Q2 * Aux
   ! Electron and ion contributions, Maxw:
   Qstar= SQRT(3./2.*Ve2C2)*Q
   VintPhipA= 0.
   VintPhipB= 0.
   DO i= 1,nph
    Phip= VPhip(i)
    Qxp= Qstar*SIN(Phip)
    Qzp= Qstar*COS(Phip)
    IF( (Qxp>=VQx(nqx)) .OR. (Qzp>=VQz(nqz)) ) THEN
     HBrL=0. ! The new effects are not relevant for large Q and are neglected
     HBrS=0.
     HGcL=0.
     HGcS=0.
     WaveT= "T"
     CALL Iwave_Init(Qxp,Qzp,Iinit,WaveT,HBrL,HBrS,HGcL,HGcS)
     Iqpp= Iinit
     Iqpm= Iinit
    ELSE
     IF (Qxp<=VQx(1)) THEN
      ires= 1
      Qxpaux= VQx(1)+EpsMin
      IF (Qzp<=VQz(1)) THEN
       kres= 1
       Iqpp= ITp(1,1)
       Iqpm= ITm(1,1)
      ELSE
       kres= IQzLLT1(iqx,kqz,i)
       Qzpaux= Qzp
       CALL Aitp2d2(nqx,nqz,VQx,VQz,ITp,Qxpaux,Qzpaux,Iqpp,ires,kres)
       CALL Aitp2d2(nqx,nqz,VQx,VQz,ITm,Qxpaux,Qzpaux,Iqpm,ires,kres)
      END IF
     ELSE
      ires= IQxLLT1(iqx,kqz,i)
      Qxpaux= Qxp
      IF (Qzp<=VQz(1)) THEN
       kres= 1
       Qzpaux= VQz(1)+EpsMin
       CALL Aitp2d2(nqx,nqz,VQx,VQz,ITp,Qxpaux,Qzpaux,Iqpp,ires,kres)
       CALL Aitp2d2(nqx,nqz,VQx,VQz,ITm,Qxpaux,Qzpaux,Iqpm,ires,kres)
      ELSE
       kres= IQzLLT1(iqx,kqz,i)
       Qzpaux= Qzp
       CALL Aitp2d2(nqx,nqz,VQx,VQz,ITp,Qxpaux,Qzpaux,Iqpp,ires,kres)
       CALL Aitp2d2(nqx,nqz,VQx,VQz,ITm,Qxpaux,Qzpaux,Iqpm,ires,kres)
      END IF
     END IF
    END IF
    DO S1= -1,1,2
     sigmap= 1.
     Aux= (SIN(Phi-S1*sigma*sigmap*Phip))**2
     VintPhipA(i)= VintPhipA(i)+Aux*Iqpp/2.
     VintPhipB(i)= VintPhipB(i)+Aux
     sigmap= -1.
     Aux= (SIN(Phip))**2
     VintPhipA(i)= VintPhipA(i)+Aux*Iqpm/2.
     VintPhipB(i)= VintPhipB(i)+Aux
    END DO
   END DO
   CALL Simpson(VPhip,VintPhipA,nph,Aux)
   CoefLLTsA= CoefLLTsA + 2.*Ve2C2/VRNeNs(m) * (Zlq)**3 * AuxSE*Geff*Aux &
    * (2.-RatioNke)
   CALL Simpson(VPhip,VintPhipB,nph,Aux)
   CoefLLTsB= CoefLLTsB + 2.*Ve2C2/VRNeNs(m) * (Zlq)**3 * AuxSE*Geff*Aux &
    * (2.-RatioNke)
   ! Induced scattering due to the ions:
   VintPhipB= 0.
   DO i= 1,nph
    Phip= VPhip(i)
    Sphip= SIN(Phip)
    Cphip= COS(Phip)
    DO sigmap= -1,1,2
     DO S1= -1,1,2
      Aux1= (Q2+Qstar**2-2.*S1*Q*Qstar*COS(Phi-S1*sigma*sigmap*Phip))
      IF (Aux1<1.E-6) THEN
       AuxNum= SQRT(1.E-6)
       Aalpha= RMiMe/VRTiTs(m)*(1./Ve2C2**2)*(Qstar/Zlq)**2/1.E-6
      ELSE
       AuxNum= SQRT(Aux1)
       Aalpha= RMiMe/VRTiTs(m)*(1./Ve2C2**2)*(Qstar/Zlq)**2/Aux1 
      END IF
      Aux= (SIN(Phi-S1*sigma*sigmap*Phip))**2*AuxNum
      DO SS= -1,1,2
       Qs= Qstar+SS/SQRT(2.*Aalpha)
       Qxp= Qs*Sphip
       Qzp= Qs*Cphip
       IF( (Qxp>=VQx(nqx)) .OR. (Qzp>=VQz(nqz)) ) THEN
        HBrL=0. ! The new effects are not relevant for large Q and are neglected
        HBrS=0.
        HGcL=0.
        HGcS=0.
        WaveT= "T"
        CALL Iwave_Init(Qxp,Qzp,Iinit,WaveT,HBrL,HBrS,HGcL,HGcS)
        Iqp= Iinit
       ELSE
        S2= sigma*sigmap
        IF (Qxp<=VQx(1)) THEN
         ires= 1
         Qxpaux= VQx(1)+EpsMin
         IF (Qzp<=VQz(1)) THEN
          kres= 1
          IF (sigmap==1) THEN
           Iqp= ITp(1,1) 
          ELSE
           Iqp= ITm(1,1) 
          END IF
         ELSE
          kres= IQzLLT2i(iqx,kqz,i,S2,S1,SS)
          Qzpaux= Qzp
          IF (sigmap==1) THEN
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ITp,Qxpaux,Qzpaux,Iqpp,ires,kres)
          ELSE
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ITm,Qxpaux,Qzpaux,Iqpm,ires,kres)
          END IF
         END IF
        ELSE
         ires= IQxLLT2i(iqx,kqz,i,S2,S1,SS)
         Qxpaux= Qxp
         IF (Qzp<=VQz(1)) THEN
          kres= 1
          Qzpaux= VQz(1)+EpsMin
          IF (sigmap==1) THEN
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ITp,Qxpaux,Qzpaux,Iqpp,ires,kres)
          ELSE
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ITm,Qxpaux,Qzpaux,Iqpm,ires,kres)
          END IF
         ELSE
          kres= IQzLLT2i(iqx,kqz,i,S2,S1,SS)
          Qzpaux= Qzp
          IF (sigmap==1) THEN
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ITp,Qxpaux,Qzpaux,Iqpp,ires,kres)
          ELSE
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ITm,Qxpaux,Qzpaux,Iqpm,ires,kres)
          END IF
         END IF
        END IF
       END IF
       VintPhipB(i)= VintPhipB(i)+Aux*Iqp/2. * SS*Qs/Qstar
      END DO !SS
     END DO !S1
    END DO !sigmap
   END DO
   CALL Simpson(VPhip,VintPhipB,nph,Aux)
   CoefLLTsB= CoefLLTsB + 2.*Ve2C2/VRNeNs(m)/SQRT(VRTiTs(m)*RMiMe) &
      * EXP(-0.5)/SQRT(2.) * (Zlq)**2 * Aux * (1.-RatioNki)
   ! Induced scattering due to the electrons:
   VintPhipB= 0.
   DO i= 1,nph
    Phip= VPhip(i)
    Sphip= SIN(Phip)
    Cphip= COS(Phip)
    DO sigmap= -1,1,2
     DO S1= -1,1,2
      Aux1= (Q2+Qstar**2-2.*S1*Q*Qstar*COS(Phi-S1*sigma*sigmap*Phip))
      IF (Aux1<1.E-6) THEN
       AuxNum= SQRT(1.E-6)
       Aalpha= 1./VRTeTs(m)*(1./Ve2C2**2)*(Qstar/Zlq)**2/1.E-6
      ELSE
       AuxNum= SQRT(Aux1)
       Aalpha= 1./VRTeTs(m)*(1./Ve2C2**2)*(Qstar/Zlq)**2/Aux1 
      END IF
      Aux= (SIN(Phi-S1*sigma*sigmap*Phip))**2*AuxNum
      DO SS= -1,1,2
       Qs= Qstar+SS/SQRT(2.*Aalpha)
       Qxp= Qs*Sphip
       Qzp= Qs*Cphip
       Ztqp= ZT(Qxp,Qzp)
       IF( (Qxp>=VQx(nqx)) .OR. (Qzp>=VQz(nqz)) ) THEN
        HBrL=0. ! The new effects are not relevant for large Q and are neglected
        HBrS=0.
        HGcL=0.
        HGcS=0.
        WaveT= "T"
        CALL Iwave_Init(Qxp,Qzp,Iinit,WaveT,HBrL,HBrS,HGcL,HGcS)
        Iqp= Iinit
       ELSE
        S2= sigma*sigmap
        IF (Qxp<=VQx(1)) THEN
         ires= 1
         Qxpaux= VQx(1)+EpsMin
         IF (Qzp<=VQz(1)) THEN
          kres= 1
          IF (sigmap==1) THEN
           Iqp= ITp(1,1) 
          ELSE
           Iqp= ITm(1,1) 
          END IF
         ELSE
          kres= IQzLLT2e(iqx,kqz,i,S2,S1,SS)
          Qzpaux= Qzp
          IF (sigmap==1) THEN
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ITp,Qxpaux,Qzpaux,Iqpp,ires,kres)
          ELSE
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ITm,Qxpaux,Qzpaux,Iqpm,ires,kres)
          END IF
         END IF
        ELSE
         ires= IQxLLT2e(iqx,kqz,i,S2,S1,SS)
         Qxpaux= Qxp
         IF (Qzp<=VQz(1)) THEN
          kres= 1
          Qzpaux= VQz(1)+EpsMin
          IF (sigmap==1) THEN
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ITp,Qxpaux,Qzpaux,Iqpp,ires,kres)
          ELSE
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ITm,Qxpaux,Qzpaux,Iqpm,ires,kres)
          END IF
         ELSE
          kres= IQzLLT2e(iqx,kqz,i,S2,S1,SS)
          Qzpaux= Qzp
          IF (sigmap==1) THEN
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ITp,Qxpaux,Qzpaux,Iqpp,ires,kres)
          ELSE
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ITm,Qxpaux,Qzpaux,Iqpm,ires,kres)
          END IF
         END IF
        END IF
       END IF
       VintPhipB(i)= VintPhipB(i)+Aux*Iqp/2. * SS*Qs/Qstar * (Zlq-Ztqp)
      END DO !SS
     END DO !S1
    END DO !sigmap
   END DO
   CALL Simpson(VPhip,VintPhipB,nph,Aux)
   CoefLLTsB= CoefLLTsB - 4.*Ve2C2/VRNeNs(m)/SQRT(VRTeTs(m)) &
      * EXP(-0.5)/SQRT(2.) * (Zlq) * Aux * (1.-RatioNke)
  ELSE
   CoefLLTsA= 0.
   CoefLLTsB= 0.
  END IF

  IF (NewEffects2== "Yes") THEN
   IF(Gcoll== "Yes")THEN
    GcollLB= 2.*GcollLp(iqx,kqz)
   ELSE
    GcollLB=0.
   END IF
   IF(Bremss== "Yes")THEN
    BremLA= BremL(iqx,kqz)
   ELSE
    BremLA= 0.
   END IF
  ELSE
   GcollLB=0.
   BremLA= 0.
  END IF

  CoefA(iqx,kqz)= CoefEA+CoefLLSdA+CoefLLTdA+CoefLSTdA+CoefLTTdA &
     +CoefLLLsA+CoefLLTsA+BremLA
  CoefB(iqx,kqz)= CoefEB+CoefLLSdB+CoefLLTdB+CoefLSTdB+CoefLTTdB &
     +CoefLLLsB+CoefLLTsB+GcollLB
 END DO
END DO

! CALL Output_Coef("Lwave",CoefA,CoefB)

RETURN
END SUBROUTINE Coef_Lwave

SUBROUTINE Coef_Swave(sigma,Dfdux,Dfduz,CoefA,CoefB)
USE Common_Params
USE Common_Arrays
USE Math_Constants
USE Phys_Constants
!USE hyp_2f1_module
IMPLICIT NONE
REAL*8 :: Ux,Qx,Qz,Qx2,Qz2,Q,Q2
REAL*8 :: Zlq,Zsq,Zsoq,Uresp,Uresm
REAL*8 :: Aux,AuxSE,D1,AuxCoef
REAL*8 :: CoefEA,CoefSLLdA,CoefSLTdA
REAL*8 :: CoefEB,CoefSLLdB,CoefSLTdB
REAL*8 :: BremSA,GcollSB
REAL*8 :: Zlqp,Zlqdif,Ztqdif
REAL*8 :: Muq
REAL*8 :: Qxp,Qzp
REAL*8 :: Qxp2,Qzp2
REAL*8 :: Qxdif,Qzdif
REAL*8 :: Iqp,Iqdif
REAL*8 :: Fesum,Femax,Febkgr,Fef,Feb
REAL*8 :: Uz
REAL*8 :: Beta
REAL*8, DIMENSION(nux) :: VintuxA,VauxuxA
REAL*8, DIMENSION(nux) :: VintuxB,VauxuxB
REAL*8, DIMENSION(nuz) :: VauxuzA,VintuzA
REAL*8, DIMENSION(nuz) :: VauxuzB,VintuzB
REAL*8, DIMENSION(nux,nuz), INTENT(in) :: Dfdux,Dfduz
REAL*8, DIMENSION(nqz) :: VauxqzA,VauxqzB
REAL*8, DIMENSION(nqx) :: Vauxqxp!,Vauxqxdif
REAL*8, DIMENSION(nqx,nqz), INTENT(out) :: CoefA,CoefB
REAL*8, DIMENSION(nqx,nqz,-1:1) :: EcalL
REAL*8, DIMENSION(nqx,nqz,-1:1) :: EcalT
REAL*8, DIMENSION(nqx) :: VQxa
REAL*8, DIMENSION(nqz) :: VQza
INTEGER, INTENT(in) :: sigma
INTEGER :: l,sigmap,sigmapp
INTEGER :: iqx,kqz,kqzp
INTEGER :: iresp,iresm
INTEGER :: ires,iresdif,iresdif2
INTEGER :: m,nqxa,nqza
INTEGER :: S1,S2,S3
INTEGER :: Asp,Aspp

m= 1   ! Auxiliary to the space profiles (just one point, for the moment)
Beta= VRTeTs(m)/VRNeNs(m)
nqxa= nqx
nqza= nqz
VQxa= VQx
VQza= VQz

IF(SpontEmis=="Yes") THEN
 AuxSE= 1.E0
ELSE
 AuxSE= 0.E0
END IF

CoefA= 0.
CoefB= 0.
DO iqx= 1,nqx
 Qx= VQx(iqx)
 DO kqz= 1,nqz
  Qz= VQz(kqz)
  Qx2= Qx**2
  Qz2= Qz**2
  Q2= Qx2+Qz2
  Q= SQRT(Q2)
  Muq= Q**3*AA/2.
  Zlq= ZL(Qx,Qz)
  Zsq= ZS(Qx,Qz)
  Zsoq= AA*SQRT(VRTeTs(m))/SQRT(1.+Q2/2.*VRTeTs(m)/VRNeNs(m))
!TODO
  IF(Semis=="Yes") THEN
   ! Contribution due to spontaneous and induced emission:
   IF (Qz>Qx) THEN
    DO l= 1,nux
     Ux= VUx(l)
     Uresp= UzrpSql(iqx,kqz,l,sigma)
     Uresm= UzrmSql(iqx,kqz,l,sigma)
     iresp= IUzrpSql(iqx,kqz,l,sigma)
     iresm= IUzrmSql(iqx,kqz,l,sigma)
     IF (iresp==0 .OR. iresp==nuz) THEN
      Uz= Uresp
      CALL Fe_Init(Ux,Uz,Fesum,Femax,Febkgr,Fef,Feb)
      VintuxA(l)= AuxSE * VRNeNs(m)*Geff*Fesum
      VintuxB(l)= -(sigma*Zlq) &
       * ( -Qx*2.*(Ux-U0x)*Femax/VRTeTs(m) &
       - Qx*2.*(Ux-U0x)*(Kappae+Alphae)/(Kappae*Uek2)*Febkgr &
       * (1+((Ux-U0x)**2+(Uz-U0z)**2)/(Kappae*Uek2))**(-1) &
       - Qx*2.*(Ux-Ufx)*Fef/VRTfperpTs(m) &
       - Qx*2.*(Ux-Ubx)*Feb/VRTbperpTs(m) &
       + Qz*2.*(Uz-U0z)*Femax/VRTeTs(m) &
       + Qz*2.*(Uz-U0z)*(Kappae+Alphae)/(Kappae*Uek2)*Febkgr &
       * (1+((Ux-U0x)**2+(Uz-U0z)**2)/(Kappae*Uek2))**(-1) &
       + Qz*2.*(Uz-Ufz)*Fef/VRTfparTs(m) &
       + Qz*2.*(Uz-Ubz)*Feb/VRTbparTs(m) )
     ELSE
      VauxuzA(:)= AuxSE * VRNeNs(m)*Geff*Fe(l,:)
      VauxuzB(:)= (sigma*Zlq)*(-Qx*Dfdux(l,:) &
        +Qz*Dfduz(l,:))
      CALL Aitp1d2(nuz,VUz,VauxuzA,Uresp,Aux,iresp)
      VintuxA(l)= Aux
      CALL Aitp1d2(nuz,VUz,VauxuzB,Uresp,Aux,iresp)
      VintuxB(l)= Aux
     END IF
     IF (iresm==0 .OR. iresm==nuz) THEN
      Uz= Uresm
      CALL Fe_Init(Ux,Uz,Fesum,Femax,Febkgr,Fef,Feb)
      VintuxA(l)= VintuxA(l) + AuxSE*VRNeNs(m)*Geff*Fesum
      VintuxB(l)= VintuxB(l)+ (sigma*Zlq) &
       * ( -Qx*2.*(Ux-U0x)*Femax/VRTeTs(m) &
       - Qx*2.*(Ux-U0x)*(Kappae+Alphae)/(Kappae*Uek2)*Febkgr &
       * (1+((Ux-U0x)**2+(Uz-U0z)**2)/(Kappae*Uek2))**(-1) &
       - Qx*2.*(Ux-Ufx)*Fef/VRTfperpTs(m) &
       - Qx*2.*(Ux-Ubx)*Feb/VRTbperpTs(m) &
       - Qz*2.*(Uz-U0z)*Femax/VRTeTs(m) &
       - Qz*2.*(Uz-U0z)*(Kappae+Alphae)/(Kappae*Uek2)*Febkgr &
       * (1+((Ux-U0x)**2+(Uz-U0z)**2)/(Kappae*Uek2))**(-1) &
       - Qz*2.*(Uz-Ufz)*Fef/VRTfparTs(m) &
       - Qz*2.*(Uz-Ubz)*Feb/VRTbparTs(m) )
     ELSE
      VauxuzA(:)= AuxSE * VRNeNs(m)*Geff*Fe(l,:)
      VauxuzB(:)= (sigma*Zlq)*(Qx*Dfdux(l,:) &
        +Qz*Dfduz(l,:))
      CALL Aitp1d2(nuz,VUz,VauxuzA,Uresm,Aux,iresm)
      VintuxA(l)= VintuxA(l) + Aux
      CALL Aitp1d2(nuz,VUz,VauxuzB,Uresm,Aux,iresm)
      VintuxB(l)= VintuxB(l) + Aux
     END IF
    END DO
    CALL Simpson(VUx,VintuxA,nux,CoefEA)
    CALL Simpson(VUx,VintuxB,nux,CoefEB)
    CoefEA= CoefEA/(ABS(Qz))
    CoefEB= CoefEB/(ABS(Qz))
   ELSE
    DO l= 1,nuz
     Uz= VUz(l)
     Uresp= UxrpSql(iqx,kqz,l,sigma)
     Uresm= UxrmSql(iqx,kqz,l,sigma)
     iresp= IUxrpSql(iqx,kqz,l,sigma)
     iresm= IUxrmSql(iqx,kqz,l,sigma)
     IF (iresm==0 .OR. iresm==nux) THEN
      IF (iresm==nux) THEN
       Ux= Uresm
       CALL Fe_Init(Ux,Uz,Fesum,Femax,Febkgr,Fef,Feb)
       VintuzA(l)= AuxSE * VRNeNs(m)*Geff*Fesum
       VintuzB(l)= -(sigma*Zlq) &
        * ( -Qx*2.*(Ux-U0x)*Femax/VRTeTs(m) &
        - Qx*2.*(Ux-U0x)*(Kappae+Alphae)/(Kappae*Uek2)*Febkgr &
        * (1+((Ux-U0x)**2+(Uz-U0z)**2)/(Kappae*Uek2))**(-1) &
        - Qx*2.*(Ux-Ufx)*Fef/VRTfperpTs(m) &
        - Qx*2.*(Ux-Ubx)*Feb/VRTbperpTs(m) &
        + Qz*2.*(Uz-U0z)*Femax/VRTeTs(m) &
        + Qz*2.*(Uz-U0z)*(Kappae+Alphae)/(Kappae*Uek2)*Febkgr &
        * (1+((Ux-U0x)**2+(Uz-U0z)**2)/(Kappae*Uek2))**(-1) &
        + Qz*2.*(Uz-Ufz)*Fef/VRTfparTs(m) &
        + Qz*2.*(Uz-Ubz)*Feb/VRTbparTs(m) )
      ELSE
       VintuzA(l)= 0.0
       VintuzB(l)= 0.0
      END IF
     ELSE
      VauxuxA(:)= AuxSE * VRNeNs(m)*Geff*Fe(:,l)
      VauxuxB(:)= (sigma*Zlq)*(-Qx*Dfdux(:,l) &
        +Qz*Dfduz(:,l))
      CALL Aitp1d2(nux,VUx,VauxuxA,Uresm,Aux,iresm)
      VintuzA(l)= Aux
      CALL Aitp1d2(nux,VUx,VauxuxB,Uresm,Aux,iresm)
      VintuzB(l)= Aux
     END IF
     IF (iresp==0 .OR. iresp==nux) THEN
      IF (iresp==nux) THEN
       Ux= Uresp
       CALL Fe_Init(Ux,Uz,Fesum,Femax,Febkgr,Fef,Feb)
       VintuzA(l)= VintuzA(l) + AuxSE * VRNeNs(m)*Geff*Fesum
       VintuzB(l)= -(sigma*Zlq) &
        * ( -Qx*2.*(Ux-U0x)*Femax/VRTeTs(m) &
        - Qx*2.*(Ux-U0x)*(Kappae+Alphae)/(Kappae*Uek2)*Febkgr &
        * (1+((Ux-U0x)**2+(Uz-U0z)**2)/(Kappae*Uek2))**(-1) &
        - Qx*2.*(Ux-Ufx)*Fef/VRTfperpTs(m) &
        - Qx*2.*(Ux-Ubx)*Feb/VRTbperpTs(m) &
        + Qz*2.*(Uz-U0z)*Femax/VRTeTs(m) &
        + Qz*2.*(Uz-U0z)*(Kappae+Alphae)/(Kappae*Uek2)*Febkgr &
        * (1+((Ux-U0x)**2+(Uz-U0z)**2)/(Kappae*Uek2))**(-1) &
        + Qz*2.*(Uz-Ufz)*Fef/VRTfparTs(m) &
        + Qz*2.*(Uz-Ubz)*Feb/VRTbparTs(m) )
      ELSE
       VintuzA(l)= VintuzA(l)+ 0.0
       VintuzB(l)= VintuzB(l)+ 0.0
      END IF
     ELSE
      VauxuxA(:)= AuxSE * VRNeNs(m)*Geff*Fe(:,l)
      VauxuxB(:)= (sigma*Zlq)*(Qx*Dfdux(:,l) &
        +Qz*Dfduz(:,l))
      CALL Aitp1d2(nux,VUx,VauxuxA,Uresp,Aux,iresp)
      VintuzA(l)= VintuzA(l) + Aux
      CALL Aitp1d2(nux,VUx,VauxuxB,Uresp,Aux,iresp)
      VintuzB(l)= VintuzB(l) + Aux
     END IF
    END DO
    CALL Simpson(VUz,VintuzA,nuz,CoefEA)
    CALL Simpson(VUz,VintuzB,nuz,CoefEB)
    CoefEA= CoefEA/(ABS(Qx))
    CoefEB= CoefEB/(ABS(Qx))
   END IF
   CoefEA= Pi * Q*AA/2. &
        * (VRTeTs(m))**2/SQRT(VRTeTs(m)*VRNeNs(m)) * CoefEA
   CoefEB= Pi * Q*AA/2. &
        * (VRTeTs(m))**2/SQRT(VRTeTs(m)*VRNeNs(m)) * CoefEB
   ! Adding the ion contribution
   Aux= SQRT(Pi)*AA/2. * VRTeTs(m)/SQRT(VRNeNs(m)) &
       * SQRT(RMiMe*RTeTi)*EXP(-RMiMe/VRTiTs(m)*(Zsoq)**2)
   CoefEA= CoefEA + Aux*(AuxSE * VRNeNs(m)*Geff)
   CoefEB= CoefEB + Aux*(- 2./VRTiTs(m)*Zlq*Zsq)
  ELSE
   CoefEA= 0.
   CoefEB= 0.
  END IF

  IF(SdecayLL=="Yes") THEN
   ! Contribution due to spontaneous and induced decay:
   EcalL(:,:,-1)= ILm(:,:)
   EcalL(:,:,+1)= ILp(:,:)
!   IF (Qz>Qx) THEN
!   IF (Qz<0.) THEN
!   ELSE
    VauxqzA= 0.   ! Initializes the integrand of the Qzp integrals
    VauxqzB= 0.
    DO kqzp= 1,nqz
     Qzp= VQz(kqzp)
     Qzp2= Qzp**2
     DO sigmap= -1,1,2
      Asp= sigmap/sigma
      Vauxqxp(:)= EcalL(:,kqzp,sigmap)
      DO sigmapp= -1,1,2
       Aspp= sigmapp/sigma
       DO S1= -1,1,2
        DO S2= -1,1,2
         Qxp= QxpSLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
         ires= IQxpSLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
         iresdif= IQxpdifSLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
         Qxp2= Qxp**2
         Qxdif= Qx-S1*Qxp
         Qzdif= Qz-S2*Qzp
         S3= INT(SIGN(1.D0,Qzdif))
         Zlqp= ZL(Qxp,Qzp)
         Zlqdif= ZL(Qxdif,ABS(Qzdif))
         D1= -S2*sigmap*1.5*Beta*Qxp &
          /SQRT(1.+1.5*Beta*(Qxp2+Qzp2)) &
          + S1*S3*sigmapp*1.5*Beta*Qxdif &
          /SQRT(1.+1.5*Beta*(Qxdif**2+Qzdif**2))
         AuxCoef= (Qxp*Qxdif+Qzp*Qzdif)**2 &
          /(Qxp2+Qzp2+EpsMin)/(Qxdif**2+Qzdif**2+EpsMin)/(ABS(D1)+EpsMin)
         iresdif2= ABS(kqz-S2*kqzp)+S2
         !Vauxqxdif(:)= EcalL(:,ABS(iresdif2+S2*SIGN(1,iresdif2)),sigmapp)
         CALL Aux_Coef_Decay_z(nqxa,nqza,sigmapp,ires,iresdif,iresdif2,&
          Qxp,Qzp,Qxdif,Qzdif,"L","L",VQxa,VQza,Vauxqxp,EcalL,&
          Iqp,Iqdif)
         VauxqzA(kqzp)= VauxqzA(kqzp)+AuxCoef*sigma*Zlq*Iqp*Iqdif
         VauxqzB(kqzp)= VauxqzB(kqzp)-AuxCoef*(S2*sigmap*Zlqp*Iqdif &
          +S3*sigmapp*Zlqdif*Iqp)
        END DO
       END DO
      END DO
     END DO
    END DO
    CALL Simpson(VQz,VauxqzA,nqz,Aux)
    CoefSLLdA= Aux
    CALL Simpson(VQz,VauxqzB,nqz,Aux)
    CoefSLLdB= Aux
!   END IF
   CoefSLLdA= (AA/2.)*(sigma*Zlq)*Q &
       /VRNeNs(m)/SQRT(VRNeNs(m)*VRTeTs(m)) * CoefSLLdA
   CoefSLLdB= (AA/2.)*(sigma*Zlq)*Q &
       /VRNeNs(m)/SQRT(VRNeNs(m)*VRTeTs(m)) * CoefSLLdB
  ELSE
   CoefSLLdA= 0.
   CoefSLLdB= 0.
  END IF

  IF(SdecayLT=="Yes") THEN
   ! Contribution due to spontaneous and induced decay involving L and T waves:
   EcalL(:,:,-1)= ILm(:,:)
   EcalL(:,:,+1)= ILp(:,:)
   EcalT(:,:,-1)= ITm(:,:)
   EcalT(:,:,+1)= ITp(:,:)
!   IF (Qz>Qx) THEN
!   IF (Qz<0.) THEN
!   ELSE
    VauxqzA= 0.   ! Initializes the integrand of the Qzp integrals
    VauxqzB= 0.
    DO kqzp= 1,nqz
     Qzp= VQz(kqzp)
     Qzp2= Qzp**2
     DO sigmap= -1,1,2
      Asp= sigmap/sigma
      Vauxqxp(:)= EcalL(:,kqzp,sigmap)
      DO sigmapp= -1,1,2
       Aspp= sigmapp/sigma
       DO S1= -1,1,2
        DO S2= -1,1,2
         Qzdif= Qz-S2*Qzp
         Qxp= QxpSLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
         ires= IQxpSLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
         iresdif= IQxpdifSLTd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
         Qxp2= Qxp**2
         Qxdif= Qx-S1*Qxp
         S3= INT(SIGN(1.D0,Qzdif))
         Zlqp= ZL(Qxp,Qzp)
         Ztqdif= ZT(Qxdif,ABS(Qzdif))
         D1= -S2*sigmap*(1.5*Qxp/SQRT(1.+1.5*(Qxp2+Qzp2)) &
          +S1*S3*sigmapp*Qxdif/Ve2C2/SQRT(1.+(Qxdif**2+Qzdif**2)/Ve2C2))
         AuxCoef= (Q2*(Qxp2+Qzp2)-(S1*Qx*Qxp+S2*Qz*Qzp)**2) &
          /(Qxp2+Qzp2+EpsMin)/(Qxdif**2+Qzdif**2+EpsMin)/(ABS(D1)+EpsMin)
         iresdif2= ABS(kqz-S2*kqzp)+S2
         !Vauxqxdif(:)= EcalT(:,ABS(iresdif2+S2*SIGN(1,iresdif2)),sigmapp)
         CALL Aux_Coef_Decay_z(nqxa,nqza,sigmapp,ires,iresdif,iresdif2,&
          Qxp,Qzp,Qxdif,Qzdif,"L","T",VQxa,VQza,Vauxqxp,EcalT,&
          Iqp,Iqdif)
         VauxqzA(kqzp)= VauxqzA(kqzp)+AuxCoef*sigma*Zlq*Iqp*Iqdif
         VauxqzB(kqzp)= VauxqzB(kqzp)-AuxCoef*(S2*sigmap*Zlqp*Iqdif &
           +S3*sigmapp*2.*Ztqdif*Iqp)
        END DO
       END DO
      END DO
     END DO
    END DO
    CALL Simpson(VQz,VauxqzA,nqz,Aux)
    CoefSLTdA= Aux
    CALL Simpson(VQz,VauxqzB,nqz,Aux)
    CoefSLTdB= Aux
!   END IF
   CoefSLTdA= (AA/2.)/VRNeNs(m)/SQRT(VRNeNs(m)*VRTeTs(m)) &
       *Q*(sigma*Zlq) * CoefSLTdA
   CoefSLTdB= (AA/2.)/VRNeNs(m)/SQRT(VRNeNs(m)*VRTeTs(m)) &
       *Q*(sigma*Zlq) * CoefSLTdB
  ELSE
   CoefSLTdA= 0.
   CoefSLTdB= 0.
  END IF

  IF (NewEffects2== "Yes") THEN
   IF(Gcoll== "Yes")THEN
    GcollSB= 2.*GcollSp(iqx,kqz)
   ELSE
    GcollSB=0.
   END IF
   IF(Bremss== "Yes")THEN
    BremSA= BremS(iqx,kqz)
   ELSE
    BremSA= 0.
   END IF
  ELSE
   GcollSB=0.
   BremSA= 0.
  END IF

  CoefA(iqx,kqz)= CoefEA+CoefSLLdA+CoefSLTdA+BremSA
  CoefB(iqx,kqz)= CoefEB+CoefSLLdB+CoefSLTdB+GcollSB
 END DO
END DO

! CALL Output_Coef("Swave",CoefA,CoefB)

RETURN
END SUBROUTINE Coef_Swave

SUBROUTINE Coef_Twave(sigma,CoefA,CoefB)
USE Common_Params
USE Common_Arrays
USE Math_Constants
USE Phys_Constants
!USE hyp_2f1_module
IMPLICIT NONE
REAL*8 :: Qx,Qz
REAL*8 :: Qx2,Qz2,Q2,Q
REAL*8 :: Zlq,Ztq,Zlqp,Zlqdif,Ztqp
REAL*8 :: Aux,AuxSE,AuxCoef,Aux0,Aux1
REAL*8 :: Auxe,Auxe2
REAL*8 :: Aux2,Aux3,Aux4
REAL*8 :: AuxScatElSpo
!REAL*8 :: AuxScatElInd
REAL*8 :: CoefTLLdA,CoefTLSdA,CoefTTLdA,CoefTLTsA
REAL*8 :: CoefTLLdB,CoefTLSdB,CoefTTLdB,CoefTLTsB
REAL*8 :: Qxp,Qzp,Qxp2,Qzp2
REAL*8 :: Qxpaux,Qzpaux
REAL*8 :: Qxdif,Qzdif
REAL*8 :: HBrL,HBrS,HGcL,HGcS
REAL*8 :: Iqp,Iqdif,Iqpp,Iqpm
REAL*8 :: D1
REAL*8 :: Iinit
REAL*8 :: Phip,Phi,Sphip,Cphip
REAL*8 :: Aalpha,AuxNum
REAL*8 :: Qstar,Qs
REAL*8, DIMENSION(nqx2) :: VintqxA,VintqxB
REAL*8, DIMENSION(nqz2) :: VintqzA,VintqzB
!REAL*8, DIMENSION(nqz) :: VauxqzA,VauxqzB
!REAL*8, DIMENSION(nqx) :: VauxqxA,VauxqxB
REAL*8, DIMENSION(nqx,nqz), INTENT(out) :: CoefA,CoefB
REAL*8, DIMENSION(nqx2,nqz2,-1:1) :: EcalL2
REAL*8, DIMENSION(nqx2,nqz2,-1:1) :: EcalS2
REAL*8, DIMENSION(nqx2,nqz2,-1:1) :: EcalT2
REAL*8, DIMENSION(nqz2) :: VauxqzA,VauxqzB
REAL*8, DIMENSION(nqx2) :: Vauxqxp!,Vauxqxdif
!REAL*8, DIMENSION(nqz) :: Vauxqzp,Vauxqzdif
REAL*8, DIMENSION(nqx2,nqz2) :: CoefAaux,CoefBaux
REAL*8, DIMENSION(nph) :: VintPhipA,VintPhipB
REAL*8, DIMENSION(nqx2) :: VQxa
REAL*8, DIMENSION(nqz2) :: VQza
INTEGER, INTENT(in) :: sigma
INTEGER :: i,k
INTEGER :: m,nqxa,nqza
INTEGER :: sigmap,sigmapp
INTEGER :: iqx,kqz,kqzp
INTEGER :: ires,iresdif,iresdif2,kres
INTEGER :: iresq,kresq
INTEGER :: S1,S2,S3,SS
INTEGER :: Asp,Aspp
!REAL*8 :: AuxEcalL
CHARACTER(LEN=1) :: WaveT

m= 1   ! Auxiliary to the space profiles (just one point, for the moment)
nqxa= nqx2
nqza= nqz2
VQxa= VQx2
VQza= VQz2
IF(SpontEmis=="Yes") THEN
 AuxSE= 1.E0
ELSE
 AuxSE= 0.E0
END IF
IF(ScatElSpo=="Yes") THEN
 AuxScatElSpo= 1.E0
ELSE
 AuxScatElSpo= 0.E0
END IF
!IF(ScatElInd=="Yes") THEN
! AuxScatElInd= 1.E0
!ELSE
! AuxScatElInd= 0.E0
!END IF

DO iqx= 1,nqx2
 Qx= VQx2(iqx)
 iresq= IQresx2(iqx)
 DO kqz= 1,nqz2
  Qz= VQz2(kqz)
  kresq= IQresz2(kqz)
  CALL Aitp2d2(nqx,nqz,VQx,VQz,ILm,Qx,Qz,Aux,iresq,kresq)
  EcalL2(iqx,kqz,-1)= Aux
  CALL Aitp2d2(nqx,nqz,VQx,VQz,ILp,Qx,Qz,Aux,iresq,kresq)
  EcalL2(iqx,kqz,+1)= Aux
  CALL Aitp2d2(nqx,nqz,VQx,VQz,ISm,Qx,Qz,Aux,iresq,kresq)
  EcalS2(iqx,kqz,-1)= Aux
  CALL Aitp2d2(nqx,nqz,VQx,VQz,ISp,Qx,Qz,Aux,iresq,kresq)
  EcalS2(iqx,kqz,+1)= Aux
  CALL Aitp2d2(nqx,nqz,VQx,VQz,ITm,Qx,Qz,Aux,iresq,kresq)
  EcalT2(iqx,kqz,-1)= Aux
  CALL Aitp2d2(nqx,nqz,VQx,VQz,ITp,Qx,Qz,Aux,iresq,kresq)
  EcalT2(iqx,kqz,+1)= Aux
 END DO
END DO

! Evaluation of the decay and scattering terms using the (nqx2,nqz2) grid, for 
! interpolation over the (nqx,nqz) grid:
CoefAaux= 0.
CoefBaux= 0.
DO iqx= 1,nqx2
 Qx= VQx2(iqx)
 DO kqz= 1,nqz2
  Qz= VQz2(kqz)
  Qx2= Qx**2
  Qz2= Qz**2
  Q2= Qx2+Qz2
  Q= SQRT(Q2)
  Zlq= ZL(Qx,Qz)
  Ztq= ZT(Qx,Qz)
  Phi= ACOS(Qz/Q)

  IF(TdecayLL=="Yes") THEN
   ! Contribution due to spontaneous and induced decay involving L waves:
!   IF (Qz>Qx) THEN
!   IF (Qz<0.) THEN
!   ELSE
    VauxqzA= 0.   ! Initializes the integrand of the Qzp integrals
    VauxqzB= 0.
    DO kqzp= 1,nqz2
     Qzp= VQz2(kqzp)
     Qzp2= Qzp**2
     DO sigmap= -1,1,2
      Asp= sigmap/sigma
      Vauxqxp(:)= EcalL2(:,kqzp,sigmap)
      DO sigmapp= -1,1,2
       Aspp= sigmapp/sigma
       DO S1= -1,1,2
        DO S2= -1,1,2
         Qxp= QxpTLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
         ires= IQxpTLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
         iresdif= IQxpdifTLLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
         Qxp2= Qxp**2
         Qxdif= Qx-S1*Qxp
         Qzdif= Qz-S2*Qzp
         S3= INT(SIGN(1.D0,Qzdif))
         Zlqp= ZL(Qxp,Qzp)
         Zlqdif= ZL(Qxdif,ABS(Qzdif))
         D1= -S2*sigmap*1.5*Qxp/SQRT(1.+1.5*(Qxp2+Qzp2)) &
          +S1*S3*sigmapp*1.5*Qxdif/SQRT(1.+1.5*(Qxdif**2+Qzdif**2))
         AuxCoef= (Q2*(Qxp2+Qzp2)-(S1*Qx*Qxp+S2*Qz*Qzp)**2) &
          /(Qxp2+Qzp2+EpsMin)/(Qxdif**2+Qzdif**2+EpsMin) &
          *((Qxp2+Qzp2)/S2/sigmap/Zlqp &
          -(Qxdif**2+Qzdif**2)/S3/sigmapp/Zlqdif)**2/(ABS(D1)+EpsMin)
         iresdif2= ABS(kqz-S2*kqzp)+S2
         !Vauxqxdif(:)= EcalL2(:,ABS(iresdif2+S2*SIGN(1,iresdif2)),sigmapp)
         CALL Aux_Coef_Decay_z(nqxa,nqza,sigmapp,ires,iresdif,iresdif2,&
          Qxp,Qzp,Qxdif,Qzdif,"L","L",VQxa,VQza,Vauxqxp,EcalL2,&
          Iqp,Iqdif)
         VauxqzA(kqzp)= VauxqzA(kqzp)+AuxCoef*sigma*2.*Ztq*Iqp*Iqdif
         VauxqzB(kqzp)= VauxqzB(kqzp)-AuxCoef*(S2*sigmap*Zlqp*Iqdif &
           +S3*sigmapp*Zlqdif*Iqp)
        END DO
       END DO
      END DO
     END DO
    END DO
    CALL Simpson(VQz2,VauxqzA,nqz2,Aux)
    CoefTLLdA= Aux
    CALL Simpson(VQz2,VauxqzB,nqz2,Aux)
    CoefTLLdB= Aux
!   END IF
   CoefTLLdA= (1./32.) * (sigma*Ztq/Q2) * CoefTLLdA
   CoefTLLdB= (1./32.) * (sigma*Ztq/Q2) * CoefTLLdB
  ELSE
   CoefTLLdA= 0.
   CoefTLLdB= 0.
  END IF

  IF(TdecayLS=="Yes") THEN
   ! Contribution due to spontaneous and induced decay involving L and S waves:
!   IF (Qz>Qx) THEN
!   IF (Qz<0.) THEN
!   ELSE
    VauxqzA= 0.   ! Initializes the integrand of the Qzp integrals
    VauxqzB= 0.
    DO kqzp= 1,nqz2
     Qzp= VQz2(kqzp)
     Qzp2= Qzp**2
     DO sigmap= -1,1,2
      Asp= sigmap/sigma
      Vauxqxp(:)= EcalL2(:,kqzp,sigmap)
      DO sigmapp= -1,1,2
       Aspp= sigmapp/sigma
       DO S1= -1,1,2
        DO S2= -1,1,2
         Qxp= QxpTLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
         ires= IQxpTLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
         iresdif= IQxpdifTLSd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
         Qxp2= Qxp**2
         Qxdif= Qx-S1*Qxp
         Qzdif= Qz-S2*Qzp
         S3= INT(SIGN(1.D0,Qzdif))
         Zlqp= ZL(Qxp,Qzp)
         Zlqdif= ZL(Qxdif,ABS(Qzdif))
         D1= -S2*sigmap*(1.5*Qxp/SQRT(1.+1.5*(Qxp2+Qzp2))) &
          + S1*S3*sigmapp*AA*Qxdif/SQRT(Qxdif**2+Qzdif**2+EpsMin) &
          /(1.+(Qxdif**2+Qzdif**2)/2.)**(1.5)
         AuxCoef= SQRT((Qxdif)**2+(Qzdif)**2) &
          *(Q2*(Qxp2+Qzp2)-(S1*Qx*Qxp+S2*Qz*Qzp)**2) &
          /(Qxp2+Qzp2+EpsMin)/(ABS(D1)+EpsMin)
         iresdif2= ABS(kqz-S2*kqzp)+S2
         !Vauxqxdif(:)= EcalS2(:,ABS(iresdif2+S2*SIGN(1,iresdif2)),sigmapp)
         CALL Aux_Coef_Decay_z(nqxa,nqza,sigmapp,ires,iresdif,iresdif2,&
          Qxp,Qzp,Qxdif,Qzdif,"L","S",VQxa,VQza,Vauxqxp,EcalS2,&
          Iqp,Iqdif)
         VauxqzA(kqzp)= VauxqzA(kqzp)+AuxCoef*sigma*2.*Ztq*Iqp*Iqdif
         VauxqzB(kqzp)= VauxqzB(kqzp)-AuxCoef*(S2*sigmap*Zlqp*Iqdif &
            +S3*sigmapp*Zlqdif*Iqp)
        END DO
       END DO
      END DO
     END DO
    END DO
    CALL Simpson(VQz2,VauxqzA,nqz2,Aux)
    CoefTLSdA= Aux
    CALL Simpson(VQz2,VauxqzB,nqz2,Aux)
    CoefTLSdB= Aux
!   END IF
   CoefTLSdA= (AA/2.) * (sigma*Ztq/Q2) * CoefTLSdA
   CoefTLSdB= (AA/2.) * (sigma*Ztq/Q2) * CoefTLSdB
  ELSE
   CoefTLSdA= 0.
   CoefTLSdB= 0.
  END IF

  IF(TdecayTL=="Yes") THEN
   ! Contribution due to spontaneous and induced decay involving L and S waves:
!   IF (Qz>Qx) THEN
!   IF (Qz<0.) THEN
!   ELSE
    VauxqzA= 0.   ! Initializes the integrand of the Qzp integrals
    VauxqzB= 0.
    DO kqzp= 1,nqz2
     Qzp= VQz2(kqzp)
     Qzp2= Qzp**2
     DO sigmap= -1,1,2
      Asp= sigmap/sigma
      Vauxqxp(:)= EcalT2(:,kqzp,sigmap)
      DO sigmapp= -1,1,2
       Aspp= sigmapp/sigma
       DO S1= -1,1,2
        DO S2= -1,1,2
         Qxp= QxpTTLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
         ires= IQxpTTLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
         iresdif= IQxpdifTTLd(iqx,kqz,kqzp,Asp,Aspp,S1,S2)
         Qxp2= Qxp**2
         Qxdif= Qx-S1*Qxp
         Qzdif= Qz-S2*Qzp
         S3= INT(SIGN(1.D0,Qzdif))
         Ztqp= ZT(Qxp,Qzp)
         Zlqdif= ZL(Qxdif,ABS(Qzdif))
         D1= -S2*sigmap*(Qxp/Ve2C2/SQRT(1.+(Qxp2+Qzp2)/Ve2C2)) &
          + S1*S3*sigmapp*1.5*Qxdif/SQRT(1.+1.5*(Qxdif**2+Qzdif**2))
         AuxCoef= ((Qxdif)**2+(Qzdif)**2)/Ztqp**2 &
          *(Q2+(S1*Qx*Qxp+S2*Qz*Qzp)**2/(Qxp2+Qzp2+EpsMin))/(ABS(D1)+EpsMin)
         iresdif2= ABS(kqz-S2*kqzp)+S2
         !Vauxqxdif(:)= EcalL2(:,ABS(iresdif2+S2*SIGN(1,iresdif2)),sigmapp)
         CALL Aux_Coef_Decay_z(nqxa,nqza,sigmapp,ires,iresdif,iresdif2,&
          Qxp,Qzp,Qxdif,Qzdif,"T","L",VQxa,VQza,Vauxqxp,EcalL2,&
          Iqp,Iqdif)
         VauxqzA(kqzp)= VauxqzA(kqzp)+AuxCoef*sigma*2.*Ztq*Iqp*Iqdif
         VauxqzB(kqzp)= VauxqzB(kqzp)-AuxCoef*(S2*sigmap*2.*Ztqp*Iqdif &
            +S3*sigmapp*Zlqdif*Iqp)
        END DO
       END DO
      END DO
     END DO
    END DO
    CALL Simpson(VQz2,VauxqzA,nqz2,Aux)
    CoefTTLdA= Aux
    CALL Simpson(VQz2,VauxqzB,nqz2,Aux)
    CoefTTLdB= Aux
!   END IF
   CoefTTLdA= (sigma*Ztq/Ztq**2/8./Q2) * CoefTTLdA
   CoefTTLdB= (sigma*Ztq/Ztq**2/8./Q2) * CoefTTLdB
  ELSE
   CoefTTLdA= 0.
   CoefTTLdB= 0.
  END IF

  IF(TscatLT=="Yes") THEN
   ! Contribution due to spontaneous and induced scattering involving L and T
   ! waves:
!   EcalL(:,:,-1)= ILm(:,:)
!   EcalL(:,:,+1)= ILp(:,:)
   ! Electron contribution, Kappa:
   VintqxA= 0.
   VintqxB= 0.
   DO i= 1,nqx2
    Qxp= VQx3(i)
    VintqzA= 0.
    VintqzB= 0.
    DO k= 1,nqz2
     Qzp= VQz3(k)
     Zlqp= ZL(Qxp,Qzp)
     IF( (Qxp>=VQx(nqx) .OR. Qxp<=VQx(1)) .OR. &
      (Qzp>=VQz(nqz) .OR. Qzp<=VQz(1)) ) THEN
      Iqpp= ILgrid3(i,k)
      Iqpm= ILgrid3(i,k)
     ELSE
      ires= IQresx3(i)
      kres= IQresz3(k)
      CALL Aitp2d2(nqx,nqz,VQx,VQz,ILp,Qxp,Qzp,Iqpp,ires,kres)
      CALL Aitp2d2(nqx,nqz,VQx,VQz,ILm,Qxp,Qzp,Iqpm,ires,kres)
     END IF
     DO sigmap= -1,1,2
      DO S1= -1,1,2
       Qxdif= Qx-S1*Qxp
       S2= sigma*sigmap  ! The other value of S2 gives negligible contribution.
       Qzdif= Qz-S2*Qzp
       IF((Qxdif**2+Qzdif**2)<=1.E-6) THEN
        Auxe= 0.
        Auxe2= 0.
       ELSE
        Aux0= (Q2*(Qxp**2+Qzp**2)-(S1*Qx*Qxp+S2*Qz*Qzp)**2) &
         /(Qxp**2+Qzp**2+EpsMin) &
         /SQRT(Qxdif**2+Qzdif**2)
        Aux= (sigma*Ztq-sigma*Zlqp)**2/(Qxdif**2+Qzdif**2)
        Aux1= Kappae*Uek2
        Aux2= 1.+Aux/Aux1
        Aux3= Aux0 * (RatioNke/SQRT(Aux1)/Aux2**(Kappae+Alphae-0.5))
        Aux4= Aux3 * 2./Aux1/Aux2
        Auxe= Aux3 * RatioGammaem051 
        Auxe2= Aux4 * RatioGammaep051 * (sigma*Ztq-sigma*Zlqp) 
       END IF
       VintqzA(k)= VintqzA(k) + (Auxe*AuxScatElSpo)*(sigma*Ztq) &
        * AuxSE*Geff*(sigma*Ztq)*((1-sigmap)*Iqpm+(1+sigmap)*Iqpp)/2.
       VintqzB(k)= VintqzB(k) &
        + (Auxe*AuxScatElSpo)*(AuxSE*Geff*(sigma*Ztq)*(-sigma*Zlqp)) &
        + Auxe2*(sigma*Ztq-sigma*Zlqp)*((1-sigmap)*Iqpm+(1+sigmap)*Iqpp)/2.
      END DO ! S1
     END DO ! sigmap
    END DO ! Qzp
    CALL Simpson(VQz3,VintqzA,nqz2,Aux)
    VintqxA(i)= Aux
    CALL Simpson(VQz3,VintqzB,nqz2,Aux)
    VintqxB(i)= Aux
   END DO   ! Qxp
   CALL Simpson(VQx3,VintqxA,nqx2,Aux)
   CoefTLTsA= 1./SQRT(Pi)/VRNeNs(m)/Q2 * Aux
   CALL Simpson(VQx3,VintqxB,nqx2,Aux)
   CoefTLTsB= 1./SQRT(Pi)/VRNeNs(m)/Q2 * Aux/2.
   ! Electron and ion contributions, Maxw:
   Qstar= SQRT(2./3./Ve2C2/VRTeTs(m))*Q
   VintPhipA= 0.
   VintPhipB= 0.
   DO i= 1,nph
    Phip= VPhip(i)
    Qxp= Qstar*SIN(Phip)
    Qzp= Qstar*COS(Phip)
    IF( (Qxp>=VQx(nqx)) .OR. (Qzp>=VQz(nqz)) ) THEN
     HBrL=0. ! The new effects are not relevant for large Q and are neglected
     HBrS=0.
     HGcL=0.
     HGcS=0.
     WaveT= "L"
     CALL Iwave_Init(Qxp,Qzp,Iinit,WaveT,HBrL,HBrS,HGcL,HGcS)
     Iqpp= Iinit
     Iqpm= Iinit
    ELSE
     IF (Qxp<=VQx(1)) THEN
      ires= 1
      Qxpaux= VQx(1)+EpsMin
      IF (Qzp<=VQz(1)) THEN
       kres= 1
       Iqpp= ILp(1,1)
       Iqpm= ILm(1,1)
      ELSE
       kres= IQzTLT1(iqx,kqz,i)
       Qzpaux= Qzp
       CALL Aitp2d2(nqx,nqz,VQx,VQz,ILp,Qxpaux,Qzpaux,Iqpp,ires,kres)
       CALL Aitp2d2(nqx,nqz,VQx,VQz,ILm,Qxpaux,Qzpaux,Iqpm,ires,kres)
      END IF
     ELSE
      ires= IQxTLT1(iqx,kqz,i)
      Qxpaux= Qxp
      IF (Qzp<=VQz(1)) THEN
       kres= 1
       Qzpaux= VQz(1)+EpsMin
       CALL Aitp2d2(nqx,nqz,VQx,VQz,ILp,Qxpaux,Qzpaux,Iqpp,ires,kres)
       CALL Aitp2d2(nqx,nqz,VQx,VQz,ILm,Qxpaux,Qzpaux,Iqpm,ires,kres)
      ELSE
       kres= IQzTLT1(iqx,kqz,i)
       Qzpaux= Qzp
       CALL Aitp2d2(nqx,nqz,VQx,VQz,ILp,Qxpaux,Qzpaux,Iqpp,ires,kres)
       CALL Aitp2d2(nqx,nqz,VQx,VQz,ILm,Qxpaux,Qzpaux,Iqpm,ires,kres)
      END IF
     END IF
    END IF
    DO S1= -1,1,2
     sigmap= 1.
     Aux= (SIN(Phi-S1*sigma*sigmap*Phip))**2
     VintPhipA(i)= VintPhipA(i)+Aux*Iqpp
     VintPhipB(i)= VintPhipB(i)+Aux
     sigmap= -1.
     Aux= (SIN(Phi-S1*sigma*sigmap*Phip))**2
     VintPhipA(i)= VintPhipA(i)+Aux*Iqpm
     VintPhipB(i)= VintPhipB(i)+Aux
    END DO
   END DO
   CALL Simpson(VPhip,VintPhipA,nph,Aux)
   CoefTLTsA= CoefTLTsA &
    + (2./3.)/VRNeNs(m) * (Ztq)**3 * AuxSE*Geff*Aux * (1.-RatioNki) &
    + (2./3.)/VRNeNs(m) * (Ztq)**3 * AuxScatElSpo*AuxSE*Geff*Aux * (1.-RatioNke)
   CALL Simpson(VPhip,VintPhipB,nph,Aux)
   CoefTLTsB= CoefTLTsB &
    - (2./3.)/VRNeNs(m) * (Ztq)**3 * AuxSE*Geff*Aux/2. * (1.-RatioNki) &
    - (2./3.)/VRNeNs(m) * (Ztq)**3 * AuxScatElSpo*AuxSE*Geff*Aux/2. &
    * (1.-RatioNke)
   ! Induced scattering due to the ions:
   VintPhipB= 0.
   DO i= 1,nph
    Phip= VPhip(i)
    Sphip= SIN(Phip)
    Cphip= COS(Phip)
    DO sigmap= -1,1,2
     DO S1= -1,1,2
      Aux1= (Q2+Qstar**2-2.*S1*Q*Qstar*COS(Phi-S1*sigma*sigmap*Phip))
      IF (Aux1<1.E-6) THEN
       AuxNum= SQRT(1.E-6)
       Aalpha= RMiMe/VRTiTs(m)*(9./4.)*(Qstar/Ztq)**2/1.E-6
      ELSE
       AuxNum= SQRT(Aux1)
       Aalpha= RMiMe/VRTiTs(m)*(9./4.)*(Qstar/Ztq)**2/Aux1 
      END IF
      Aux= (SIN(Phi-S1*sigma*sigmap*Phip))**2*AuxNum
      DO SS= -1,1,2
       Qs= Qstar+SS/SQRT(2.*Aalpha)
       Qxp= Qs*Sphip
       Qzp= Qs*Cphip
       IF( (Qxp>=VQx(nqx)) .OR. (Qzp>=VQz(nqz)) ) THEN
        HBrL=0. ! The new effects are not relevant for large Q and are neglected
        HBrS=0.
        HGcL=0.
        HGcS=0.
        WaveT= "L"
        CALL Iwave_Init(Qxp,Qzp,Iinit,WaveT,HBrL,HBrS,HGcL,HGcS)
        Iqp= Iinit
       ELSE
        S2= sigma*sigmap
        IF (Qxp<=VQx(1)) THEN
         ires= 1
         Qxpaux= VQx(1)+EpsMin
         IF (Qzp<=VQz(1)) THEN
          kres= 1
          IF (sigmap==1) THEN
           Iqp= ILp(1,1) 
          ELSE
           Iqp= ILm(1,1) 
          END IF
         ELSE
          kres= IQzTLT2i(iqx,kqz,i,S2,S1,SS)
          Qzpaux= Qzp
          IF (sigmap==1) THEN
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ILp,Qxpaux,Qzpaux,Iqpp,ires,kres)
          ELSE
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ILm,Qxpaux,Qzpaux,Iqpm,ires,kres)
          END IF
         END IF
        ELSE
         ires= IQxTLT2i(iqx,kqz,i,S2,S1,SS)
         Qxpaux= Qxp
         IF (Qzp<=VQz(1)) THEN
          kres= 1
          Qzpaux= VQz(1)+EpsMin
          IF (sigmap==1) THEN
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ILp,Qxpaux,Qzpaux,Iqpp,ires,kres)
          ELSE
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ILm,Qxpaux,Qzpaux,Iqpm,ires,kres)
          END IF
         ELSE
          kres= IQzTLT2i(iqx,kqz,i,S2,S1,SS)
          Qzpaux= Qzp
          IF (sigmap==1) THEN
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ILp,Qxpaux,Qzpaux,Iqpp,ires,kres)
          ELSE
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ILm,Qxpaux,Qzpaux,Iqpm,ires,kres)
          END IF
         END IF
        END IF
       END IF
       VintPhipB(i)= VintPhipB(i)+Aux*Iqp * SS*Qs/Qstar 
      END DO ! SS
     END DO ! S1
    END DO ! sigmap
   END DO
   CALL Simpson(VPhip,VintPhipB,nph,Aux)
   CoefTLTsB= CoefTLTsB + (4./3.)/VRNeNs(m)/SQRT(VRTiTs(m)*RMiMe) &
      * EXP(-0.5)/SQRT(2.) * (Ztq)**2 * Aux/2. * (1.-RatioNki)
   ! Induced scattering due to the electrons:
   VintPhipB= 0.
   DO i= 1,nph
    Phip= VPhip(i)
    Sphip= SIN(Phip)
    Cphip= COS(Phip)
    DO sigmap= -1,1,2
     DO S1= -1,1,2
      Aux1= (Q2+Qstar**2-2.*S1*Q*Qstar*COS(Phi-S1*sigma*sigmap*Phip))
      IF (Aux1<1.E-6) THEN
       AuxNum= SQRT(1.E-6)
       Aalpha= 1.D0/VRTeTs(m)*(9./4.)*(Qstar/Ztq)**2/1.E-6
      ELSE
       AuxNum= SQRT(Aux1)
       Aalpha= 1.D0/VRTeTs(m)*(9./4.)*(Qstar/Ztq)**2/Aux1 
      END IF
      Aux= (SIN(Phi-S1*sigma*sigmap*Phip))**2*AuxNum
      DO SS= -1,1,2
       Qs= Qstar+SS/SQRT(2.*Aalpha)
       Qxp= Qs*Sphip
       Qzp= Qs*Cphip
       IF( (Qxp>=VQx(nqx)) .OR. (Qzp>=VQz(nqz)) ) THEN
        HBrL=0. ! The new effects are not relevant for large Q and are neglected
        HBrS=0.
        HGcL=0.
        HGcS=0.
        WaveT= "L"
        CALL Iwave_Init(Qxp,Qzp,Iinit,WaveT,HBrL,HBrS,HGcL,HGcS)
        Iqp= Iinit
       ELSE
        S2= sigma*sigmap
        IF (Qxp<=VQx(1)) THEN
         ires= 1
         Qxpaux= VQx(1)+EpsMin
         IF (Qzp<=VQz(1)) THEN
          kres= 1
          IF (sigmap==1) THEN
           Iqp= ILp(1,1) 
          ELSE
           Iqp= ILm(1,1) 
          END IF
         ELSE
          kres= IQzTLT2e(iqx,kqz,i,S2,S1,SS)
          Qzpaux= Qzp
          IF (sigmap==1) THEN
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ILp,Qxpaux,Qzpaux,Iqpp,ires,kres)
          ELSE
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ILm,Qxpaux,Qzpaux,Iqpm,ires,kres)
          END IF
         END IF
        ELSE
         ires= IQxTLT2e(iqx,kqz,i,S2,S1,SS)
         Qxpaux= Qxp
         IF (Qzp<=VQz(1)) THEN
          kres= 1
          Qzpaux= VQz(1)+EpsMin
          IF (sigmap==1) THEN
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ILp,Qxpaux,Qzpaux,Iqpp,ires,kres)
          ELSE
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ILm,Qxpaux,Qzpaux,Iqpm,ires,kres)
          END IF
         ELSE
          kres= IQzTLT2e(iqx,kqz,i,S2,S1,SS)
          Qzpaux= Qzp
          IF (sigmap==1) THEN
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ILp,Qxpaux,Qzpaux,Iqpp,ires,kres)
          ELSE
           CALL Aitp2d2(nqx,nqz,VQx,VQz,ILm,Qxpaux,Qzpaux,Iqpm,ires,kres)
          END IF
         END IF
        END IF
       END IF
       VintPhipB(i)= VintPhipB(i)+Aux*Iqp * SS*Qs/Qstar * (Ztq-Zlqp)
      END DO ! SS
     END DO ! S1
    END DO ! sigmap
   END DO
   CALL Simpson(VPhip,VintPhipB,nph,Aux)
   CoefTLTsB= CoefTLTsB - (4./3.)/VRNeNs(m)/SQRT(VRTeTs(m)) &
      * EXP(-0.5)/SQRT(2.) * (Ztq) * Aux/2. * (1.-RatioNke)
  ELSE
   CoefTLTsA= 0.
   CoefTLTsB= 0.
  END IF

  CoefAaux(iqx,kqz)= CoefTLLdA+CoefTLSdA+CoefTTLdA+CoefTLTsA
  CoefBaux(iqx,kqz)= CoefTLLdB+CoefTLSdB+CoefTTLdB+CoefTLTsB
 END DO
END DO

CoefA= 0.
CoefB= 0.
DO iqx= 1,nqx
 Qx= VQx(iqx)
 iresq= IQresx(iqx)
 DO kqz= 1,nqz
  Qz= VQz(kqz)
  kresq= IQresz(kqz)

  ! Interpolation of the contribution due to decay terms (TLLd, TLSd, TTLd),
  ! and the contribution of spontaneous and induced scattering involving L and
  ! T waves (TLTs):
  ! Routine "Aitp2d2b" corrects the interpolation at the extreme points.
  CALL Aitp2d2b(nqx2,nqz2,VQx2,VQz2,CoefAaux,Qx,Qz,Aux,iresq,kresq)
  CoefA(iqx,kqz)= Aux
  CALL Aitp2d2b(nqx2,nqz2,VQx2,VQz2,CoefBaux,Qx,Qz,Aux,iresq,kresq)
  CoefB(iqx,kqz)= Aux
 END DO
END DO

! CALL Output_Coef("Twave",CoefA,CoefB)

RETURN
END SUBROUTINE Coef_Twave
  
SUBROUTINE Aux_Coef_Decay_z(nqxa,nqza,sigmapp,ires,iresdif,iresdif2,&
  Qxp,Qzp,Qxdif,Qzdif,Wavep,Wavedif,VQxa,VQza,Vauxqxp,EcalAux,&
  Iqp,Iqdif)
USE Common_Params
USE Common_Arrays
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
INTEGER, INTENT(in) :: nqxa,nqza,sigmapp,ires,iresdif,iresdif2
REAL*8, INTENT(in) :: Qxp,Qzp,Qxdif,Qzdif
REAL*8, INTENT(out) :: Iqp,Iqdif
REAL*8 :: Iinit
REAL*8 :: HBrL,HBrS,HGcL,HGcS
REAL*8 :: Aux
REAL*8, DIMENSION(nqxa), INTENT(in) :: VQxa
REAL*8, DIMENSION(nqza), INTENT(in) :: VQza
REAL*8, DIMENSION(nqxa), INTENT(in) :: Vauxqxp
REAL*8, DIMENSION(nqxa) :: Vauxqxdif
REAL*8, DIMENSION(nqxa,nqza,-1:1), INTENT(in) :: EcalAux
CHARACTER(LEN=1), INTENT(in) :: Wavep,Wavedif

Aux= VQza(1)   ! Only to avoid a warning message at compilation.

SELECT CASE(ires)
 CASE(-1)   ! Qxp<0
  Iqp= 0.
 CASE(0)
  Iqp= Vauxqxp(1)
 CASE(1:)   ! ires >= 1
  IF (ires<nqxa) THEN
   CALL Aitp1d2(nqxa,VQxa,Vauxqxp,Qxp,Iqp,ires)
  ELSE
   HBrL=0. ! The new effects are not relevant for large Q and are neglected
   HBrS=0.
   HGcL=0.
   HGcS=0.
   CALL Iwave_Init(Qxp,Qzp,Iinit,Wavep,HBrL,HBrS,HGcL,HGcS)
   Iqp= Iinit
  END IF
 CASE DEFAULT
 OPEN(98,FILE='Warning_Select_Case_Aux.wt')
 WRITE(98,*) ' ires= ',ires
 WRITE(98,*) ' ires must be integer greater or equal -1 !!'
 CLOSE(98)
END SELECT

IF (ires==(-1)) THEN
 Iqdif= 0.   ! Qxp<0
ELSE
 SELECT CASE(iresdif)
  CASE(-1)   ! Qxdif<0
   Iqdif= 0.
  CASE(0)
   SELECT CASE(iresdif2)
    CASE(-1)   ! Qzdif<0
     Iqdif= 0.
    CASE(0)
     Iqdif= EcalAux(1,1,sigmapp)
    CASE(1:)
     IF (iresdif2<nqza) THEN
      Iqdif= EcalAux(1,iresdif2,sigmapp)
     ELSE
      HBrL=0. ! The new effects are not relevant for large Q and are neglected
      HBrS=0.
      HGcL=0.
      HGcS=0.
      CALL Iwave_Init(ABS(VQx(1)),ABS(Qzdif),Iinit,Wavedif,HBrL,HBrS,HGcL,HGcS)
      Iqdif= Iinit
     END IF
    CASE DEFAULT
    OPEN(98,FILE='Warning_Select_Case_Aux.wt')
    WRITE(98,*) ' iresdif2= ',iresdif2
    WRITE(98,*) ' iresdif2 must be integer greater or equal -1 !!'
    CLOSE(98)
   END SELECT
  CASE(1:)   ! iresdif >= 1
   IF (iresdif<nqxa) THEN
    SELECT CASE(iresdif2)
     CASE(-1)   ! Qzdif<0
      Iqdif= 0.
     CASE(0)
      Vauxqxdif(:)= EcalAux(:,1,sigmapp)
      CALL Aitp1d2(nqxa,VQxa,Vauxqxdif,ABS(Qxdif),Iqdif,iresdif)
     CASE(1:)
      IF (iresdif2<nqza) THEN
       Vauxqxdif(:)= EcalAux(:,iresdif2,sigmapp)
       CALL Aitp1d2(nqxa,VQxa,Vauxqxdif,ABS(Qxdif),Iqdif,iresdif)
      ELSE
       HBrL=0. ! The new effects are not relevant for large Q and are neglected
       HBrS=0.
       HGcL=0.
       HGcS=0.
       CALL Iwave_Init(ABS(Qxdif),ABS(Qzdif),Iinit,Wavedif,HBrL,HBrS,HGcL,HGcS)
       Iqdif= Iinit
      END IF
     CASE DEFAULT
     OPEN(98,FILE='Warning_Select_Case_Aux.wt')
     WRITE(98,*) ' iresdif2= ',iresdif2
     WRITE(98,*) ' iresdif2 must be integer greater or equal -1 !!'
     CLOSE(98)
    END SELECT
   ELSE
    HBrL=0. ! The new effects are not relevant for large Q and are neglected
    HBrS=0.
    HGcL=0.
    HGcS=0.
    CALL Iwave_Init(ABS(Qxdif),ABS(Qzdif),Iinit,Wavedif,HBrL,HBrS,HGcL,HGcS)
    Iqdif= Iinit
   END IF
  CASE DEFAULT
  OPEN(98,FILE='Warning_Select_Case_Aux.wt')
  WRITE(98,*) ' iresdif= ',iresdif
  WRITE(98,*) ' iresdif must be integer greater or equal -1 !!'
  CLOSE(98)
 END SELECT
END IF

RETURN
END SUBROUTINE Aux_Coef_Decay_z

REAL*8 FUNCTION ZL(Qx,Qz)
USE Common_Params
USE Common_Arrays
IMPLICIT NONE
REAL*8, INTENT(in) :: Qx,Qz
REAL*8 :: Q2
INTEGER :: m

m= 1
Q2= Qx**2+Qz**2
ZL= SQRT(VRNeNs(m))*SQRT(1.+1.5*Q2*VRTeTs(m)/VRNeNs(m))
IF (Qz < 0.) THEN
 ZL= -ZL
ELSE
END IF
RETURN
END FUNCTION ZL

REAL*8 FUNCTION ZS(Qx,Qz)
USE Common_Params
USE Common_Arrays
IMPLICIT NONE
REAL*8, INTENT(in) :: Qx,Qz
REAL*8 :: Q,Q2
INTEGER :: m

m= 1   ! Auxiliary quantity, for the spatial profile.
Q2= Qx**2+Qz**2
Q= SQRT(Q2)
ZS= Q*AA*SQRT(VRTeTs(m))/SQRT(1.+Q2/2.*VRTeTs(m)/VRNeNs(m))
IF (Qz < 0.) THEN
 ZS= -ZS
ELSE
END IF
RETURN
END FUNCTION ZS

REAL*8 FUNCTION ZT(Qx,Qz)
USE Common_Params
USE Common_Arrays
IMPLICIT NONE
REAL*8, INTENT(in) :: Qx,Qz
REAL*8 :: Q2

Q2= Qx**2+Qz**2
ZT= SQRT(1.+Q2/Ve2C2)
IF (Qz < 0.) THEN
 ZT= -ZT
ELSE
END IF
RETURN
END FUNCTION ZT

      SUBROUTINE ZBRAC2(FUNC,X1,X2,SUCCES) 
      IMPLICIT NONE
      REAL*8, INTENT(INOUT) :: X1,X2
      REAL*8, PARAMETER :: FACTOR= 1.6
      REAL*8 F1,F2 
      INTEGER, PARAMETER :: NTRY= 10
      INTEGER J
      LOGICAL, INTENT(OUT) :: SUCCES 
!      REAL*8X1,X2,FUNC,FACTOR
!      INTEGER NTRY,J
!      LOGICAL SUCCES 
!      EXTERNAL FUNC
!      PARAMETER (FACTOR=1.6,NTRY=10)
!        Given a function FUNC and an initial guessed range X1 to X2, the
!        routine expands the range geometrically until a root is bracketed by
!        the returned values X1 and X2 (in which case succes returns as .true.)
!        or until the range becomes unacceptably large (in which case succes
!        returns as .false.).
!        This is a version of ZBRAC. It expands only the upper limit of the 
!        range.
      INTERFACE 
       REAL*8 FUNCTION FUNC(X)
       REAL*8, INTENT(in) :: X
       END FUNCTION FUNC
      END INTERFACE
      !IF(X1.EQ.X2)THEN
      IF(ABS(X1-X2).LT.1.E-16)THEN
       OPEN(98,FILE='Warning_Zbrac2.wt')
       WRITE(98,*) 'You have to guess an initial range in zbrac'
       CLOSE(98)
       STOP
      ELSE
      END IF 
      F1=FUNC(X1) 
      F2=FUNC(X2) 
      SUCCES=.TRUE. 
      DO J=1,NTRY 
       IF(F1*F2.LT.0.)RETURN 
        X2=X2+FACTOR*(X2-X1) 
        F2=FUNC(X2) 
      END DO 
      SUCCES=.FALSE. 
      RETURN 
      END SUBROUTINE ZBRAC2

      REAL*8 FUNCTION RTBIS(FUNC,X1,X2,XACC) 
      IMPLICIT NONE
      INTEGER JMAX,J
      REAL*8 X1,X2,XACC
      REAL*8 F,FMID,DX,XMID
!      REAL*8RTBIS,X1,X2,XACC,FUNC
!      EXTERNAL FUNC
      PARAMETER (JMAX=40)     ! Maximum allowed number of bisections.
!         Using bisection, find the root of a function FUNC known to lie
!         between x1 and x2. The root, returned as RTBIS, will be refined
!         until its accuracy is +-XACC.
      INTERFACE 
       REAL*8 FUNCTION FUNC(X)
       REAL*8, INTENT(in) :: X
       END FUNCTION FUNC
      END INTERFACE
      FMID=FUNC(X2) 
      F=FUNC(X1) 
      IF(F*FMID.GE.0.)THEN
       OPEN(98,FILE='Warning_Rtbis.wt')
       WRITE(98,*)'Root must be bracketed in rtbis.'
       CLOSE(98)
       STOP
      ELSE
      END IF 
      IF(F.LT.0.)THEN         ! Orient the search so that f>0 lies at x+dx
        RTBIS=X1 
        DX=X2-X1 
      ELSE 
        RTBIS=X2 
        DX=X1-X2 
      ENDIF 
      DO J=1,JMAX           ! Bisection loop
        DX=DX*.5 
        XMID=RTBIS+DX 
        FMID=FUNC(XMID) 
        IF(FMID.LT.0.)RTBIS=XMID 
        !IF(ABS(DX).LT.XACC .OR. FMID.EQ.0.) RETURN 
        IF(ABS(DX).LT.XACC .OR. ABS(FMID).LT.1.e-16) RETURN 
      END DO
      OPEN(98,FILE='Warning_Rtbis.wt')
      WRITE(98,*)'Too many bisections in rtbis' 
      CLOSE(98)
      STOP
      END FUNCTION RTBIS

REAL*8 FUNCTION Funcx_RcdLLS(Qxp)
USE Common_Params
USE Common_Arrays
USE Math_Constants
IMPLICIT NONE
REAL*8, INTENT(in) :: Qxp
REAL*8 :: Qx,Qz,Qzp,Qxdif,Qzdif,Beta
INTEGER :: Asp,Aspp,S1,S2,S3 
INTEGER :: m

Qx= Aux1_Rcd(1)
Qz= Aux1_Rcd(2)
Qzp= Aux1_Rcd(3)
Asp= Aux2_Rcd(1)
Aspp= Aux2_Rcd(2)
S1= Aux2_Rcd(3)
S2= Aux2_Rcd(4)
Qxdif= Qx-S1*Qxp
Qzdif= Qz-S2*Qzp
S3= INT(SIGN(1.D0,Qzdif))
m= 1   ! Auxiliary to the space profiles (just one point, for the moment)
Beta= VRTeTs(m)/VRNeNs(m)

Funcx_RcdLLS= ZL(Qx,ABS(Qz))-Asp*S2*ZL(Qxp,ABS(Qzp)) - Aspp*S3*ZS(Qxdif,ABS(Qzdif))
RETURN
END FUNCTION Funcx_RcdLLS

REAL*8 FUNCTION Funcx_RcdLLT(Qxp)
USE Common_Params
USE Common_Arrays
USE Math_Constants
IMPLICIT NONE
REAL*8, INTENT(in) :: Qxp
REAL*8 :: Qx,Qz,Qzp,Qxdif,Qzdif,Beta
INTEGER :: Asp,Aspp,S1,S2,S3 
INTEGER :: m

Qx= Aux1_Rcd(1)
Qz= Aux1_Rcd(2)
Qzp= Aux1_Rcd(3)
Asp= Aux2_Rcd(1)
Aspp= Aux2_Rcd(2)
S1= Aux2_Rcd(3)
S2= Aux2_Rcd(4)
Qxdif= Qx-S1*Qxp
Qzdif= Qz-S2*Qzp
S3= INT(SIGN(1.D0,Qzdif))
m= 1   ! Auxiliary to the space profiles (just one point, for the moment)
Beta= VRTeTs(m)/VRNeNs(m)

Funcx_RcdLLT= ZL(Qx,ABS(Qz))-Asp*S2*ZL(Qxp,ABS(Qzp)) &
    -Aspp*S3*ZT(Qxdif,ABS(Qzdif))
RETURN
END FUNCTION Funcx_RcdLLT

REAL*8 FUNCTION Funcx_RcdLST(Qxp)
USE Common_Params
USE Common_Arrays
USE Math_Constants
IMPLICIT NONE
REAL*8, INTENT(in) :: Qxp
REAL*8 :: Qx,Qz,Qzp,Qxdif,Qzdif,Beta
INTEGER :: Asp,Aspp,S1,S2,S3 
INTEGER :: m

Qx= Aux1_Rcd(1)
Qz= Aux1_Rcd(2)
Qzp= Aux1_Rcd(3)
Asp= Aux2_Rcd(1)
Aspp= Aux2_Rcd(2)
S1= Aux2_Rcd(3)
S2= Aux2_Rcd(4)
Qxdif= Qx-S1*Qxp
Qzdif= Qz-S2*Qzp
S3= INT(SIGN(1.D0,Qzdif))
m= 1   ! Auxiliary to the space profiles (just one point, for the moment)
Beta= VRTeTs(m)/VRNeNs(m)

Funcx_RcdLST= ZL(Qx,ABS(Qz))-Asp*S2*ZS(Qxp,ABS(Qzp)) &
    -Aspp*S3*ZT(Qxdif,ABS(Qzdif))
RETURN
END FUNCTION Funcx_RcdLST

REAL*8 FUNCTION Funcx_RcdLTT(Qxp)
USE Common_Params
USE Common_Arrays
USE Math_Constants
IMPLICIT NONE
REAL*8, INTENT(in) :: Qxp
REAL*8 :: Qx,Qz,Qzp,Qxdif,Qzdif,Beta
INTEGER :: Asp,Aspp,S1,S2,S3 
INTEGER :: m

Qx= Aux1_Rcd(1)
Qz= Aux1_Rcd(2)
Qzp= Aux1_Rcd(3)
Asp= Aux2_Rcd(1)
Aspp= Aux2_Rcd(2)
S1= Aux2_Rcd(3)
S2= Aux2_Rcd(4)
Qxdif= Qx-S1*Qxp
Qzdif= Qz-S2*Qzp
S3= INT(SIGN(1.D0,Qzdif))
m= 1   ! Auxiliary to the space profiles (just one point, for the moment)
Beta= VRTeTs(m)/VRNeNs(m)

Funcx_RcdLTT= ZL(Qx,ABS(Qz))-Asp*S2*ZT(Qxp,ABS(Qzp)) &
    -Aspp*S3*ZT(Qxdif,ABS(Qzdif))
RETURN
END FUNCTION Funcx_RcdLTT

REAL*8 FUNCTION Funcx_RcdSLL(Qxp)
USE Common_Params
USE Common_Arrays
USE Math_Constants
IMPLICIT NONE
REAL*8, INTENT(in) :: Qxp
REAL*8 :: Qx,Qz,Qzp,Qxdif,Qzdif,Beta
INTEGER :: Asp,Aspp,S1,S2,S3 
INTEGER :: m

Qx= Aux1_Rcd(1)
Qz= Aux1_Rcd(2)
Qzp= Aux1_Rcd(3)
Asp= Aux2_Rcd(1)
Aspp= Aux2_Rcd(2)
S1= Aux2_Rcd(3)
S2= Aux2_Rcd(4)
Qxdif= Qx-S1*Qxp
Qzdif= Qz-S2*Qzp
S3= INT(SIGN(1.D0,Qzdif))
m= 1   ! Auxiliary to the space profiles (just one point, for the moment)
Beta= VRTeTs(m)/VRNeNs(m)

Funcx_RcdSLL= ZS(Qx,ABS(Qz))-Asp*S2*ZL(Qxp,ABS(Qzp)) &
    -Aspp*S3*ZL(Qxdif,ABS(Qzdif))
RETURN
END FUNCTION Funcx_RcdSLL

REAL*8 FUNCTION Funcx_RcdSLT(Qxp)
USE Common_Params
USE Common_Arrays
USE Math_Constants
IMPLICIT NONE
REAL*8, INTENT(in) :: Qxp
REAL*8 :: Qx,Qz,Qzp,Qxdif,Qzdif,Beta
INTEGER :: Asp,Aspp,S1,S2,S3 
INTEGER :: m

Qx= Aux1_Rcd(1)
Qz= Aux1_Rcd(2)
Qzp= Aux1_Rcd(3)
Asp= Aux2_Rcd(1)
Aspp= Aux2_Rcd(2)
S1= Aux2_Rcd(3)
S2= Aux2_Rcd(4)
Qxdif= Qx-S1*Qxp
Qzdif= Qz-S2*Qzp
S3= INT(SIGN(1.D0,Qzdif))
m= 1   ! Auxiliary to the space profiles (just one point, for the moment)
Beta= VRTeTs(m)/VRNeNs(m)

Funcx_RcdSLT= ZS(Qx,ABS(Qz))-Asp*S2*ZL(Qxp,ABS(Qzp)) &
    -Aspp*S3*ZT(Qxdif,ABS(Qzdif))
RETURN
END FUNCTION Funcx_RcdSLT

REAL*8 FUNCTION Funcx_RcdTLL(Qxp)
USE Common_Params
USE Common_Arrays
USE Math_Constants
IMPLICIT NONE
REAL*8, INTENT(in) :: Qxp
REAL*8 :: Qx,Qz,Qzp,Qxdif,Qzdif,Beta
INTEGER :: Asp,Aspp,S1,S2,S3 
INTEGER :: m

Qx= Aux1_Rcd(1)
Qz= Aux1_Rcd(2)
Qzp= Aux1_Rcd(3)
Asp= Aux2_Rcd(1)
Aspp= Aux2_Rcd(2)
S1= Aux2_Rcd(3)
S2= Aux2_Rcd(4)
Qxdif= Qx-S1*Qxp
Qzdif= Qz-S2*Qzp
S3= INT(SIGN(1.D0,Qzdif))
m= 1   ! Auxiliary to the space profiles (just one point, for the moment)
Beta= VRTeTs(m)/VRNeNs(m)

Funcx_RcdTLL= ZT(Qx,ABS(Qz))-Asp*S2*ZL(Qxp,ABS(Qzp)) &
    -Aspp*S3*ZL(Qxdif,ABS(Qzdif))
RETURN
END FUNCTION Funcx_RcdTLL

REAL*8 FUNCTION Funcx_RcdTLS(Qxp)
USE Common_Params
USE Common_Arrays
USE Math_Constants
IMPLICIT NONE
REAL*8, INTENT(in) :: Qxp
REAL*8 :: Qx,Qz,Qzp,Qxdif,Qzdif,Beta
INTEGER :: Asp,Aspp,S1,S2,S3 
INTEGER :: m

Qx= Aux1_Rcd(1)
Qz= Aux1_Rcd(2)
Qzp= Aux1_Rcd(3)
Asp= Aux2_Rcd(1)
Aspp= Aux2_Rcd(2)
S1= Aux2_Rcd(3)
S2= Aux2_Rcd(4)
Qxdif= Qx-S1*Qxp
Qzdif= Qz-S2*Qzp
S3= INT(SIGN(1.D0,Qzdif))
m= 1   ! Auxiliary to the space profiles (just one point, for the moment)
Beta= VRTeTs(m)/VRNeNs(m)

Funcx_RcdTLS= ZT(Qx,ABS(Qz))-Asp*S2*ZL(Qxp,ABS(Qzp)) &
    -Aspp*S3*ZS(Qxdif,ABS(Qzdif))
RETURN
END FUNCTION Funcx_RcdTLS

REAL*8 FUNCTION Funcx_RcdTTL(Qxp)
USE Common_Params
USE Common_Arrays
USE Math_Constants
IMPLICIT NONE
REAL*8, INTENT(in) :: Qxp
REAL*8 :: Qx,Qz,Qzp,Qxdif,Qzdif,Beta
INTEGER :: Asp,Aspp,S1,S2,S3 
INTEGER :: m

Qx= Aux1_Rcd(1)
Qz= Aux1_Rcd(2)
Qzp= Aux1_Rcd(3)
Asp= Aux2_Rcd(1)
Aspp= Aux2_Rcd(2)
S1= Aux2_Rcd(3)
S2= Aux2_Rcd(4)
Qxdif= Qx-S1*Qxp
Qzdif= Qz-S2*Qzp
S3= INT(SIGN(1.D0,Qzdif))
m= 1   ! Auxiliary to the space profiles (just one point, for the moment)
Beta= VRTeTs(m)/VRNeNs(m)

Funcx_RcdTTL= ZT(Qx,ABS(Qz))-Asp*S2*ZT(Qxp,ABS(Qzp)) &
    -Aspp*S3*ZL(Qxdif,ABS(Qzdif))
RETURN
END FUNCTION Funcx_RcdTTL

SUBROUTINE Fnorm(Anorm)
USE Common_Params
USE Common_Arrays
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
INTEGER :: l,m
REAL*8, INTENT(out) :: Anorm
REAL*8, DIMENSION(nux) :: Vintux
REAL*8, DIMENSION(nuz) :: Vintuz

 DO l= 1,nux
  FORALL(m=1:nuz)
   Vintuz(m)= Fe(l,m)
  END FORALL
  CALL Simpson(VUz,Vintuz,nuz,Anorm)
  Vintux(l)= Anorm
 END DO
 CALL Simpson(VUx,Vintux,nux,Anorm)
 Anorm= 2.*Anorm  ! Symmetry in Ux.
IF (RenormFe=="Yes") THEN 
 Fe= Fe/Anorm
ELSE
END IF
RETURN
END SUBROUTINE Fnorm

SUBROUTINE Energy(Epart,Ewave,EppEw)
USE Common_Params
USE Common_Arrays
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
INTEGER :: i,k
REAL*8, INTENT(out) :: Ewave,Epart,EppEw
REAL*8 :: Qx2,Uz2,Aux
REAL*8, DIMENSION(nux) :: Vintux
REAL*8, DIMENSION(nuz) :: Vintuz
REAL*8, DIMENSION(nqx) :: Vintqx
REAL*8, DIMENSION(nqz) :: Vintqz

 DO k= 1,nuz
  Uz2= (VUz(k))**2
  FORALL(i=1:nux)
   Vintux(i)= (Uz2+(VUx(i))**2)*Fe(i,k)
  END FORALL
  CALL Simpson(VUx,Vintux,nux,Aux)
  Vintuz(k)= Aux
 END DO
 CALL Simpson(VUz,Vintuz,nuz,Aux)
 Epart= Aux/2.

 DO i= 1,nqx
  Qx2= (VQx(i))**2
  FORALL(k=1:nqz)
   Vintqz(k)= ILp(i,k)+ILm(i,k) &
        + (SQRT(Qx2+(VQz(k))**2))**3*AA/2.*(ISp(i,k)+ISm(i,k)) &
        + ITp(i,k)+ITm(i,k)
  END FORALL
  CALL Simpson(VQz,Vintqz,nqz,Aux)
  Vintqx(i)= Aux
 END DO
 CALL Simpson(VQx,Vintqx,nqx,Aux)
 Ewave= Aux

EppEw= Epart+Ewave

RETURN
END SUBROUTINE Energy

SUBROUTINE Energy2(EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3)
USE Common_Params
USE Common_Arrays
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
INTEGER :: i,k
REAL*8, INTENT(out) :: EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3
REAL*8 :: Qx2,Aux
REAL*8 :: Q,Dq,Qxaux,Qzaux
REAL*8, DIMENSION(nqx) :: Vintqx
REAL*8, DIMENSION(nqz) :: Vintqz
REAL*8, DIMENSION(nph) :: VQ,Vintq,VZtq,VintZq

 DO i= 1,nqx
  FORALL(k=1:nqz)
   Vintqz(k)= ILp(i,k)+ILm(i,k)
  END FORALL
  CALL Simpson(VQz,Vintqz,nqz,Aux)
  Vintqx(i)= Aux
 END DO
 CALL Simpson(VQx,Vintqx,nqx,Aux)
 EwaveL= Aux

 DO i= 1,nqx
  Qx2= (VQx(i))**2
  FORALL(k=1:nqz)
   Vintqz(k)= (SQRT(Qx2+(VQz(k))**2))**3*AA/2.*(ISp(i,k)+ISm(i,k))
  END FORALL
  CALL Simpson(VQz,Vintqz,nqz,Aux)
  Vintqx(i)= Aux
 END DO
 CALL Simpson(VQx,Vintqx,nqx,Aux)
 EwaveS= Aux

 DO i= 1,nqx
  FORALL(k=1:nqz)
   Vintqz(k)= ITp(i,k)+ITm(i,k)
  END FORALL
  CALL Simpson(VQz,Vintqz,nqz,Aux)
  Vintqx(i)= Aux
 END DO
 CALL Simpson(VQx,Vintqx,nqx,Aux)
 EwaveT= Aux

 Dq= (Qzf-Qmin)/(nph-1)
 DO i= 1,nph
  Q= Qmin+(i-1)*Dq
  VQ(i)= Q
  VZtq(i)= SQRT(1.+Q**2/Ve2C2)
 END DO
 DO i= 1,nph
  Q= VQ(i)
  DO k= 1,nph
   Qxaux= Q*SIN(VPhip(k))
   Qzaux= Q*COS(VPhip(k))
   CALL Aitp2d(nqx,nqz,VQx,VQz,ITp,Qxaux,Qzaux,Aux)
   Vintq(k)= Aux
   CALL Aitp2d(nqx,nqz,VQx,VQz,ITm,Qxaux,Qzaux,Aux)
   Vintq(k)= Vintq(k)+Aux
  END DO 
  CALL Simpson(VPhip,Vintq,nph,Aux)
  VintZq(i)= Aux*VQ(i)
 END DO
 EwaveTF= 0.
 EwaveTH= 0.
 EwaveT3= 0.
 DO i= 1,nph-1
  IF (VZtq(i)<1.5) THEN
   EwaveTF= EwaveTF+Dq*(VintZq(i)+VintZq(i+1))/2.   
  ELSE
   IF (VZtq(i)<2.5) THEN
    EwaveTH= EwaveTH+Dq*(VintZq(i)+VintZq(i+1))/2.   
   ELSE
    IF (VZtq(i)<3.5) THEN
     EwaveT3= EwaveT3+Dq*(VintZq(i)+VintZq(i+1))/2.   
    ELSE
    END IF
   END IF
  END IF
 END DO

RETURN
END SUBROUTINE Energy2

SUBROUTINE Output_Coef(CoefChoice,CoefA,CoefB)
USE Common_Params
USE Common_Arrays
USE Math_Constants
IMPLICIT NONE
INTEGER :: i,k
INTEGER, PARAMETER :: Step=1
REAL*8 :: Qx,Qz
REAL*8, DIMENSION(nqx,nqz), INTENT(in) :: CoefA,CoefB
CHARACTER(LEN=5), INTENT(in) :: CoefChoice

SELECT CASE(CoefChoice)
 CASE("Lwave")
  OPEN(1,FILE='Lwave.wt')
 CASE("Swave")
  OPEN(1,FILE='Swave.wt')
 CASE("Twave")
  OPEN(1,FILE='Twave.wt')
 CASE DEFAULT
END SELECT

  DO i= 1,nqx,Step
   Qx= VQx(nqx+1-i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    WRITE(1,*) -Qx,-Qz,CoefA(nqx+1-i,nqz+1-k),CoefB(nqx+1-i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    WRITE(1,*) -Qx,Qz,CoefA(nqx+1-i,k),CoefB(nqx+1-i,k)
   END DO
  END DO
  DO i= 1,nqx,Step
   Qx= VQx(i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    WRITE(1,*) Qx,-Qz,CoefA(i,nqz+1-k),CoefB(i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    WRITE(1,*) Qx,Qz,CoefA(i,k),CoefB(i,k)
   END DO
  END DO
CLOSE(1)

RETURN
END SUBROUTINE Output_Coef

SUBROUTINE Output(WriteChoice)
USE Common_Params
USE Common_Arrays
USE Math_Constants
IMPLICIT NONE
INTEGER :: i,k
INTEGER, PARAMETER :: Step=1
REAL*8 :: Aux
REAL*8 :: Qx,Qz,Q,Muq
REAL*8, DIMENSION(nqx) :: Vintqx
REAL*8, DIMENSION(nqz) :: Iqzm,Iqzp,Iqzm2,Iqzp2
REAL*8, DIMENSION(nux) :: Vintux
REAL*8, DIMENSION(nuz) :: Feuz
CHARACTER(LEN=3), INTENT(in) :: WriteChoice

SELECT CASE(WriteChoice)

 CASE("Fe0")
  OPEN(1,FILE='Fe0.wt')
  DO i= 1,nux,Step
   WRITE(1,*)' '
   DO k= 1,nuz,Step
    IF(ABS(Fe0(nux+1-i,k))<=EpsMin) THEN
     Aux= EpsMin
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
     Aux= EpsMin
    ELSE
     Aux= Fe0(i,k)
    END IF
    WRITE(1,*) VUx(i),VUz(k),Aux
   END DO
  END DO
  CLOSE(1)
  OPEN(1,FILE='Fe1D0.wt')
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

 CASE("Fe ")
  OPEN(1,FILE='Fe.wt')
  DO i= 1,nux,Step
   WRITE(1,*)' '
   DO k= 1,nuz,Step
    IF(ABS(Fe(nux+1-i,k))<=EpsMin) THEN
     Aux= EpsMin
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
     Aux= EpsMin
    ELSE
     Aux= Fe(i,k)
    END IF
    WRITE(1,*) VUx(i),VUz(k),Aux
   END DO
  END DO
  CLOSE(1)

  CASE("DDz")
  OPEN(1,FILE='DDzzduz.wt')
  DO i= 1,nux,Step
   WRITE(1,*)' '
   DO k= 1,nuz,Step
    WRITE(1,*) -VUx(nux+1-i),VUz(k), DDzzduz(nux+1-i, k)
   END DO
  END DO
  DO i= 1,nux,Step
   WRITE(1,*)' '
   DO k= 1,nuz,Step
    WRITE(1,*) VUx(i),VUz(k),DDzzduz(i,k)
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
   Qx= VQx(nqx+1-i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    ! IF(ILm(nqx+1-i, nqz+1-k) .le. 4.6E-006) THEN
    !  ILm(nqx+1-i, nqz+1-k) = 4.6E-006
    ! END IF
    WRITE(1,*) -Qx,-Qz,ILm(nqx+1-i,nqz+1-k),ZL(-Qx,-Qz)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    ! IF(ILp(nqx+1-i, k) .le. 4.6E-006) THEN
    !  ILp(nqx+1-i, k) = 4.6E-006
    ! END IF
    WRITE(1,*) -Qx,Qz,ILp(nqx+1-i,k),ZL(-Qx,Qz)
   END DO
  END DO
  DO i= 1,nqx,Step
   Qx= VQx(i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    ! IF(ILm(i, nqz+1-k) .le. 4.6E-006) THEN
    !  ILm(i, nqz+1-k) = 4.6E-006
    ! END IF
    WRITE(1,*) Qx,-Qz,ILm(i,nqz+1-k),ZL(Qx,-Qz)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    ! IF(ILp(i, k) .le. 4.6E-006) THEN
    !  ILp(i, k) = 4.6E-006
    ! END IF
    WRITE(1,*) Qx,Qz,ILp(i,k),ZL(Qx,Qz)
   END DO
  END DO
  CLOSE(1)


 CASE("IS ")
  OPEN(1,FILE='IS.wt')
  DO i= 1,nqx,Step
   Qx= VQx(nqx+1-i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    Q= SQRT(Qx**2+Qz**2)
    Muq= Q**3*AA/2.
    WRITE(1,*) -Qx,-Qz,Muq*ISm(nqx+1-i,nqz+1-k),ZS(-Qx,-Qz)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    Q= SQRT(Qx**2+Qz**2)
    Muq= Q**3*AA/2.
    WRITE(1,*) -Qx,Qz,Muq*ISp(nqx+1-i,k),ZS(-Qx,Qz)
   END DO
  END DO
  DO i= 1,nqx,Step
   Qx= VQx(i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    Q= SQRT(Qx**2+Qz**2)
    Muq= Q**3*AA/2.
    WRITE(1,*) Qx,-Qz,Muq*ISm(i,nqz+1-k),ZS(Qx,-Qz)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    Q= SQRT(Qx**2+Qz**2)
    Muq= Q**3*AA/2.
    WRITE(1,*) Qx,Qz,Muq*ISp(i,k),ZS(Qx,Qz)
   END DO
  END DO
  CLOSE(1)

 CASE("IT ")
  OPEN(1,FILE='IT.wt')
  DO i= 1,nqx,Step
   Qx= VQx(nqx+1-i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    WRITE(1,*) -Qx,-Qz,ITm(nqx+1-i,nqz+1-k),ZT(-Qx,-Qz)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    WRITE(1,*) -Qx,Qz,ITp(nqx+1-i,k),ZT(-Qx,Qz)
   END DO
  END DO
  DO i= 1,nqx,Step
   Qx= VQx(i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    WRITE(1,*) Qx,-Qz,ITm(i,nqz+1-k),ZT(Qx,-Qz)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    WRITE(1,*) Qx,Qz,ITp(i,k),ZT(Qx,Qz)
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
  OPEN(1,FILE='IL0Q.wt')
  DO k= 1,nqz
   WRITE(1,*) VQz(k),ILp(1,k)
  END DO
  CLOSE(1)
  OPEN(1,FILE='IL1D0.wt')
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
    IF(Q<=5.E-3) THEN
     Muq= (SQRT(Qx**2+(VQz(nqz+1-k+Step))**2))**3*AA/2.
     WRITE(1,*) Qx,-Qz,Muq*ISm(nqx+1-i,nqz+1-k+Step)
    ELSE
     Muq= Q**3*AA/2.
     WRITE(1,*) Qx,-Qz,Muq*ISm(nqx+1-i,nqz+1-k)
    END IF
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    Q= SQRT(Qx**2+Qz**2)
    IF(Q<=5.E-3) THEN
     Muq= (SQRT(Qx**2+(VQz(k+Step))**2))**3*AA/2.
     WRITE(1,*) Qx,Qz,Muq*ISp(nqx+1-i,k+Step)
    ELSE
     Muq= Q**3*AA/2.
     WRITE(1,*) Qx,Qz,Muq*ISp(nqx+1-i,k)
    END IF
   END DO
  END DO
  DO i= 1,nqx,Step
   Qx= VQx(i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    Q= SQRT(Qx**2+Qz**2)
    IF(Q<=5.E-3) THEN
     Muq= (SQRT(Qx**2+(VQz(nqz+1-k+Step))**2))**3*AA/2.
     WRITE(1,*) Qx,-Qz,Muq*ISm(i,nqz+1-k+Step)
    ELSE
     Muq= Q**3*AA/2.
     WRITE(1,*) Qx,-Qz,Muq*ISm(i,nqz+1-k)
    END IF
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    Q= SQRT(Qx**2+Qz**2)
    IF(Q<=5.E-3) THEN
     Muq= (SQRT(Qx**2+(VQz(k+Step))**2))**3*AA/2.
     WRITE(1,*) Qx,Qz,Muq*ISp(i,k+Step)
    ELSE
     Muq= Q**3*AA/2.
     WRITE(1,*) Qx,Qz,Muq*ISp(i,k)
    END IF
   END DO
  END DO
  CLOSE(1)
  OPEN(1,FILE='IS0Q.wt')
  DO k= 1,nqz
   Qz= VQz(k)
   IF(Qz<=5.E-3) THEN
    Muq= (SQRT((VQz(k+Step))**2))**3*AA/2.
    WRITE(1,*) Qz,Muq*ISp(1,k+Step)
   ELSE
    Muq= Qz**3*AA/2.
    WRITE(1,*) Qz,Muq*ISp(1,k)
   END IF
  END DO
  CLOSE(1)
  OPEN(1,FILE='IS1D0.wt')
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
   WRITE(1,*) -VQz(nqz+1-k),Iqzm2(nqz+1-k),Iqzm(nqz+1-k)
  END DO
  DO k= 1,nqz,Step
   WRITE(1,*) VQz(k),Iqzp2(k),Iqzp(k)
  END DO
  CLOSE(1)

 CASE("IT0")
  OPEN(1,FILE='IT0.wt')
  DO i= 1,nqx,Step
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    WRITE(1,*) -VQx(nqx+1-i),-VQz(nqz+1-k),ITm(nqx+1-i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    WRITE(1,*) -VQx(nqx+1-i),VQz(k),ITp(nqx+1-i,k)
   END DO
  END DO
  DO i= 1,nqx,Step
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    WRITE(1,*) VQx(i),-VQz(nqz+1-k),ITm(i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    WRITE(1,*) VQx(i),VQz(k),ITp(i,k)
   END DO
  END DO
  CLOSE(1)
  OPEN(1,FILE='IT0Q.wt')
  DO k= 1,nqz
   WRITE(1,*) VQz(k),ITp(1,k)
  END DO
  CLOSE(1)
  OPEN(1,FILE='IT1D0.wt')
  DO k= 1,nqz
   DO i= 1,nqx
    Vintqx(i)= ITm(i,k)
   END DO
   CALL Simpson(VQx,Vintqx,nqx,Aux)
   Iqzm(k)= Aux
   DO i= 1,nqx
    Vintqx(i)= ITp(i,k)
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
   WRITE(1,*) -VQz(nqz+1-k),Iqzm2(nqz+1-k),Iqzm(nqz+1-k)
  END DO
  DO k= 1,nqz,Step
   WRITE(1,*) VQz(k),Iqzp2(k),Iqzp(k)
  END DO
  CLOSE(1)

 CASE("IT1")
  OPEN(1,FILE='IT1D.wt')
  DO k= 1,nqz
   DO i= 1,nqx
    Vintqx(i)= ITm(i,k)
   END DO
   CALL Simpson(VQx,Vintqx,nqx,Aux)
   Iqzm(k)= Aux
   DO i= 1,nqx
    Vintqx(i)= ITp(i,k)
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

 CASE("GcL")
  OPEN(1,FILE='GcollL.wt')
  DO i= 1,nqx,Step
   Qx= VQx(nqx+1-i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    WRITE(1,*) -Qx,-Qz,GcollLm(nqx+1-i,nqz+1-k),-GcollLm(nqx+1-i,nqz+1-k),&
     GcollLm(nqx+1-i,nqz+1-k)/(GqlLm(nqx+1-i,nqz+1-k)+EpsMin)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    WRITE(1,*) -Qx,Qz,GcollLp(nqx+1-i,k),-GcollLp(nqx+1-i,k),&
     GcollLp(nqx+1-i,k)/(GqlLp(nqx+1-i,k)+EpsMin)
   END DO
  END DO
  DO i= 1,nqx,Step
   Qx= VQx(i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    WRITE(1,*) Qx,-Qz,GcollLm(i,nqz+1-k),-GcollLm(i,nqz+1-k),&
     GcollLm(i,nqz+1-k)/(GqlLm(i,nqz+1-k)+EpsMin)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    WRITE(1,*) Qx,Qz,GcollLp(i,k),-GcollLp(i,k),&
     GcollLp(i,k)/(GqlLp(i,k)+EpsMin)
   END DO
  END DO
  CLOSE(1)
  OPEN(1,FILE='GcollL1D.wt')
  DO i= 1,nqmod
   WRITE(1,*) VQQ(i),GcollL1D(i)
  END DO
  CLOSE(1)

 CASE("GqL")
  OPEN(1,FILE='GqlL.wt')
  DO i= 1,nqx,Step
   Qx= VQx(nqx+1-i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    WRITE(1,*) -Qx,-Qz,GqlLm(nqx+1-i,nqz+1-k),-GqlLm(nqx+1-i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    WRITE(1,*) -Qx,Qz,GqlLp(nqx+1-i,k),-GqlLp(nqx+1-i,k)
   END DO
  END DO
  DO i= 1,nqx,Step
   Qx= VQx(i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    WRITE(1,*) Qx,-Qz,GqlLm(i,nqz+1-k),-GqlLm(i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    WRITE(1,*) Qx,Qz,GqlLp(i,k),-GqlLp(i,k)
   END DO
  END DO
  CLOSE(1)
  OPEN(1,FILE='GqlL1D.wt')
  DO i= 1,nqmod
   WRITE(1,*) VQQ(i),GqlL1D(i)
  END DO
  CLOSE(1)

 CASE("GcS")
  OPEN(1,FILE='GcollS.wt')
  DO i= 1,nqx,Step
   Qx= VQx(nqx+1-i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    WRITE(1,*) -Qx,-Qz,GcollSm(nqx+1-i,nqz+1-k),-GcollSm(nqx+1-i,nqz+1-k),&
     GcollSm(nqx+1-i,nqz+1-k)/(GqlSm(nqx+1-i,nqz+1-k)+EpsMin)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    WRITE(1,*) -Qx,Qz,GcollSp(nqx+1-i,k),-GcollSp(nqx+1-i,k),&
     GcollSp(nqx+1-i,k)/(GqlSp(nqx+1-i,k)+EpsMin)
   END DO
  END DO
  DO i= 1,nqx,Step
   Qx= VQx(i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    WRITE(1,*) Qx,-Qz,GcollSm(i,nqz+1-k),-GcollSm(i,nqz+1-k),&
     GcollSm(i,nqz+1-k)/(GqlSm(i,nqz+1-k)+EpsMin)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    WRITE(1,*) Qx,Qz,GcollSp(i,k),-GcollSp(i,k),&
     GcollSp(i,k)/(GqlSp(i,k)+EpsMin)
   END DO
  END DO
  CLOSE(1)
  OPEN(1,FILE='GcollS1D.wt')
  DO i= 1,nqmod
   WRITE(1,*) VQQ(i),GcollS1D(i)
  END DO
  CLOSE(1)

 CASE("GqS")
  OPEN(1,FILE='GqlS.wt')
  DO i= 1,nqx,Step
   Qx= VQx(nqx+1-i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    WRITE(1,*) -Qx,-Qz,GqlSm(nqx+1-i,nqz+1-k),-GqlSm(nqx+1-i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    WRITE(1,*) -Qx,Qz,GqlSp(nqx+1-i,k),-GqlSp(nqx+1-i,k)
   END DO
  END DO
  DO i= 1,nqx,Step
   Qx= VQx(i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    WRITE(1,*) Qx,-Qz,GqlSm(i,nqz+1-k),-GqlSm(i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    WRITE(1,*) Qx,Qz,GqlSp(i,k),-GqlSp(i,k)
   END DO
  END DO
  CLOSE(1)
  OPEN(1,FILE='GqlS1D.wt')
  DO i= 1,nqmod
   WRITE(1,*) VQQ(i),GqlS1D(i)
  END DO
  CLOSE(1)

 CASE("PbL")
  OPEN(1,FILE='PbrL.wt')
  DO i= 1,nqx,Step
   Qx= VQx(nqx+1-i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    WRITE(1,*) -Qx,-Qz,BremL(nqx+1-i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    WRITE(1,*) -Qx,Qz,BremL(nqx+1-i,k)
   END DO
  END DO
  DO i= 1,nqx,Step
   Qx= VQx(i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    WRITE(1,*) Qx,-Qz,BremL(i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    WRITE(1,*) Qx,Qz,BremL(i,k)
   END DO
  END DO
  CLOSE(1)
  OPEN(1,FILE='PbrL1D.wt')
  DO i= 1,nqmod
   WRITE(1,*) VQQ(i),BremL1D(i)
  END DO
  CLOSE(1)

 CASE("PbS")
  OPEN(1,FILE='PbrS.wt')
  DO i= 1,nqx,Step
   Qx= VQx(nqx+1-i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    WRITE(1,*) -Qx,-Qz,BremS(nqx+1-i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    WRITE(1,*) -Qx,Qz,BremS(nqx+1-i,k)
   END DO
  END DO
  DO i= 1,nqx,Step
   Qx= VQx(i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    WRITE(1,*) Qx,-Qz,BremS(i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    WRITE(1,*) Qx,Qz,BremS(i,k)
   END DO
  END DO
  CLOSE(1)
  OPEN(1,FILE='PbrS1D.wt')
  DO i= 1,nqmod
   WRITE(1,*) VQQ(i),BremS1D(i)
  END DO
  CLOSE(1)

  CASE DEFAULT

END SELECT

RETURN

END SUBROUTINE Output

SUBROUTINE Output2(WriteChoice)
! Subroutine generating output in the format suitable for MatLab (Nov/2013)
! Nov. 28, 2013
USE Common_Params
USE Common_Arrays
USE Math_Constants
IMPLICIT NONE
INTEGER :: i,k
INTEGER, PARAMETER :: Step=1
REAL*8 :: Aux
REAL*8 :: Qx,Qz,Q,Muq
REAL*8 :: Dq,Qxaux,Qzaux
REAL*8, DIMENSION(nqx) :: Vintqx
REAL*8, DIMENSION(nqz) :: Iqzm,Iqzp,Iqzm2,Iqzp2
REAL*8, DIMENSION(nux) :: Vintux
REAL*8, DIMENSION(nuz) :: Feuz
REAL*8, DIMENSION(nph) :: VQ,Vintq
CHARACTER(LEN=3), INTENT(in) :: WriteChoice

SELECT CASE(WriteChoice)

 CASE("Ui ")
  OPEN(1,FILE='Ux')
  DO i= 1,nux,Step
   WRITE(1,*) -VUx(nux+1-i)
  END DO
  DO i= 1,nux,Step
   WRITE(1,*) VUx(i)
  END DO
  CLOSE(1)
  OPEN(1,FILE='Uz')
  DO i= 1,nuz,Step
   WRITE(1,*) VUz(i)
  END DO
  CLOSE(1)

 CASE("Qi ")
  OPEN(1,FILE='Qx')
  DO i= 1,nqx,Step
   WRITE(1,*) -VQx(nqx+1-i)
  END DO
  DO i= 1,nqx,Step
   WRITE(1,*) VQx(i)
  END DO
  CLOSE(1)
  OPEN(1,FILE='Qz')
  DO i= 1,nqz,Step
   WRITE(1,*) -VQz(nqz+1-i)
  END DO
  DO i= 1,nqz,Step
   WRITE(1,*) VQz(i)
  END DO
  CLOSE(1)

 CASE("Q  ")
  OPEN(1,FILE='Q')
  Dq= (Qzf-Qmin)/(nph-1)
  DO i= 1,nph
   Q= Qmin+(i-1)*Dq
   WRITE(1,*) Q
  END DO
  CLOSE(1)
  OPEN(1,FILE='ZTQ')
  DO i= 1,nph
   Q= Qmin+(i-1)*Dq
   WRITE(1,*) SQRT(1.+Q**2/Ve2C2)
  END DO
  CLOSE(1)

 CASE("Fe0")
  OPEN(1,FILE='Fe02D')
  DO i= 1,nux,Step
   WRITE(1,*)' '
   DO k= 1,nuz,Step
    IF(ABS(Fe0(nux+1-i,k))<=EpsMin) THEN
     Aux= EpsMin
    ELSE
     Aux= Fe0(nux+1-i,k)
    END IF
    WRITE(1,*) Aux
   END DO
  END DO
  DO i= 1,nux,Step
   WRITE(1,*)' '
   DO k= 1,nuz,Step
    IF(ABS(Fe0(i,k))<=EpsMin) THEN
     Aux= EpsMin
    ELSE
     Aux= Fe0(i,k)
    END IF
    WRITE(1,*) Aux
   END DO
  END DO
  CLOSE(1)
  OPEN(1,FILE='Fe01D')
  DO k= 1,nuz
   DO i= 1,nux
    Vintux(i)= Fe(i,k)
   END DO
   CALL Simpson(VUx,Vintux,nux,Aux)
   Feuz(k)= Aux
  END DO
  DO k= 1,nuz,Step
   WRITE(1,*) Feuz(k)
  END DO
  CLOSE(1)

 CASE("Fe ")
  OPEN(1,FILE='Fe2D')
  DO i= 1,nux,Step
   WRITE(1,*)' '
   DO k= 1,nuz,Step
    IF(ABS(Fe(nux+1-i,k))<=EpsMin) THEN
     Aux= EpsMin
    ELSE
     Aux= Fe(nux+1-i,k)
    END IF
    WRITE(1,*) Aux
   END DO
  END DO
  DO i= 1,nux,Step
   WRITE(1,*)' '
   DO k= 1,nuz,Step
    IF(ABS(Fe(i,k))<=EpsMin) THEN
     Aux= EpsMin
    ELSE
     Aux= Fe(i,k)
    END IF
    WRITE(1,*) Aux
   END DO
  END DO
  CLOSE(1)

 CASE("IL ")
  OPEN(1,FILE='IL2D')
  DO i= 1,nqx,Step
   Qx= VQx(nqx+1-i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    WRITE(1,*) ILm(nqx+1-i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    WRITE(1,*) ILp(nqx+1-i,k)
   END DO
  END DO
  DO i= 1,nqx,Step
   Qx= VQx(i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    WRITE(1,*) ILm(i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    WRITE(1,*) ILp(i,k)
   END DO
  END DO
  CLOSE(1)

 CASE("IS ")
  OPEN(1,FILE='IS2D')
  DO i= 1,nqx,Step
   Qx= VQx(nqx+1-i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    Q= SQRT(Qx**2+Qz**2)
    Muq= Q**3*AA/2.
    WRITE(1,*) Muq*ISm(nqx+1-i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    Q= SQRT(Qx**2+Qz**2)
    Muq= Q**3*AA/2.
    WRITE(1,*) Muq*ISp(nqx+1-i,k)
   END DO
  END DO
  DO i= 1,nqx,Step
   Qx= VQx(i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    Q= SQRT(Qx**2+Qz**2)
    Muq= Q**3*AA/2.
    WRITE(1,*) Muq*ISm(i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    Q= SQRT(Qx**2+Qz**2)
    Muq= Q**3*AA/2.
    WRITE(1,*) Muq*ISp(i,k)
   END DO
  END DO
  CLOSE(1)

 CASE("IT ")
  OPEN(1,FILE='IT2D')
  DO i= 1,nqx,Step
   Qx= VQx(nqx+1-i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    WRITE(1,*) ITm(nqx+1-i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    WRITE(1,*) ITp(nqx+1-i,k)
   END DO
  END DO
  DO i= 1,nqx,Step
   Qx= VQx(i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    WRITE(1,*) ITm(i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    WRITE(1,*) ITp(i,k)
   END DO
  END DO
  CLOSE(1)

 CASE("IL0")
  OPEN(1,FILE='IL02D')
  DO i= 1,nqx,Step
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    WRITE(1,*) ILm(nqx+1-i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    WRITE(1,*) ILp(nqx+1-i,k)
   END DO
  END DO
  DO i= 1,nqx,Step
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    WRITE(1,*) ILm(i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    WRITE(1,*) ILp(i,k)
   END DO
  END DO
  CLOSE(1)
  OPEN(1,FILE='IL01D')
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
   WRITE(1,*) Iqzm(nqz+1-k)
  END DO
  DO k= 1,nqz,Step
   WRITE(1,*) Iqzp(k)
  END DO
  CLOSE(1)

 CASE("IS0")
  OPEN(1,FILE='IS02D')
  DO i= 1,nqx,Step
   Qx= -VQx(nqx+1-i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    Q= SQRT(Qx**2+Qz**2)
    Muq= Q**3*AA/2.
    WRITE(1,*) Muq*ISm(nqx+1-i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    Q= SQRT(Qx**2+Qz**2)
    Muq= Q**3*AA/2.
    WRITE(1,*) Muq*ISp(nqx+1-i,k)
   END DO
  END DO
  DO i= 1,nqx,Step
   Qx= VQx(i)
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    Qz= VQz(nqz+1-k)
    Q= SQRT(Qx**2+Qz**2)
    Muq= Q**3*AA/2.
    WRITE(1,*) Muq*ISm(i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    Qz= VQz(k)
    Q= SQRT(Qx**2+Qz**2)
    Muq= Q**3*AA/2.
    WRITE(1,*) Muq*ISp(i,k)
   END DO
  END DO
  CLOSE(1)
  OPEN(1,FILE='IS01D')
  DO k= 1,nqz
   Qz= VQz(k)
   DO i= 1,nqx
    Qx= VQx(i)
    Q= SQRT(Qx**2+Qz**2)
    Muq= Q**3*AA/2.
    Vintqx(i)= Muq*ISm(i,k)
   END DO
   CALL Simpson(VQx,Vintqx,nqx,Aux)
   Iqzm2(k)= Aux
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
   WRITE(1,*) Iqzm2(nqz+1-k)
  END DO
  DO k= 1,nqz,Step
   WRITE(1,*) Iqzp2(k)
  END DO
  CLOSE(1)

 CASE("IT0")
  OPEN(1,FILE='IT02D')
  DO i= 1,nqx,Step
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    WRITE(1,*) ITm(nqx+1-i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    WRITE(1,*) ITp(nqx+1-i,k)
   END DO
  END DO
  DO i= 1,nqx,Step
   WRITE(1,*)' '
   DO k= 1,nqz,Step
    WRITE(1,*) ITm(i,nqz+1-k)
   END DO
   DO k= 1,nqz,Step
    WRITE(1,*) ITp(i,k)
   END DO
  END DO
  CLOSE(1)
  OPEN(1,FILE='IT01D')
  DO k= 1,nqz
   DO i= 1,nqx
    Vintqx(i)= ITm(i,k)
   END DO
   CALL Simpson(VQx,Vintqx,nqx,Aux)
   Iqzm(k)= Aux
   DO i= 1,nqx
    Vintqx(i)= ITp(i,k)
   END DO
   CALL Simpson(VQx,Vintqx,nqx,Aux)
   Iqzp(k)= Aux
  END DO
  DO k= 1,nqz,Step
   WRITE(1,*) Iqzm(nqz+1-k)
  END DO
  DO k= 1,nqz,Step
   WRITE(1,*) Iqzp(k)
  END DO
  CLOSE(1)

 CASE("Fe1")
  OPEN(1,FILE='Fe1D')
  DO k= 1,nuz
   DO i= 1,nux
    Vintux(i)= Fe(i,k)
   END DO
   CALL Simpson(VUx,Vintux,nux,Aux)
   Feuz(k)= Aux
  END DO
  DO k= 1,nuz,Step
   WRITE(1,*) Feuz(k)
  END DO
  CLOSE(1)

 CASE("IL1")
  OPEN(1,FILE='IL1D')
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
   WRITE(1,*) Iqzm(nqz+1-k)
  END DO
  DO k= 1,nqz,Step
   WRITE(1,*) Iqzp(k)
  END DO
  CLOSE(1)

 CASE("IS1")
  OPEN(1,FILE='IS1D')
  DO k= 1,nqz
   Qz= VQz(k)
   DO i= 1,nqx
    Qx= VQx(i)
    Q= SQRT(Qx**2+Qz**2)
    Muq= Q**3*AA/2.
    Vintqx(i)= Muq*ISm(i,k)
   END DO
   CALL Simpson(VQx,Vintqx,nqx,Aux)
   Iqzm2(k)= Aux
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
   WRITE(1,*) Iqzm2(nqz+1-k)
  END DO
  DO k= 1,nqz,Step
   WRITE(1,*) Iqzp2(k)
  END DO
  CLOSE(1)

 CASE("IT1")
  OPEN(1,FILE='IT1D')
  DO k= 1,nqz
   DO i= 1,nqx
    Vintqx(i)= ITm(i,k)
   END DO
   CALL Simpson(VQx,Vintqx,nqx,Aux)
   Iqzm(k)= Aux
   DO i= 1,nqx
    Vintqx(i)= ITp(i,k)
   END DO
   CALL Simpson(VQx,Vintqx,nqx,Aux)
   Iqzp(k)= Aux
  END DO
  DO k= 1,nqz,Step
   WRITE(1,*) Iqzm(nqz+1-k)
  END DO
  DO k= 1,nqz,Step
   WRITE(1,*) Iqzp(k)
  END DO
  CLOSE(1)

 CASE("ILQ")
  OPEN(1,FILE='ILQ1')
  OPEN(2,FILE='ILQ2')
  Dq= (Qzf-Qmin)/(nph-1)
  DO i= 1,nph
   Q= Qmin+(i-1)*Dq
   VQ(i)= Q
  END DO
  DO i= 1,nph
   Q= VQ(i)
   DO k= 1,nph
    Qxaux= Q*SIN(VPhip(k))
    Qzaux= Q*COS(VPhip(k))
    CALL Aitp2d(nqx,nqz,VQx,VQz,ILp,Qxaux,Qzaux,Aux)
    Vintq(k)= Aux
    CALL Aitp2d(nqx,nqz,VQx,VQz,ILm,Qxaux,Qzaux,Aux)
    Vintq(k)= Vintq(k)+Aux
   END DO 
   CALL Simpson(VPhip,Vintq,nph,Aux)
   WRITE(1,*) Aux
   WRITE(2,*) Aux*Q
  END DO
  CLOSE(1)
  CLOSE(2)

 CASE("ISQ")
  OPEN(1,FILE='ISQ1')
  OPEN(2,FILE='ISQ2')
  Dq= (Qzf-Qmin)/(nph-1)
  DO i= 1,nph
   Q= Qmin+(i-1)*Dq
   VQ(i)= Q
  END DO
  DO i= 1,nph
   Q= VQ(i)
   DO k= 1,nph
    Qxaux= Q*SIN(VPhip(k))
    Qzaux= Q*COS(VPhip(k))
    CALL Aitp2d(nqx,nqz,VQx,VQz,ISp,Qxaux,Qzaux,Aux)
    Vintq(k)= Aux
    CALL Aitp2d(nqx,nqz,VQx,VQz,ISm,Qxaux,Qzaux,Aux)
    Vintq(k)= Vintq(k)+Aux
   END DO 
   CALL Simpson(VPhip,Vintq,nph,Aux)
   WRITE(1,*) Aux
   WRITE(2,*) Aux*Q
  END DO
  CLOSE(1)
  CLOSE(2)

 CASE("ITQ")
  OPEN(1,FILE='ITQ1')
  OPEN(2,FILE='ITQ2')
  Dq= (Qzf-Qmin)/(nph-1)
  DO i= 1,nph
   Q= Qmin+(i-1)*Dq
   VQ(i)= Q
  END DO
  DO i= 1,nph
   Q= VQ(i)
   DO k= 1,nph
    Qxaux= Q*SIN(VPhip(k))
    Qzaux= Q*COS(VPhip(k))
    CALL Aitp2d(nqx,nqz,VQx,VQz,ITp,Qxaux,Qzaux,Aux)
    Vintq(k)= Aux
    CALL Aitp2d(nqx,nqz,VQx,VQz,ITm,Qxaux,Qzaux,Aux)
    Vintq(k)= Vintq(k)+Aux
   END DO 
   CALL Simpson(VPhip,Vintq,nph,Aux)
   WRITE(1,*) Aux
   WRITE(2,*) Aux*Q
  END DO
  CLOSE(1)
  CLOSE(2)

  CASE DEFAULT

END SELECT

RETURN

END SUBROUTINE Output2

!---------------------------------------------------------------------
! Mathematical Routines:
!---------------------------------------------------------------------

SUBROUTINE Split(Dux,Duz,DTau,it)

! See Notes D-10.

USE Common_Params
USE Common_Arrays
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
INTEGER, INTENT(in) :: it
REAL*8 :: Auxx,Auxz,Auxxx,Auxzz
REAL*8, INTENT(in) :: DTau,Dux,Duz
REAL*8, DIMENSION(nux) :: Alpha1,Beta1,Gamma1,Psi1
REAL*8, DIMENSION(nuz) :: Alpha2,Beta2,Gamma2,Psi2
REAL*8, DIMENSION(nux) :: Faux1,Faux1b
REAL*8, DIMENSION(nuz) :: Faux2,Faux2b
REAL*8, DIMENSION(nux) :: Baux1,Gaux1
REAL*8, DIMENSION(nuz) :: Baux2,Gaux2
REAL*8, DIMENSION(nux,nuz) :: Dfdux,Dfduz
REAL*8, DIMENSION(nux,nuz) :: Dfduxdux,Dfduxduz
REAL*8, DIMENSION(nux,nuz) :: Dfduzdux,Dfduzduz
REAL*8, DIMENSION(nux,nuz) :: DAxdux
REAL*8, DIMENSION(nux,nuz) :: DAzduz
REAL*8, DIMENSION(nux,nuz) :: DDxxdux,DDxzdux
REAL*8, DIMENSION(nux,nuz) :: DDzxduz!,DDzzduz
REAL*8, DIMENSION(nux,nuz) :: array
REAL*8, DIMENSION(nux,nuz) :: CAx,CAz
REAL*8, DIMENSION(nux,nuz) :: CDxx,CDxz,CDzx,CDzz
REAL*8, DIMENSION(nux) :: Auxi
REAL*8, DIMENSION(nuz) :: Auxk
REAL*8, DIMENSION(nux,nuz) :: Fex,Fez
INTEGER :: i,k

! Limitador
! FeOld = Fe

! CALL Coef_A
CALL Coef_D
Auxx= DTau/4./Dux
Auxz= DTau/4./Duz
Auxxx= DTau/2./Dux/Dux
Auxzz= DTau/2./Duz/Duz

CAx= Ax+ColAx
CAz= Az+ColAz
CDxx= Dxx+ColDxx
CDxz= Dxz+ColDxz
CDzx= Dzx+ColDxz
CDzz= Dzz+ColDzz

CALL Derivx5p2d(nux,nuz,VUx,CDxx,DDxxdux)
CALL Derivx5p2d(nux,nuz,VUx,CDxz,DDxzdux)
CALL Derivy5p2d(nux,nuz,VUz,CDzx,DDzxduz)
CALL Derivy5p2d(nux,nuz,VUz,CDzz,DDzzduz)
CALL Derivx5p2d(nux,nuz,VUx,CAx,DAxdux)
CALL Derivy5p2d(nux,nuz,VUz,CAz,DAzduz)

! CALL BezierDDzzduz
! CALL BezierDDzzduz


SELECT CASE(DerivLn)
 CASE("Yes")
  CALL Derivxy5pln2d(nux,nuz,VUx,VUz,Fe,Dfdux,Dfduz)
  Dfdux(1,:)= 0.
 CASE("No ")
  CALL Derivxy5p2d(nux,nuz,VUx,VUz,Fe,Dfdux,Dfduz)
  Dfdux(1,:)= 0.
END SELECT
CALL Derivxy5p2d(nux,nuz,VUx,VUz,Dfdux,Dfduxdux,Dfduxduz)
CALL Derivxy5p2d(nux,nuz,VUx,VUz,Dfduz,Dfduzdux,Dfduzduz)
Dfduxdux(1,:)= 0.
Dfduxduz(1,:)= 0.
Dfduzdux(1,:)= 0.

! Operator L-ux:
 DO k= 1,nuz
  FORALL(i=2:nux-1)
   Auxi(i)= (TANH(SQRT((VUx(i))**2+(VUz(k))**2)))**2
   Alpha1(i)= ( ( CAx(i,k)+DDxxdux(i,k) + DDzxduz(i,k)) * Auxx - CDxx(i,k)*Auxxx ) *Auxi(i)
   Beta1(i)= 1. + (- DAxdux(i,k)*DTau/2. + 2.*CDxx(i,k)*Auxxx ) *Auxi(i)
   Gamma1(i)= (- ( CAx(i,k)+DDxxdux(i,k)+DDzxduz(i,k) ) * Auxx - CDxx(i,k)*Auxxx ) *Auxi(i)
   Psi1(i)= Fe(i,k) + ( DAxdux(i,k)*Fe(i,k)*(DTau/2.) &
        + ( CDxz(i,k)*Dfduxduz(i,k) ) * (DTau) &
        + ( CAx(i,k)+DDxxdux(i,k)+DDzxduz(i,k) ) * Auxx * (Fe(i+1,k)-Fe(i-1,k)) &
        + CDxx(i,k)*Auxxx * (Fe(i+1,k)-2.*Fe(i,k)+Fe(i-1,k)) ) *Auxi(i)
   Faux1b(i)= Fe(i,k)
  END FORALL
  i= 1
  Psi1(i)= Psi1(i)-Alpha1(i)*Fe(i-1,k)
  Alpha1(i)= 0.
  Faux1b(i)= Fe(i,k)
  i= nux
  Auxi(i)= (TANH(SQRT((VUx(i))**2+(VUz(k))**2)))**2
  Alpha1(i)= ( CAx(i,k)+DDxxdux(i,k) +DDzxduz(i,k)) * (DTau/2./Dux) *Auxi(i)
  Beta1(i)= 1. + (- DAxdux(i,k)*DTau/2. &
       - ( CAx(i,k)+DDxxdux(i,k)+DDzxduz(i,k) ) * (DTau/2./Dux) ) *Auxi(i)
  Gamma1(i)= 0.
  Psi1(i)= Fe(i,k) + ( DAxdux(i,k)*Fe(i,k)*(DTau/2.) &
       + ( CDxz(i,k)*Dfduxduz(i,k) ) * (DTau) &
       + ( CAx(i,k)+DDxxdux(i,k)+DDzxduz(i,k) ) * (DTau/2./Dux) * (Fe(i,k)-Fe(i-1,k)) &
       + CDxx(i,k)*(DTau/Dux/Dux) * (Fe(i,k)-2.*Fe(i-1,k)+Fe(i-2,k)) ) *Auxi(i)
  Faux1b(i)= Fe(i,k)

  CALL Tridag(2,nux,Alpha1,Beta1,Gamma1,Psi1,Faux1,Baux1,Gaux1)

  ! Temporary boundary condition:
  Faux1(1)= Fe(1,k)
  Faux1b(1)= Fe(1,k)

  ! Correcting instabilities:
  CALL Cor_Ampli_2(Faux1b,Faux1,nux)

  ! Boundary condition:
  ! Faux1(1)= Faux1(2)

  FORALL(i=1:nux)
   Fex(i,k)= Faux1(i)
  END FORALL
 END DO

! Operator L-uz:
 DO i= 1,nux
  FORALL(k=2:nuz-1)
   Auxk(k)= (TANH(SQRT((VUx(i))**2+(VUz(k))**2)))**2
   Alpha2(k)= ( ( CAz(i,k)+DDzzduz(i,k)+DDxzdux(i,k) ) * Auxz - CDzz(i,k)*Auxzz ) *Auxk(k)
   Beta2(k)= 1. + (- DAzduz(i,k)*DTau/2. + 2.*CDzz(i,k)*Auxzz) *Auxk(k)
   Gamma2(k)= (-( CAz(i,k)+DDzzduz(i,k)+DDxzdux(i,k) ) * Auxz - CDzz(i,k)*Auxzz ) *Auxk(k)
   Psi2(k)= Fe(i,k) + ( DAzduz(i,k)*Fe(i,k)*(DTau/2.) &
        + ( CDzx(i,k)*Dfduzdux(i,k) ) * (DTau) &
        + ( CAz(i,k)+DDzzduz(i,k)+DDxzdux(i,k) ) * Auxz * (Fe(i,k+1)-Fe(i,k-1)) &
        + CDzz(i,k)*Auxzz * (Fe(i,k+1)-2.*Fe(i,k)+Fe(i,k-1)) )*Auxk(k)
   Faux2b(k)= Fe(i,k)
  END FORALL
  k= 1
  Auxk(k)= (TANH(SQRT((VUx(i))**2+(VUz(k))**2)))**2
  Alpha2(k)= 0.
  Beta2(k)= 1. + ( - DAzduz(i,k)*DTau/2. &
       + ( CAz(i,k)+DDzzduz(i,k)+DDxzdux(i,k) ) * (DTau/2./Duz) ) *Auxk(k)
  Gamma2(k)= - ( CAz(i,k)+DDzzduz(i,k)+DDxzdux(i,k) ) * (DTau/2./Duz)*Auxk(k)
  Psi2(k)= Fe(i,k) + ( DAzduz(i,k)*Fe(i,k)*(DTau/2.) &
       + ( CDzx(i,k)*Dfduzdux(i,k) ) * (DTau) &
       + ( CAz(i,k)+DDzzduz(i,k)+DDxzdux(i,k) ) * (DTau/2./Duz) * (Fe(i,k+1)-Fe(i,k)) &
       + CDzz(i,k)*(DTau/Duz/Duz) * (Fe(i,k+2)-2.*Fe(i,k+1)+Fe(i,k)) ) *Auxk(k)
  Faux2b(k)= Fe(i,k)
  k= nuz
  Auxk(k)= (TANH(SQRT((VUx(i))**2+(VUz(k))**2)))**2
  Alpha2(k)= ( CAz(i,k)+DDzzduz(i,k)+DDxzdux(i,k) ) * (DTau/2./Duz)*Auxk(k)
  Beta2(k)= 1. + (- DAzduz(i,k)*DTau/2. &
       - ( CAz(i,k)+DDzzduz(i,k)+DDxzdux(i,k) ) * (DTau/2./Duz) ) *Auxk(k)
  Gamma2(k)= 0.
  Psi2(k)= Fe(i,k) + ( DAzduz(i,k)*Fe(i,k)*(DTau/2.) &
       + ( CDzx(i,k)*Dfduzdux(i,k) ) * (DTau) &
       + ( CAz(i,k)+DDzzduz(i,k)+DDxzdux(i,k) ) * (DTau/2./Duz) * (Fe(i,k)-Fe(i,k-1)) &
       + CDzz(i,k)*(DTau/Duz/Duz) * (Fe(i,k)-2.*Fe(i,k-1)+Fe(i,k-2)) ) *Auxk(k)
  Faux2b(k)= Fe(i,k)

  CALL Tridag(1,nuz,Alpha2,Beta2,Gamma2,Psi2,Faux2,Baux2,Gaux2)

  ! Correcting instabilities:
  CALL Cor_Ampli_2(Faux2b,Faux2,nuz)

  FORALL(k=1:nuz)
   Fez(i,k)= Faux2(k)
  END FORALL
 END DO

Fe= Fex+Fez-Fe
! Boundary condition at Ux=0 (zero derivative):
DO k= 1,nuz
 IF (SQRT((VUx(1))**2+(VUz(k))**2)>Ucrit) THEN
  Fe(1,k)= Fe(2,k)
 ELSE
 END IF
END DO

! Limitador 2

! it=0
! Fe0
! FeOld = Fe0
! Calcula Fe1
! Do nothing
! DeltaFe = Fe1-FeOld

! it=1
! Fe1
! FeOld = Fe1
! Calcula Fe2
! Testa se Fe2-Fe1 > DeltaFe
! 	Se sim, Fe2 = Fe1+DeltaFe
! 	Se nao, Fe2 = Fe2 calculada
! DeltaFe = Fe2-FeOld

! IF(it > 1) THEN
! 	DO k=1, nuz
! 		DO i=1, nux
! 			IF ( ABS(Fe(i,k)-FeOld(i,k)) > 1.01*ABS(DeltaFe(i,k)) ) THEN
! 				IF (Fe(i,k)-FeOld(i,k)>0) THEN
! 					Fe(i,k) = FeOld(i,k) + 1.01*ABS(DeltaFe(i,k))
! 				ELSE
! 					Fe(i,k) = FeOld(i,k) - 1.01*ABS(DeltaFe(i,k))
! 				END IF
! 			END IF
! 		END DO
! 	END DO
! END IF

! DeltaFe = Fe - FeOld


RETURN
END SUBROUTINE Split

SUBROUTINE Evol_Iwave(DTau,Tau,WaveType)

! See Notes B-11.

USE Common_Params
USE Common_Arrays
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL*8, INTENT(in) :: DTau,Tau
INTEGER :: i,k,sigma
INTEGER :: Nvar,ind
REAL*8, DIMENSION(nux,nuz) :: Dfdux,Dfduz
REAL*8, DIMENSION(nqx,nqz) :: Iwave,Iwavenew
REAL*8, DIMENSION(nqx,nqz) :: CoefA,CoefB
REAL*8, DIMENSION(nqx*nqz) :: Ft,Dfdt,Ftout
CHARACTER(LEN=6), INTENT(in) :: WaveType

SELECT CASE(DerivLn)
 CASE("Yes")
  CALL Derivxy5pln2d(nux,nuz,VUx,VUz,Fe,Dfdux,Dfduz)
 CASE("No ")
  CALL Derivxy5p2d(nux,nuz,VUx,VUz,Fe,Dfdux,Dfduz)
END SELECT
Dfdux(1,:)= 0.

SELECT CASE(WaveType)
 CASE("Lwavep")
  sigma= 1
  Iwave= ILp
  CALL Coef_Lwave(sigma,Dfdux,Dfduz,CoefA,CoefB)
 CASE("Lwavem")
  sigma= -1
  Iwave= ILm
  CALL Coef_Lwave(sigma,Dfdux,Dfduz,CoefA,CoefB)
 CASE("Swavep")
  sigma= 1
  Iwave= ISp
  CALL Coef_Swave(sigma,Dfdux,Dfduz,CoefA,CoefB)
 CASE("Swavem")
  sigma= -1
  Iwave= ISm
  CALL Coef_Swave(sigma,Dfdux,Dfduz,CoefA,CoefB)
 CASE("Twavep")
  sigma= 1
  Iwave= ITp
  CALL Coef_Twave(sigma,CoefA,CoefB)
 CASE("Twavem")
  sigma= -1
  Iwave= ITm
  CALL Coef_Twave(sigma,CoefA,CoefB)
END SELECT
CoefARK= CoefA
CoefBRK= CoefB

!H= DTau
Nvar= nqx*nqz
ind= 0
DO i= 1,nqx
 DO k= 1,nqz
  ind= ind+1
  Ft(ind)= Iwave(i,k)
 END DO
END DO
!IF( (Tau+H-Tau2)*(Tau+H-Tau1) > 0. ) H=Tau2-Tau
CALL DERIVS(Tau+DTau,Ft,Dfdt,Nvar)
CALL RK4(Ft,Dfdt,Nvar,Tau,DTau,Ftout)
CALL Cor_Ampli_2(Ft,Ftout,Nvar)
Ft= Ftout
ind= 0
DO i= 1,nqx
 DO k= 1,nqz
  ind= ind+1
  IF (Ft(ind)<Infinity) THEN
   Iwavenew(i,k)= Ft(ind)
  ELSE 
   Iwavenew(i,k)= Infinity
  END IF
 END DO
END DO

SELECT CASE(WaveType)
 CASE("Lwavep")
  ILp= Iwavenew
 CASE("Lwavem")
  ILm= Iwavenew
 CASE("Swavep")
  ISp= Iwavenew
 CASE("Swavem")
  ISm= Iwavenew
 CASE("Twavep")
  ITp= Iwavenew
 CASE("Twavem")
  ITm= Iwavenew
END SELECT

RETURN
END SUBROUTINE Evol_Iwave

!CHECKPOINT
SUBROUTINE Correct_ILp(it)
USE Common_Params
USE Common_Arrays
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
INTEGER, INTENT(in) :: it
INTEGER :: i, k
REAL*8 :: increment

increment = 1.1

IF(it > 1) THEN
	DO k=1, nuz
		DO i=1, nux
			IF ( ABS(ILp(i,k)-ILpOld(i,k)) > increment*ABS(DeltaILp(i,k)) ) THEN
				IF (ILp(i,k)-ILpOld(i,k)>0) THEN
					ILp(i,k) = ILpOld(i,k) + increment*ABS(DeltaILp(i,k))
				ELSE
					ILp(i,k) = ILpOld(i,k) - increment*ABS(DeltaILp(i,k))
				END IF
			END IF
		END DO
	END DO
END IF

DeltaILp = ILp - ILpOld
ILpOld = ILp

END SUBROUTINE Correct_ILp

SUBROUTINE Correct_ILm(it)
USE Common_Params
USE Common_Arrays
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
INTEGER, INTENT(in) :: it
INTEGER :: i, k
REAL*8 :: increment

increment = 1.1

IF(it > 1) THEN
	DO k=1, nuz
		DO i=1, nux
			IF ( ABS(ILm(i,k)-ILmOld(i,k)) > increment*ABS(DeltaILm(i,k)) ) THEN
				IF (ILm(i,k)-ILmOld(i,k)>0) THEN
					ILm(i,k) = ILmOld(i,k) + increment*ABS(DeltaILm(i,k))
				ELSE
					ILm(i,k) = ILmOld(i,k) - increment*ABS(DeltaILm(i,k))
				END IF
			END IF
		END DO
	END DO
END IF

DeltaILm = ILm - ILmOld
ILmOld = ILm

END SUBROUTINE Correct_ILm



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
!* Mod. Aug/10: Fortran 95, FORALL instead of DO.
!
SUBROUTINE Tridag(If,k,A,B,C,D,V,Beta,Gamma)
! 	PARAMETER(K=100)
IMPLICIT NONE
INTEGER, INTENT(IN) :: If,k 
INTEGER :: Lim,Ifp1,Last,i,j
REAL*8, INTENT(in) :: A(k),B(k),C(k),D(k)
REAL*8, INTENT(out) :: V(k)
REAL*8, INTENT(out) :: Beta(k),Gamma(k)
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

!... Compute final solution vector V .....
V(Lim)=Gamma(Lim)
Last=Lim-If
DO j=1,Last
 i=Lim-j
 V(i)=Gamma(i)-C(i)*V(i+1)/Beta(i)
END DO
RETURN
END SUBROUTINE Tridag

SUBROUTINE DERIVS(Tau,Ft,Dfdt,N)
USE Common_Params
USE Common_Arrays
IMPLICIT NONE
!INTEGER, PARAMETER :: NMAX=6561    ! NMAX= 81*81
INTEGER, INTENT(in) :: N
INTEGER :: i,k,ind
REAL*8 :: Tau
REAL*8, DIMENSION(N) :: Ft,DFdt
ind= 0
DO i= 1,nqx
 DO k= 1,nqz
  ind= ind+1
  Dfdt(ind)= CoefARK(i,k)+CoefBRK(i,k)*Ft(ind)
 END DO
END DO
Tau= 0.    ! Just to avoid compilation warning
RETURN
END SUBROUTINE DERIVS

SUBROUTINE RK4(Y,DYDX,N,X,H,YOUT)
! Given values for N variables Y and their derivatives DYDX
! known at X, use the fourth-order Runge-Kutta method to advan-
! ce the solution over an interval H and return the incremen-
! ted variables as YOUT, which need not be a distinct array from
! Y. The user supplies the subroutine DERIVS(X,Y,DYDX,N) which
! returns derivatives DYDX at X.
! See NUMERICAL RECIPES.
! Version Fortran 95, Feb. 2007.

IMPLICIT NONE
INTEGER, PARAMETER :: NMAX=6561    ! NMAX= 81*81
INTEGER, INTENT(in) :: N
INTEGER :: I
REAL*8, DIMENSION(N) :: Y,DYDX,YOUT
REAL*8, DIMENSION(NMAX) :: YT,DYT,DYM
REAL*8 :: X,H,HH,H6,XH
        
HH=H*0.5
H6=H/6.
XH=X+HH
! CALL  DERIVS(X,Y,DYDX)
DO I=1,N
 YT(I)=Y(I)+HH*DYDX(I)
END DO
CALL DERIVS(XH,YT,DYT,N)
DO I=1,N
 YT(I)=Y(I)+HH*DYT(I)
END DO
CALL DERIVS(XH,YT,DYM,N)
DO I=1,N
 YT(I)=Y(I)+H*DYM(I)
 DYM(I)=DYT(I)+DYM(I)
END DO
CALL DERIVS(X+H,YT,DYT,N)
DO I=1,N
 YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.*DYM(I))
END DO
RETURN
END SUBROUTINE RK4

SUBROUTINE Cor_Ampli(Vf,n)
USE Math_Constants
IMPLICIT NONE
INTEGER :: n,i
REAL*8, DIMENSION (n) :: Vf(n)

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

SUBROUTINE Cor_Ampli_2(Vf,Vfnew,n)
USE Math_Constants
IMPLICIT NONE
INTEGER :: n,i
REAL*8, DIMENSION (n) :: Vf(n),Vfnew(n)

! Corrects negative amplitudes when they are not acceptable.
! For a given real*8vector "Vfnew", replaces negative values by the previous 
! value ("Vf").

DO i= 1,n
! IF ( Vfnew(i) < EpsMin ) THEN
!  Vfnew(i)= Vf(i) 
! ELSE
! END IF
 IF ( Vfnew(i) < 0. ) THEN
!  Vfnew(i)= 0. 
  Vfnew(i)= Vf(i) 
 ELSE
 END IF
END DO

RETURN
END SUBROUTINE Cor_Ampli_2

SUBROUTINE Aitp1d2(nx,Vx,Fx,Xp,Fp,i)
! Uses linear interpolation in order to obtain the value of a function
! F, at the point Xp. The function F is given as a set of points Fx(x).
! Uses subroutine Locate (Numerical Recipes, P. 96)
! Version Fortran 95, Aug, 2006.
! Version of Aitp1d, without the call to Locate, which is called outside.

IMPLICIT NONE
INTEGER :: nx,i
REAL*8, DIMENSION(nx) :: Vx,Fx
REAL*8 :: Xp,Fp,Aux

!CALL Locate(Vx,nx,Xp,i)

IF ( i<=0 .OR. i>=nx ) THEN
 FP= 0.
ELSE
 Aux= ( Xp-Vx(i) )/( Vx(i+1)-Vx(i) )
 Fp= Fx(i) + ( Fx(i+1)-Fx(i) )*Aux
END IF 
RETURN
END SUBROUTINE Aitp1d2

SUBROUTINE Locate(Xx,N,X,J)
! Given an array Xx of lenght N, and given a value X, returns a value 
! J such that X is between Xx(J) and Xx(J+1).
! Xx must be monotonic, either increasing or decreasing.
! J=0 or J=N is returned to indicate that X is out of range.
! See NUMERICAL RECIPES.
! Version Fortran 95, Aug 2006.

IMPLICIT NONE
INTEGER :: N,J,JL,JU,JM
REAL*8, DIMENSION (N) :: Xx
REAL*8 :: X

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
 
SUBROUTINE Aitp2d(nx,ny,Vx,Vy,Fxy,Xp,Yp,Fp)
! Uses linear interpolation in order to obtain the value of a function
! F, at the point (Xp,Yp). 
! The function F is given as a set of points Fxy(x,y).
! Uses subroutine Locate (Numerical Recipes, P. 96)
! Version Fortran 95, Sep., 2006.

IMPLICIT NONE
REAL*8, DIMENSION(nx) :: Vx
REAL*8, DIMENSION(ny) :: Vy
REAL*8, DIMENSION(nx,ny) :: Fxy
REAL*8 :: Xp,Yp,Fp
REAL*8 :: F1,F2,F3,F4,T,U
INTEGER :: nx,ny,i,j

CALL Locate(Vx,nx,Xp,i)
CALL Locate(Vy,ny,Yp,j)

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

SUBROUTINE Aitp2d2(nx,ny,Vx,Vy,Fxy,Xp,Yp,Fp,i,j)
! Uses linear interpolation in order to obtain the value of a function
! F, at the point (Xp,Yp). 
! The function F is given as a set of points Fxy(x,y).
! Uses subroutine Locate (Numerical Recipes, P. 96)
! Version Fortran 95, Sep., 2006.
! Version of Aitp2d, without the calls to Locate, which is called outside.

IMPLICIT NONE
REAL*8, DIMENSION(nx) :: Vx
REAL*8, DIMENSION(ny) :: Vy
REAL*8, DIMENSION(nx,ny) :: Fxy
REAL*8 :: Xp,Yp,Fp
REAL*8 :: F1,F2,F3,F4,T,U
INTEGER :: nx,ny,i,j

!CALL Locate(Vx,nx,Xp,i)
!CALL Locate(Vy,ny,Yp,j)

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
END SUBROUTINE Aitp2d2

SUBROUTINE Aitp2d2b(nx,ny,Vx,Vy,Fxy,Xp,Yp,Fp,i,j)
! Uses linear interpolation in order to obtain the value of a function
! F, at the point (Xp,Yp). 
! The function F is given as a set of points Fxy(x,y).
! Uses subroutine Locate (Numerical Recipes, P. 96)
! Version Fortran 95, Sep., 2006.
! Version of Aitp2d, without the calls to Locate, which is called outside.
! Corrects the interpolation at the extreme points, making the value at the
! points outside of the grid equal to the nearest point inside the grid.
! It is only for interpolation using similar grids, to avoid errors at the edge.

IMPLICIT NONE
REAL*8, DIMENSION(nx) :: Vx
REAL*8, DIMENSION(ny) :: Vy
REAL*8, DIMENSION(nx,ny) :: Fxy
REAL*8 :: Xp,Yp,Fp
REAL*8 :: F1,F2,F3,F4,T,U
INTEGER :: nx,ny,i,j

!CALL Locate(Vx,nx,Xp,i)
!CALL Locate(Vy,ny,Yp,j)

IF (i==0) THEN
 IF (j==0) THEN
  FP= Fxy(1,1)
 ELSE
  IF (j==ny) THEN
   FP= Fxy(1,ny)
  ELSE
   FP= Fxy(1,j)+ (Yp-Vy(j))/(Vy(j+1)-Vy(j))*(Fxy(1,j+1)-Fxy(1,j)) 
  END IF
 END IF
ELSE
 IF (i==nx) THEN
  IF (j==0) THEN
   FP= Fxy(nx,1)
  ELSE
   IF (j==ny) THEN
    FP= Fxy(nx,ny)
   ELSE
    FP= Fxy(nx,j)+ (Yp-Vy(j))/(Vy(j+1)-Vy(j))*(Fxy(nx,j+1)-Fxy(nx,j)) 
   END IF
  END IF
 ELSE
  IF (j==0) THEN
   FP= Fxy(i,1)+ (Xp-Vx(i))/(Vx(i+1)-Vx(i))*(Fxy(i+1,1)-Fxy(i,1)) 
  ELSE
   IF (j==ny) THEN
    FP= Fxy(i,ny)+ (Xp-Vx(i))/(Vx(i+1)-Vx(i))*(Fxy(i+1,ny)-Fxy(i,ny)) 
   ELSE
    F1= Fxy(i,j)
    F2= Fxy(i+1,j)
    F3= Fxy(i+1,j+1)
    F4= Fxy(i,j+1)
    T= ( Xp-Vx(i) )/( Vx(i+1)-Vx(i) )
    U= ( Yp-Vy(j) )/( Vy(j+1)-Vy(j) )
    FP= (1.-T)*(1.-U)*F1 + T*(1.-U)*F2 + T*U*F3 + (1.-T)*U*F4
   END IF
  END IF
 END IF
END IF
RETURN
END SUBROUTINE Aitp2d2b

SUBROUTINE Simpson(Vx,F,N,Res)
! Version Fortran 95, Aug 2006.
IMPLICIT NONE
INTEGER :: N,Nm1,Nm2,i
REAL*8, DIMENSION(N) :: F,Vx
REAL*8 :: Res,H

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

      SUBROUTINE DQsimp(func,a,b,s)
! Version 22/Aug/2015
IMPLICIT NONE
REAL*8 :: a,b,s
REAL*8 :: EPS= 1.d-4
INTEGER :: j
INTEGER :: jmax= 11
REAL*8 :: os,ost,st
!      EXTERNAL func
INTERFACE 
 REAL*8 FUNCTION FUNC(X)
  REAL*8, INTENT(in) :: X
 END FUNCTION FUNC
END INTERFACE
      ost=-1.d30
      os= -1.d30
      do j=1,JMAX
        call dtrapzd(func,a,b,st,j)
        s=(4.d0*st-ost)/3.d0
        if (j.gt.5) then
          if (dabs(s-os).lt.EPS*dabs(os).or. &
              !(s.eq.0.d0.and.os.eq.0.d0)) return
              (dabs(s).lt.1.d-16.and.dabs(os).lt.1.d-16)) return
        endif
        os=s
        ost=st
      enddo
      OPEN(98,file="Warning_qsimp")
      write(98,*) " Too many steps in qsimp! "
      close(98)
!      STOP
      END subroutine DQsimp


      SUBROUTINE DTrapzd(func,a,b,s,n)
! Version 22/Aug/2015
IMPLICIT NONE
INTEGER :: n
INTEGER :: it,j
REAL*8 :: a,b,s
REAL*8 :: del,suma,tnm,x
!      EXTERNAL func
INTERFACE 
 REAL*8 FUNCTION FUNC(X)
  REAL*8, INTENT(in) :: X
 END FUNCTION FUNC
END INTERFACE
      if (n.eq.1) then
        s=0.5d0*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5d0*del
        suma=0.d0
        do j=1,it
          suma=suma+func(x)
          x=x+del
        enddo
        s=0.5d0*(s+(b-a)*suma/tnm)
      endif
      return
      END subroutine DTrapzd

      SUBROUTINE DQsimpb(func,a,b,s)
! Version 22/Aug/2015
! To be used if "a' is finite, positive, and "b" goes to infinity, or if
! "b" is finite, negative, and "a" goes to -infinity.
IMPLICIT NONE
REAL*8 :: a,b,s
REAL*8 :: EPS= 1.d-4
INTEGER :: j
INTEGER :: jmax= 11
REAL*8 :: os,ost,st
!      EXTERNAL func
INTERFACE 
 REAL*8 FUNCTION FUNC(X)
  REAL*8, INTENT(in) :: X
 END FUNCTION FUNC
END INTERFACE
      ost=-1.d30
      os= -1.d30
      do j=1,JMAX
        call dmidinf(func,a,b,st,j)
        s=(4.d0*st-ost)/3.d0
        if (j.gt.5) then
          if (dabs(s-os).lt.EPS*dabs(os).or. &
              !(s.eq.0.d0.and.os.eq.0.d0)) return
              (dabs(s).lt.1.d-16.and.dabs(os).lt.1.d-16)) return
        endif
        os=s
        ost=st
      enddo
      OPEN(98,file="Warning_qsimpb")
      write(98,*) " Too many steps in qsimpb! "
      close(98)
!      STOP
      END subroutine DQsimpb

      SUBROUTINE DMIDINF(FUNK,AA,BB,S,N) 
!	Numerical Recipes, Cambridge, 1989, p. 118:
!	This routine is an exact replacement for MIDPNT, i.e. return as S 
!	the Nth stage of refinement of the integral of FUNK from AA to BB,
!	except that the function is evaluated at evenly spaced points in 1/x
!	rather than in x. This allows the upper limit BB to be as large and
!	positive as the computer allows, or the lower limit AA to be as large
!	and negative, but not both. AA and BB must have the same sign.
! Version 22/Aug/2015
IMPLICIT NONE
INTEGER :: n,j
INTEGER, SAVE :: it
REAL*8 :: aa,bb,s,a,b
REAL*8 :: del,sum,tnm,x,ddel,func
!      EXTERNAL func
INTERFACE 
 REAL*8 FUNCTION FUNK(X)
  REAL*8, INTENT(in) :: X
 END FUNCTION FUNK
END INTERFACE
      FUNC(X)=FUNK(1.D0/X)/X**2   ! This is a statement function which effects
				! the change of variable. 
      B=1.D0/AA                 ! These two statements change the limits of
				! integration accordingly.
      A=1.D0/BB 
      IF (N.EQ.1) THEN          ! From this point on, the routine is exactly
				! identical to MIDPNT.
        S=(B-A)*FUNC(0.5D0*(A+B)) 
        IT=1 
      ELSE 
        TNM=IT 
        DEL=(B-A)/(3.D0*TNM) 
        DDEL=DEL+DEL 
        X=A+0.5D0*DEL 
        SUM=0.D0 
        DO 11 J=1,IT 
          SUM=SUM+FUNC(X) 
          X=X+DDEL 
          SUM=SUM+FUNC(X) 
          X=X+DEL 
11      CONTINUE 
        S=(S+(B-A)*SUM/TNM)/3.D0 
        IT=3*IT 
      ENDIF 
      RETURN 
END subroutine DMIDINF 

SUBROUTINE Derivxy5p2d(nx,ny,Vx,Vy,Fxy,Dfdx,Dfdy)
! Version Jan. 2010: 
! Evaluates the x-derivative and the y-derivative of two-dimensional
! function Fxy(x,y), given by an array on (nx,ny) elements

! Version Fortran 95, Sept. 2006.
! Uses 5-point derivative for the internal points.
! Version Fortran 95, Aug. 2010 (uses FORALL instead of DO).
! See Abramowitz & Stegun, 1970, Eq. 25.3.6, for 'p=0'
! Based on "deriv5p.f", College Park, 16/Feb/2001 (L. F. Ziebell)
! Needs equally spaced points!

IMPLICIT NONE
INTEGER :: i,j,nx,ny
INTEGER :: nxm1,nym1,nxm2,nym2
REAL*8, DIMENSION(nx) :: Vx
REAL*8, DIMENSION(ny) :: Vy
REAL*8, DIMENSION(nx,ny) :: Fxy(nx,ny),Dfdx(nx,ny),Dfdy(nx,ny)
REAL*8 :: Dx,Dy

Dx= Vx(2)-Vx(1)
Dy= Vy(2)-Vy(1)

nxm1= nx-1
nxm2= nx-2
FORALL(j=1:ny)
  Dfdx(1,j)= (Fxy(2,j)-Fxy(1,j))/Dx
  Dfdx(2,j)= (Fxy(3,j)-Fxy(1,j))/(2.*Dx)
  Dfdx(nx,j)= (Fxy(nx,j)-Fxy(nxm1,j))/Dx
  Dfdx(nxm1,j)= (Fxy(nx,j)-Fxy(nxm2,j))/(2.*Dx)
  FORALL(i=3:nxm2)
   Dfdx(i,j)= ( (Fxy(i+1,j)-Fxy(i-1,j))*2./3. &
        - (Fxy(i+2,j)-Fxy(i-2,j))/12. ) / Dx
  END FORALL
END FORALL

nym1= ny-1
nym2= ny-2
FORALL(i=1:nx)
 Dfdy(i,1)= (Fxy(i,2)-Fxy(i,1))/Dy
 Dfdy(i,2)= (Fxy(i,3)-Fxy(i,1))/(2.*Dy)
 Dfdy(i,ny)= (Fxy(i,ny)-Fxy(i,nym1))/Dy
 Dfdy(i,nym1)= (Fxy(i,ny)-Fxy(i,nym2))/(2.*Dy)
 FORALL(j=3:nym2)
  Dfdy(i,j)= ( (Fxy(i,j+1)-Fxy(i,j-1))*2./3. &
       - (Fxy(i,j+2)-Fxy(i,j-2))/12. ) / Dy
 END FORALL
END FORALL

RETURN
END SUBROUTINE Derivxy5p2d
!
!
SUBROUTINE Derivxy5pln2d(nx,ny,Vx,Vy,Fxy,Dfdx,Dfdy)
! Version Jan. 2010: Introduces logarithm approach to improve derivatives
! of functions with small values.

! Version Fortran 95, Sept. 2006.
! Uses 5-point derivative for the internal points.
! See Abramowitz & Stegun, 1970, Eq. 25.3.6, for 'p=0'
! Based on "deriv5p.f", College Park, 16/Feb/2001 (L. F. Ziebell)
! Needs equally spaced points!

! Evaluates the x-derivative and the y derivative of two-dimensional
! function Fxy(x,y), given by an array on (nx,ny) elements

IMPLICIT NONE
INTEGER :: i,j,nx,ny
INTEGER :: nxm1,nym1,nxm2,nym2
REAL*8, DIMENSION(nx) :: Vx
REAL*8, DIMENSION(ny) :: Vy
REAL*8, DIMENSION(nx,ny) :: Fxy(nx,ny),Dfdx(nx,ny),Dfdy(nx,ny)
REAL*8, DIMENSION(nx,ny) :: Fold(nx,ny),Faux(nx,ny)
INTEGER, DIMENSION(2) :: Imin
REAL*8 :: Dx,Dy
REAL*8 :: Fxymin

Fold= Fxy
Imin= MINLOC(Fxy)
Fxymin= Fxy(Imin(1),Imin(2))
Fxy= Fxy+ABS(Fxymin)
Fxy= Fxy+1.0E-30
Faux= Fxy
Fxy= LOG(Fxy)
 
Dx= Vx(2)-Vx(1)
Dy= Vy(2)-Vy(1)

nxm1= nx-1
nxm2= nx-2
FORALL(j=1:ny)
  Dfdx(1,j)= (Fxy(2,j)-Fxy(1,j))/Dx * Faux(1,j)
  Dfdx(2,j)= (Fxy(3,j)-Fxy(1,j))/(2.*Dx) * Faux(2,j)
  Dfdx(nx,j)= (Fxy(nx,j)-Fxy(nxm1,j))/Dx * Faux(nx,j)
  Dfdx(nxm1,j)= (Fxy(nx,j)-Fxy(nxm2,j))/(2.*Dx) *Faux(nxm1,j)
  FORALL(i=3:nxm2)
   Dfdx(i,j)= ( (Fxy(i+1,j)-Fxy(i-1,j))*2./3. &
        - (Fxy(i+2,j)-Fxy(i-2,j))/12. ) / Dx * Faux(i,j)
  END FORALL
END FORALL

nym1= ny-1
nym2= ny-2
FORALL(i=1:nx)
  Dfdy(i,1)= (Fxy(i,2)-Fxy(i,1))/Dy * Faux(i,1)
  Dfdy(i,2)= (Fxy(i,3)-Fxy(i,1))/(2.*Dy) * Faux(i,2)
  Dfdy(i,ny)= (Fxy(i,ny)-Fxy(i,nym1))/Dy * Faux(i,ny)
  Dfdy(i,nym1)= (Fxy(i,ny)-Fxy(i,nym2))/(2.*Dy) * Faux(i,nym1)
  FORALL(j=3:nym2)
   Dfdy(i,j)= ( (Fxy(i,j+1)-Fxy(i,j-1))*2./3. &
        - (Fxy(i,j+2)-Fxy(i,j-2))/12. ) / Dy * Faux(i,j)
  END FORALL
END FORALL

Fxy= Fold
RETURN
END SUBROUTINE Derivxy5pln2d
!

SUBROUTINE Derivx5p2d(nx,ny,Vx,Fxy,Dfdx)
! Version Jan. 2010: 
! Evaluates the x-derivative of two-dimensional
! function Fxy(x,y), given by an array on (nx,ny) elements

! Based on "Derivxy5p2d", version Fortran 95, Sept. 2006.
! Uses 5-point derivative for the internal points.
! Version Fortran 95, Aug. 2010 (uses FORALL instead of DO).
! See Abramowitz & Stegun, 1970, Eq. 25.3.6, for 'p=0'
! Based on "deriv5p.f", College Park, 16/Feb/2001 (L. F. Ziebell)
! Needs equally spaced points!

IMPLICIT NONE
INTEGER :: i,j,nx,ny
INTEGER :: nxm1,nxm2
REAL*8, DIMENSION(nx) :: Vx
REAL*8, DIMENSION(nx,ny) :: Fxy(nx,ny),Dfdx(nx,ny)
REAL*8 :: Dx

Dx= Vx(2)-Vx(1)

nxm1= nx-1
nxm2= nx-2
FORALL(j=1:ny)
  Dfdx(1,j)= (Fxy(2,j)-Fxy(1,j))/Dx
  Dfdx(2,j)= (Fxy(3,j)-Fxy(1,j))/(2.*Dx)
  Dfdx(nx,j)= (Fxy(nx,j)-Fxy(nxm1,j))/Dx
  Dfdx(nxm1,j)= (Fxy(nx,j)-Fxy(nxm2,j))/(2.*Dx)
  FORALL(i=3:nxm2)
   Dfdx(i,j)= ( (Fxy(i+1,j)-Fxy(i-1,j))*2./3. &
        - (Fxy(i+2,j)-Fxy(i-2,j))/12. ) / Dx
  END FORALL
END FORALL
RETURN
END SUBROUTINE Derivx5p2d

!
SUBROUTINE Derivy5p2d(nx,ny,Vy,Fxy,Dfdy)
! Version Jan. 2010: 
! Evaluates the y-derivative of two-dimensional
! function Fxy(x,y), given by an array on (nx,ny) elements

! Based on "Derivxy5p2d", version Fortran 95, Sept. 2006.
! Uses 5-point derivative for the internal points.
! Version Fortran 95, Aug. 2010 (uses FORALL instead of DO).
! See Abramowitz & Stegun, 1970, Eq. 25.3.6, for 'p=0'
! Based on "deriv5p.f", College Park, 16/Feb/2001 (L. F. Ziebell)
! Needs equally spaced points!

IMPLICIT NONE
INTEGER :: i,j,nx,ny
INTEGER :: nym1,nym2
REAL*8, DIMENSION(ny) :: Vy
REAL*8, DIMENSION(nx,ny) :: Fxy(nx,ny),Dfdy(nx,ny)
REAL*8 :: Dy

Dy= Vy(2)-Vy(1)

nym1= ny-1
nym2= ny-2
FORALL(i=1:nx)
 Dfdy(i,1)= (Fxy(i,2)-Fxy(i,1))/Dy
 Dfdy(i,2)= (Fxy(i,3)-Fxy(i,1))/(2.*Dy)
 Dfdy(i,ny)= (Fxy(i,ny)-Fxy(i,nym1))/Dy
 Dfdy(i,nym1)= (Fxy(i,ny)-Fxy(i,nym2))/(2.*Dy)
 FORALL(j=3:nym2)
  Dfdy(i,j)= ( (Fxy(i,j+1)-Fxy(i,j-1))*2./3. &
       - (Fxy(i,j+2)-Fxy(i,j-2))/12. ) / Dy
 END FORALL
END FORALL

RETURN
END SUBROUTINE Derivy5p2d
!
!
REAL*8 FUNCTION PERF(X)
! Error function of argument X, based on PERF
! (from R. L. Meyer, Universite de Nancy, Franca)
! Version Feb. 18, 2009.

IMPLICIT NONE
!REAL*8 :: PERFC
REAL*8 :: X,X2,SNUM,SDEN
REAL*8, DIMENSION(5) :: P,Q
DATA P /3.209377589138469472562E+03,3.774852376853020208137E+02,&
1.138641541510501556495E+02,3.161123743870565596947E+00,&
1.857777061846031526730E-01/
DATA Q /2.844236833439170622273E+03,1.282616526077372275645E+03,&
2.440246379344441733056E+02,2.360129095234412093499E+01,1./
IF (X <= 0.5) THEN
 X2=X*X
 SNUM=X*(P(1)+X2*(P(2)+X2*(P(3)+X2*(P(4)+X2*P(5)))))
 SDEN=Q(1)+X2*(Q(2)+X2*(Q(3)+X2*(Q(4)+X2*Q(5))))
 PERF=SNUM/SDEN
ELSE
 PERF=1.-PERFC(X)
END IF
RETURN
END FUNCTION PERF

REAL*8 FUNCTION PERFC(X)
! Complementary error function of argument X, based on PERF 
! (from R. L. Meyer, Universite de Nancy, Franca)
! Version Feb. 18, 2009.

IMPLICIT NONE
!REAL*8 :: PERF
REAL*8 :: X,X2,Y2,SNUM,SDEN
REAL*8, DIMENSION(9) :: P,Q
REAL*8, DIMENSION(6) :: R,S
DATA P /1.23033935479799725272E+03,2.05107837782607146532E+03,&
1.71204761263407058314E+03,8.81952221241769090411E+02,&
2.98635138197400131132E+02,6.61191906371416294775E+01,&
8.88314979438837594118E+00,5.64188496988670089180E-01,&
2.15311535474403846343E-08/
DATA Q /1.23033935480374942043E+03,3.43936767414372163696E+03,&
4.36261909014324715820E+03,3.29079923573345962678E+03,&
1.62138957456669018874E+03,5.37181101862009857509E+02,&
1.17693950891312499305E+02,1.57449261107098347253E+01,1./
DATA R /-6.58749161529837803157E-04,-1.60837851487422766278E-02,&
-1.25781726111229246204E-01,-3.60344899949804439429E-01,&
-3.05326634961232344035E-01,-1.63153871373020978498E-02/
DATA S /2.33520497626869185443E-03,6.05183413124413191178E-02,&
5.27905102951428412248E-01,1.87295284992346047209E+00,&
2.56852019228982242072E+00,1./
IF (X > 90.) THEN
 PERFC=0.
ELSE
 IF (X <= 0.5) THEN
  PERFC=1.-PERF(X)
 ELSE
  X2=X*X
  IF (X <= 4.) THEN
   SNUM=EXP(-X2)*(P(1)+X*(P(2)+X*(P(3)+X*(P(4)+X*(P(5)+X*(P(6)&
        +X*(P(7)+X*(P(8)+X*P(9)))))))))
   SDEN=Q(1)+X*(Q(2)+X*(Q(3)+X*(Q(4)+X*(Q(5)+X*(Q(6)+X*(Q(7)&
        +X*(Q(8)+X*Q(9))))))))
   PERFC=SNUM/SDEN
  ELSE
   Y2=1./X2
   SNUM= R(1)+Y2*(R(2)+Y2*(R(3)+Y2*(R(4)+Y2*(R(5)+Y2*R(6)))))
   SDEN=S(1)+Y2*(S(2)+Y2*(S(3)+Y2*(S(4)+Y2*(S(5)+Y2*S(6)))))
   PERFC=EXP(-X2)*(0.56418958354775628694+SNUM/SDEN*Y2)/X
  END IF
 END IF
END IF
RETURN
END FUNCTION PERFC

REAL*8 FUNCTION PERFCE(X)
! Evaluates EXP(X**2)*ERFC(X
! Based on DERFCE(X)
! (from R. L. Meyer, Universite de Nancy, Franca)
! Version Oct. 09, 2015

IMPLICIT NONE
REAL*8 :: X,X2,SNUM,SDEN,Y2
REAL*8, DIMENSION(9) :: P,Q
REAL*8, DIMENSION(6) :: R,S
DATA P/1.23033935479799725272E+03,2.05107837782607146532E+03,&
 1.71204761263407058314E+03,8.81952221241769090411E+02,&
 2.98635138197400131132E+02,6.61191906371416294775E+01,&
 8.88314979438837594118E+00,5.64188496988670089180E-01,&
 2.15311535474403846343E-08/
DATA Q/1.23033935480374942043E+03,3.43936767414372163696E+03,&
 4.36261909014324715820E+03,3.29079923573345962678E+03,&
 1.62138957456669018874E+03,5.37181101862009857509E+02,&
 1.17693950891312499305E+02,1.57449261107098347253E+01,&
 1./
DATA R/-6.58749161529837803157E-04,-1.60837851487422766278E-02,&
 -1.25781726111229246204E-01,-3.60344899949804439429E-01,&
 -3.05326634961232344035E-01,-1.63153871373020978498E-02/
DATA S/2.33520497626869185443E-03,6.05183413124413191178E-02,&
 5.27905102951428412248E-01,1.87295284992346047209E+00,&
 2.56852019228982242072E+00,1./
IF (X <= 0.5) THEN
 PERFCE=(1.-PERF(X))*EXP(X*X)
ELSE
 X2=X*X
 IF (X <=4.) THEN
  SNUM= (P(1)+X*(P(2)+X*(P(3)+X*(P(4)+X*(P(5)+X*(P(6)+X*(P( &
     7)+X*(P(8)+X*P(9)))))))))
  SDEN= Q(1)+X*(Q(2)+X*(Q(3)+X*(Q(4)+X*(Q(5)+X*(Q(6)+X*(Q(7)+X*(Q(8)+ &
     X*Q(9))))))))
  PERFCE=SNUM/SDEN
 ELSE
  Y2=1./X2
  SNUM= R(1)+Y2*(R(2)+Y2*(R(3)+Y2*(R(4)+Y2*(R(5)+Y2*R(6)))))
  SDEN=S(1)+Y2*(S(2)+Y2*(S(3)+Y2*(S(4)+Y2*(S(5)+Y2*S(6)))))
  PERFCE= (0.56418958354775628694+SNUM/SDEN*Y2)/X
 END IF
END IF
RETURN
END FUNCTION PERFCE

SUBROUTINE Coef_Coll
USE Common_Params
USE Common_Arrays
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL*8 :: Ux,Uz,U,U2,U3,U5
REAL*8 :: Cee,Cei,RMeMe,RVeVi,RVeVe,RViVe
REAL*8 :: PsiArg,PhiArg,PhipArg,Arg,Arg2
REAL*8 :: ColAxe,ColAze,ColAxi,ColAzi
REAL*8 :: ColDxxe,ColDxze,ColDzze,ColDxxi,ColDxzi,ColDzzi
!REAL*8 :: PERF
REAL*8 :: Aux,ZZ
INTEGER :: l,m

ZZ= 1  ! Version which works only for ions with unit charge.
Cee= (2.*Pi)*Geff * LOG(3./(8*Pi)/SQRT(2.)/Geff)
Cei= Cee*ZZ**2
RMeMe= 1.
RVeVi= SQRT(RTeTi*RMiMe)
RVeVe= 1.
RViVe= 1./RVeVi

SELECT CASE(CollTerm)

 CASE("Yes")
  SELECT CASE(CollTermForm)

   CASE("Complete")
     DO m= 1,nuz
      Uz= VUz(m)
      DO l= 1,nux
       Ux= VUx(l)
       U2= Ux**2+Uz**2
       U= SQRT(U2)
       U3= U2*U
       U5= U3*U2
       U3= U3+EpsMin
       U5= U5+EpsMin
       Aux= (TANH(U))**2

       Arg= U*RVeVe
       Arg2= Arg**2
       PhiArg= PERF(Arg)
       PhipArg= 2.*EXP(-Arg2)/SQRT(Pi)
       PsiArg= PhiArg-Arg*PhipArg
       ColAxe= Cee*(2./RMeMe*PsiArg*Ux/U3)
       ColAze= Cee*(2./RMeMe*PsiArg*Uz/U3)
       ColDxxe= Cee*((U2*PhiArg-RVeVe**2*PsiArg/2.)*Uz*Uz/U5 &
        + RVeVe**2*PsiArg*Ux*Ux/U5)
       ColDxze= - Cee*((U2*PhiArg-RVeVe**2*PsiArg/2.)*Ux*Uz/U5 &
        - RVeVe**2*PsiArg*Ux*Uz/U5)
       ColDzze= Cee*((U2*PhiArg-RVeVe**2*PsiArg/2.)*Ux*Ux/U5 &
        + RVeVe**2*PsiArg*Uz*Uz/U5)

       Arg= U*RVeVi
       Arg2= Arg**2
       PhiArg= PERF(Arg)
       PhipArg= 2.*EXP(-Arg2)/SQRT(Pi)
       PsiArg= PhiArg-Arg*PhipArg
       ColAxi= Cei*(2./RMiMe*PsiArg*Ux/U3)
       ColAzi= Cei*(2./RMiMe*PsiArg*Uz/U3)
       ColDxxi= Cei*((U2*PhiArg-RViVe**2*PsiArg/2.)*Uz*Uz/U5 &
        + RViVe**2*PsiArg*Ux*Ux/U5)
       ColDxzi= - Cei*((U2*PhiArg-RViVe**2*PsiArg/2.)*Ux*Uz/U5 &
        - RViVe**2*PsiArg*Ux*Uz/U5)
       ColDzzi= Cei*((U2*PhiArg-RViVe**2*PsiArg/2.)*Ux*Ux/U5 &
        + RViVe**2*PsiArg*Uz*Uz/U5)

       ColAx(l,m)= Aux*(ColAxe+ColAxi)
       ColAz(l,m)= Aux*(ColAze+ColAzi)
       ColDxx(l,m)= Aux*(ColDxxe+ColDxxi)
       ColDxz(l,m)= Aux*(ColDxze+ColDxzi)
       ColDzz(l,m)= Aux*(ColDzze+ColDzzi)
      END DO
     END DO

   CASE("Expanded")
     DO m= 1,nuz
      Uz= VUz(m)
      DO l= 1,nux
       Ux= VUx(l)
       U2= Ux**2+Uz**2
       U= SQRT(U2)
       U3= U2*U
       U5= U3*U2
       U3= U3+EpsMin
       U5= U5+EpsMin
       Aux= (TANH(U))**2
      
       Arg= U*RVeVe
       Arg2= Arg**2
       PhiArg= 1.
       PsiArg= 1.
       ColAxe= Cee*(2./RMeMe*PsiArg*Ux/U3)
       ColAze= Cee*(2./RMeMe*PsiArg*Uz/U3)
       ColDxxe= Cee*((U2*PhiArg-RVeVe**2*PsiArg/2.)*Uz*Uz/U5 &
        + RVeVe**2*PsiArg*Ux*Ux/U5)
       ColDxze= - Cee*((U2*PhiArg-RVeVe**2*PsiArg/2.)*Ux*Uz/U5 &
        - RVeVe**2*PsiArg*Ux*Uz/U5)
       ColDzze= Cee*((U2*PhiArg-RVeVe**2*PsiArg/2.)*Ux*Ux/U5 &
        + RVeVe**2*PsiArg*Uz*Uz/U5)

       Arg= U*RVeVi
       Arg2= Arg**2
       PhiArg= 1.
       PsiArg= 1.
       ColAxi= Cei*(2./RMiMe*PsiArg*Ux/U3)
       ColAzi= Cei*(2./RMiMe*PsiArg*Uz/U3)
       ColDxxi= Cei*((U2*PhiArg-RViVe**2*PsiArg/2.)*Uz*Uz/U5 &
        + RViVe**2*PsiArg*Ux*Ux/U5)
       ColDxzi= - Cei*((U2*PhiArg-RViVe**2*PsiArg/2.)*Ux*Uz/U5 &
        - RViVe**2*PsiArg*Ux*Uz/U5)
       ColDzzi= Cei*((U2*PhiArg-RViVe**2*PsiArg/2.)*Ux*Ux/U5 &
        + RViVe**2*PsiArg*Uz*Uz/U5)

       ColAx(l,m)= Aux*(ColAxe+ColAxi)
       ColAz(l,m)= Aux*(ColAze+ColAzi)
       ColDxx(l,m)= Aux*(ColDxxe+ColDxxi)
       ColDxz(l,m)= Aux*(ColDxze+ColDxzi)
       ColDzz(l,m)= Aux*(ColDzze+ColDzzi)
      END DO
     END DO

  CASE DEFAULT
   OPEN(98,FILE='Warning_Coef_Coll.wt')
   WRITE(98,*) ' CollTerm= ',CollTerm
   WRITE(98,*) ' CollTerm must be (Yes) or (No ) !!'
   CLOSE(98)
   STOP

  END SELECT

 CASE("No ")
   DO m= 1,nuz
    DO l= 1,nux
     ColAx(l,m)= 0.
     ColAz(l,m)= 0.
     ColDxx(l,m)= 0.
     ColDxz(l,m)= 0.
     ColDzz(l,m)= 0.
    END DO
   END DO

 CASE DEFAULT
  OPEN(98,FILE='Warning_Coef_Coll.wt')
  WRITE(98,*) ' CollTerm= ',CollTerm
  WRITE(98,*) ' CollTerm must be (Yes) or (No ) !!'
  CLOSE(98)
  STOP

END SELECT

RETURN
END SUBROUTINE Coef_Coll

SUBROUTINE Coll_Damping
USE Common_Params
USE Common_Arrays
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL*8 :: Qx,Qz,Q2,Q
REAL*8 :: Zlq,Zsq,Zsoq
REAL*8 :: RTemp0,Res
REAL*8 :: Res1L,Res2L,Res1S,Res2S
INTEGER :: i,j,m,sigma

REAL*8 :: Aux
REAL*8 :: RTempf
INTEGER :: ires

REAL*8, DIMENSION(nqx,nqz) :: GqlL,GcollL
REAL*8, DIMENSION(nqx,nqz) :: GqlS,GcollS

m= 1   ! Auxiliary to the space profiles (just one point, for the moment)
RTemp0= VRTeTs(m)
RTempf= SQRT(VRTfparTs(m)**2 + VRTfperpTs(m)**2)

! Evaluation of the 1D expression:
GqlL1D= 0.
GcollL1D= 0.
GqlS1D= 0.
GcollS1D= 0.
DO i= 1,nqmod
 ! For the L waves:
 Q= VQQ(i)
 IF (Lsaturation=="Yes" .AND. Q<1./SQRT(1./Ve2C2-1.5)) THEN
  Q= 1./SQRT(1./Ve2C2-1.5)
 ELSE
  Q= Q
 END IF
 Q2= Q**2
 Zlq= SQRT(VRNeNs(m))*SQRT(1.+1.5*Q2*VRTeTs(m)/VRNeNs(m))
! Contribution of the background distribution (aproximated, U0=0)
 Aux1_Gcoll(1)= Rtemp0
 Aux1_Gcoll(2)= Q
 Aux1_Gcoll(3)= Q2
 sigma= 1   ! Evaluated only for sigma=1, due to the symmetry
 Aux2_Gcoll(1)= sigma
 Aux1_Gcoll(4)= Zlq
 CALL DQsimp(Aux_GcollL,1.D-4,0.28D0,Res1L)
 CALL DQsimpb(Aux_GcollL,0.28D0,1.D30,Res2L)
 Res= Res1L+Res2L
 GcollL1D(i)= -Zlq**2*(16.*SQRT(Pi))/(SQRT(VRTeTs(m)))**3/Q2 &
  * (VRNeNs(m)/VRTeTs(m))**2*Geff * Res
 HGcollL1D(i)= -(8.)*Ue*Q/(SQRT(VRTeTs(m)))**3 &
  * (VRNeNs(m)/VRTeTs(m))*Geff * Res
! Evaluates the quasilinear damping:
  GqlL1D(i)= -SQRT(Pi)/(SQRT(VRTeTs(m)))**3*Zlq**2/Q**3 &
    * EXP(-Zlq**2/Q2/VRTeTs(m))
 ! For the S waves:
 Q= VQQ(i)
 Q2= Q**2
 Zsq= Q*AA*SQRT(VRTeTs(m))/SQRT(1.+Q2/2.*VRTeTs(m)/VRNeNs(m))
 Zsoq= AA*SQRT(VRTeTs(m))/SQRT(1.+Q2/2.*VRTeTs(m)/VRNeNs(m))
! Contribution of the background distribution (aproximated, U0=0)
 Aux1_Gcoll(1)= Rtemp0
 Aux1_Gcoll(2)= Q
 Aux1_Gcoll(3)= Q2
 sigma= 1   ! Evaluated only for sigma=1, due to the symmetry
 Aux2_Gcoll(1)= sigma
 Aux1_Gcoll(4)= Zsoq
 CALL DQsimp(Aux_GcollS,1.D-4,1.D0,Res1S)
 CALL DQsimpb(Aux_GcollS,1.D0,1.D30,Res2S)
 Res= Res1S+Res2S
 GcollS1D(i)= -Zlq*Zsq*(16.*SQRT(Pi)) * (Q*AA/2.) &
  * (VRNeNs(m)/VRTeTs(m))**(5./2.)*Geff * Res
 HGcollS1D(i)= -(8.)*Ue*Q *(VRNeNs(m))**2/(VRTeTs(m))**3*Geff &
  * Res
! Evaluates the quasilinear damping:
 GqlS1D(i)= -SQRT(Pi)*(AA/2.)*Zlq*Zsq &
   * (EXP(-Zsoq**2/VRTeTs(m))/(SQRT(VRTeTs(m)))**3 &
   + EXP(-Zsoq**2/VRTeTs(m)*RTeTi*RMiMe)*SQRT(RMiMe*RTeTi/VRTeTs(m)) &
   * RTeTi/VRTeTs(m) )

! Contribution from the beam distribution:
! Not included!
END DO     ! Q

! Evaluation of the 2D expression:
GqlL= 0.
GcollL= 0.
GqlS= 0.
GcollS= 0.
DO i= 1,nqx
 Qx= VQx(i)
 DO j= 1,nqz
  Qz= VQz(j)
  Q2= Qx**2+Qz**2
  Q= SQRT(Q2)
  !IF(Q<=5.E-3) Q=5.E-3
  CALL Locate(VQQ,nqmod,Q,ires)
  CALL Aitp1d2(nqmod,VQQ,GcollL1D,Q,Aux,ires)
  GcollLp(i,j)= Aux
  GcollLm(i,j)= Aux
  CALL Aitp1d2(nqmod,VQQ,GqlL1D,Q,Aux,ires)
  GqlLp(i,j)= Aux
  GqlLm(i,j)= Aux
  CALL Aitp1d2(nqmod,VQQ,GcollS1D,Q,Aux,ires)
  GcollSp(i,j)= Aux
  GcollSm(i,j)= Aux
  CALL Aitp1d2(nqmod,VQQ,GqlS1D,Q,Aux,ires)
  GqlSp(i,j)= Aux
  GqlSm(i,j)= Aux
  CALL Aitp1d2(nqmod,VQQ,HGcollL1D,Q,Aux,ires)
  HGcollLp(i,j)= Aux
  CALL Aitp1d2(nqmod,VQQ,HGcollS1D,Q,Aux,ires)
  HGcollSp(i,j)= Aux
 END DO    ! Qz
END DO     ! Qx

RETURN
END SUBROUTINE Coll_Damping

REAL*8 FUNCTION Aux_GcollL(Qp)
USE Common_Params
USE Common_Arrays
USE Math_Constants
IMPLICIT NONE
REAL*8, INTENT(in) :: Qp
REAL*8 :: Rtemp,Q,Q2,Zlq,Qp2EpsQpSq,Qp2
REAL*8 :: AuxA,AuxB,AuxA2,AuxB2,Aux
INTEGER :: sigma,m
!INTEGER :: i,j
COMPLEX*16 :: Qp2EpsQp
COMPLEX*16 :: Zeta

m= 1
Rtemp= Aux1_Gcoll(1)
Q= Aux1_Gcoll(2)
Q2= Aux1_Gcoll(3)
Zlq= Aux1_Gcoll(4) 
sigma= Aux2_Gcoll(1) 

IF(Q > 2.D-4) THEN
 IF(Qp > 2.D-4) THEN
  Qp2= Qp**2
  Zeta= CMPLX(sigma*Zlq/Qp/SQRT(VRTeTs(m)),0.E0,8)
  Qp2EpsQp= Qp2+2.D0*VRNeNs(m)/VRTeTs(m) *(1.+Zeta*Zfn(Zeta))
  !EpsQp2= ((Zlq**2-VRNeNs(m)-1.5*VRTeTs(m)*Qp2)**2+4.D0*Pi*(VRNeNs(m))**5 &
  !  /(VRTeTs(m))**3/Qp**6*EXP(-2.D0*Zlq**2/Qp2))/(Zlq)**4
  Qp2EpsQpSq= (CDABS(Qp2EpsQp))**2

  AuxA= -2.D0*Q*Qp
  AuxB= 2.D0*(1.D0+RTeTi)*VRNeNs(m)/VRTeTs(m)+Q2+Qp2
  AuxA2= AuxA**2
  AuxB2= AuxB**2
  Aux= (Qp*(-AuxA2+2.d0*AuxB2)/(AuxB2-AuxA2)+(AuxB/2.D0/Q) &
    * DLOG(1.D0+2.D0*AuxA/(AuxB-AuxA)))
  Aux_GcollL= DEXP(-Zlq**2/Qp2/VRTeTs(m))*Aux/Qp2EpsQpSq
 ELSE
  Aux_GcollL= 0.D0
 END IF
ELSE
 Aux_GcollL= 0.D0
END IF
RETURN
END FUNCTION Aux_GcollL

REAL*8 FUNCTION Aux_GcollS(Qp)
USE Common_Params
USE Common_Arrays
USE Math_Constants
IMPLICIT NONE
REAL*8, INTENT(in) :: Qp
REAL*8 :: Rtemp,Q,Q2,Zsoq,EpsQp2,Qp2
REAL*8 :: AuxA,AuxB,AuxA2,AuxB2,Aux
INTEGER :: sigma,m
!INTEGER :: i,j
COMPLEX*16 :: EpsQp
COMPLEX*16 :: Zetae,Zetai

m= 1
Rtemp= Aux1_Gcoll(1)
Q= Aux1_Gcoll(2)
Q2= Aux1_Gcoll(3)
Zsoq= Aux1_Gcoll(4) 
sigma= Aux2_Gcoll(1) 

IF(Q > 2.D-4) THEN
 IF(Qp > 2.D-4) THEN
  Qp2= Qp**2
  Zetae= sigma*Zsoq/SQRT(VRTeTs(m))
  Zetai= sigma*Zsoq/SQRT(VRTeTs(m))*SQRT(RTeTi*RMiMe)
  EpsQp= 1.D0+2.D0*VRNeNs(m)/VRTeTs(m)/Qp2 *(1.D0+Zetae*Zfn(Zetae)) &
     + 2.D0*VRNeNs(m)/VRTeTs(m)*RTeTi/Qp2 *(1.D0+Zetai*Zfn(Zetai))
  EpsQp2= (CDABS(EpsQp))**2

  AuxA= -2.D0*Q*Qp
  AuxB= 2.D0*(1.D0+RTeTi)*VRNeNs(m)/VRTeTs(m)+Q2+Qp2
  AuxA2= AuxA**2
  AuxB2= AuxB**2
  Aux= (4.D0/(AuxB2-AuxA2)+RTeTi*Qp/Q/Q2/Qp2*(-2.D0*AuxA*AuxB/(AuxB2-AuxA2) &
    + DLOG((AuxA+AuxB)/(AuxB-AuxA))))
  Aux_GcollS= (DEXP(-Zsoq**2/VRTeTs(m))/VRTeTs(m)/SQRT(VRTeTs(m)) &
   + EXP(-Zsoq**2/VRTeTs(m)*RTeTi*RMiMe)*RTeTi/VRTeTs(m) &
   * SQRT(RTeTi/VRTeTs(m)*RMiMe) ) * Aux/Qp**3/EpsQp2
 ELSE
  Aux_GcollS= 0.D0
 END IF
ELSE
 Aux_GcollS= 0.D0
END IF
RETURN
END FUNCTION Aux_GcollS

SUBROUTINE Bremsstrahlung
USE Common_Params
USE Common_Arrays
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL*8 :: Qx,Qz,Q2,Q,Muq
REAL*8 :: Zlq,Zsq
REAL*8 :: RTemp0,Res
REAL*8 :: Res1L,Res2L,Res1S,Res2S
INTEGER :: i,j,m

REAL*8 :: Aux
REAL*8 :: RTempf
INTEGER :: ires

INTEGER, PARAMETER :: nmu= 51
REAL*8, PARAMETER :: Mupi= -1.E0+1.E-6, Mupf= 1.E0-1.E-6
REAL*8 :: Mup,Dmup
REAL*8, DIMENSION(nmu) :: Vmup, VintL, VintS

m= 1   ! Auxiliary to the space profiles (just one point, for the moment)
RTemp0= VRTeTs(m)
RTempf= SQRT(VRTfparTs(m)**2 + VRTfperpTs(m)**2)

DMup= (Mupf-Mupi)/(nmu-1)
DO i= 1,nmu
 VMup(i)= Mupi+(i-1)*DMup
END DO
VMup(nmu)= Mupf

! Evaluation of the 1D expression:
BremL1D= 0.
BremS1D= 0.
HBremL1D= 0.
HBremS1D= 0.
DO i= 1,nqmod
 ! For the L waves:
 Q= VQQ(i)
 IF (Lsaturation=="Yes" .AND. Q<1./SQRT(1./Ve2C2-1.5)) THEN
  Q= 1./SQRT(1./Ve2C2-1.5)
 ELSE
  Q= Q
 END IF
 Q2= Q**2
 Zlq= SQRT(VRNeNs(m))*SQRT(1.+1.5*Q2*VRTeTs(m)/VRNeNs(m))
! Contribution of the background distribution (aproximated, U0=0)
 Aux1_Brem(1)= Rtemp0
 Aux1_Brem(2)= Q
 Aux1_Brem(3)= Q2
 DO j= 1,nmu
  Mup= Vmup(j)
  Aux1_Brem(5)= Mup
  Aux1_Brem(4)= Zlq
  CALL DQsimp(Aux_BremL,1.D-4,1.D0,Res1L)
  CALL DQsimpb(Aux_BremL,1.D0,1.D30,Res2L)
  VintL(j)= Res1L+Res2L
 END DO
 CALL Simpson(VMup,VintL,nmu,Res)
 BremL1D(i)= 384*SQRT(Pi)/(VRTeTs(m))**2*(VRNeNs(m))**4*Geff**2 &
  * (1.-RTeTi**2/RMiMe)**2/Q2 * Res
 HBremL1D(i)= 384*Ue*Q/(VRTeTs(m))**2*(VRNeNs(m))**2*Geff &
  * (1.-RTeTi**2/RMiMe)**2 * Res
 ! For the S waves:
 Q= VQQ(i)
 Q2= Q**2
 Zsq= Q*AA*SQRT(VRTeTs(m))/SQRT(1.+Q2/2.*VRTeTs(m)/VRNeNs(m))
 Muq= Q**3*AA/2.
! Contribution of the background distribution (aproximated, U0=0)
 Aux1_Brem(1)= Rtemp0
 Aux1_Brem(2)= Q
 Aux1_Brem(3)= Q2
 DO j= 1,nmu
  Mup= Vmup(j)
  Aux1_Brem(5)= Mup
  Aux1_Brem(4)= Zsq
  CALL DQsimp(Aux_Brem,1.D-4,1.D0,Res1S)
  CALL DQsimpb(Aux_Brem,1.D0,1.D30,Res2S)
  VintS(j)= Res1S+Res2S
 END DO
 CALL Simpson(VMup,VintS,nmu,Res)
 BremS1D(i)= 96*SQRT(Pi)/(VRTeTs(m))**4*(VRNeNs(m))**5*Geff**2 &
  !* (1.-RTeTi**2)**2*Muq/Q2 * Res
  * (1.-RTeTi**2)**2*AA*Q/2. * Res
 HBremS1D(i)= 96*Ue*Q/(VRTeTs(m))**5*(VRNeNs(m))**3*Geff &
  * (1.-RTeTi**2)**2*SQRT(VRNeNs(m)/VRTeTs(m)) * Res
END DO

! Contribution from the beam distribution:
! Not included!

! Evaluation of the 2D expression:
BremL= 0.
BremS= 0.
DO i= 1,nqx
 Qx= VQx(i)
 DO j= 1,nqz
  Qz= VQz(j)
  Q2= Qx**2+Qz**2
  Q= SQRT(Q2)
  !IF(Q<=5.E-3) Q=5.E-3
  CALL Locate(VQQ,nqmod,Q,ires)
  CALL Aitp1d2(nqmod,VQQ,BremL1D,Q,Aux,ires)
  BremL(i,j)= Aux
  CALL Aitp1d2(nqmod,VQQ,BremS1D,Q,Aux,ires)
  BremS(i,j)= Aux
  CALL Aitp1d2(nqmod,VQQ,HBremL1D,Q,Aux,ires)
  HBremL(i,j)= Aux
  CALL Aitp1d2(nqmod,VQQ,HBremS1D,Q,Aux,ires)
  HBremS(i,j)= Aux
 END DO    ! Qz
END DO     ! Qx

RETURN
END SUBROUTINE Bremsstrahlung

REAL*8 FUNCTION Aux_Brem(Qp)
! Expression utilized for S waves
USE Common_Params
USE Common_Arrays
USE Math_Constants
IMPLICIT NONE
REAL*8, INTENT(in) :: Qp
REAL*8 :: Rtemp,Q,Q2,Zq,Qp2,Mup
REAL*8 :: Aux,Aux1,Aux2,Aux3
INTEGER :: m

m= 1
Rtemp= Aux1_Brem(1)
Q= Aux1_Brem(2)
Q2= Aux1_Brem(3)
Zq= Aux1_Brem(4) 
Mup= Aux1_Brem(5)

Qp2= Qp**2
Aux1= 2.*VRNeNs(m)/VRTeTs(m)*(1.+RTeTi)+Qp2
Aux2= Aux1+Q2-2.*Q*Qp*Mup

Aux= 0.
Aux3= (Qp2+(Q2+Qp2-2.*Q*Qp*Mup))/VRTeTs(m)
Aux= Aux+DEXP(-Zq**2/Aux3)/DSQRT(Aux3)
Aux3= (1./RTeTi/RMiMe*Qp2+(Q2+Qp2-2.*Q*Qp*Mup))/VRTeTs(m)
Aux= Aux+DEXP(-Zq**2/Aux3)/DSQRT(Aux3)
Aux3= (1./RTeTi/RMiMe*Qp2+1./RTeTi/RMiMe*(Q2+Qp2-2.*Q*Qp*Mup))/VRTeTs(m)
Aux= Aux+DEXP(-Zq**2/Aux3)/DSQRT(Aux3)
Aux3= (Qp2+1./RTeTi/RMiMe*(Q2+Qp2-2.*Q*Qp*Mup))/VRTeTs(m)
Aux= Aux+DEXP(-Zq**2/Aux3)/DSQRT(Aux3)

Aux_Brem= Qp2/Aux1**2/Aux2**2*Aux
RETURN
END FUNCTION Aux_Brem

REAL*8 FUNCTION Aux_BremL(Qp)
! Expression utilized for L waves
USE Common_Params
USE Common_Arrays
USE Math_Constants
IMPLICIT NONE
REAL*8, INTENT(in) :: Qp
REAL*8 :: Rtemp,Q,Q2,Zq,Qp2,Mup
REAL*8 :: Aux,Aux1,Aux2,Aux3
INTEGER :: m

m= 1
Rtemp= Aux1_Brem(1)
Q= Aux1_Brem(2)
Q2= Aux1_Brem(3)
Zq= Aux1_Brem(4) 
Mup= Aux1_Brem(5)

Qp2= Qp**2
Aux1= 2.*VRNeNs(m)/VRTeTs(m)*(1.+RTeTi)+Qp2
Aux2= Aux1+Q2-2.*Q*Qp*Mup

Aux= 0.
Aux3= (Qp2+(Q2+Qp2-2.*Q*Qp*Mup))/VRTeTs(m)
Aux= Aux+DEXP(-Zq**2/Aux3)/DSQRT(Aux3)
Aux3= (1./RTeTi/RMiMe*Qp2+(Q2+Qp2-2.*Q*Qp*Mup))/VRTeTs(m)
Aux= Aux+DEXP(-Zq**2/Aux3)/DSQRT(Aux3)
Aux3= (1./RTeTi/RMiMe*Qp2+1./RTeTi/RMiMe*(Q2+Qp2-2.*Q*Qp*Mup))/VRTeTs(m)
Aux= Aux+DEXP(-Zq**2/Aux3)/DSQRT(Aux3)
Aux3= (Qp2+1./RTeTi/RMiMe*(Q2+Qp2-2.*Q*Qp*Mup))/VRTeTs(m)
Aux= Aux+DEXP(-Zq**2/Aux3)/DSQRT(Aux3)

Aux_BremL= Qp2**2*(Q2+Qp2-2.*Q*Qp*Mup)/Aux1**2/Aux2**2*Aux
RETURN
END FUNCTION Aux_BremL

SUBROUTINE Rebuild_L
! Based on the subroutine developed by J. Pavan, April 2009.

USE Common_Params
USE Common_Arrays
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL*8 :: P0,P1,P2,P3,t
REAL*8, DIMENSION(nqx,nqz) :: ILp1,ILm1
INTEGER :: i,k,j,isqr,ksqr

CALL Locate(VQx,nqx,QxSqr,isqr)
CALL Locate(VQz,nqz,QzSqr,ksqr)
!
! Uses Four-Points Bezier Method.
!==============================
! ILp,ILm: z=0 axis.
!
DO i=1,isqr-1   ! Until the beginning of the 'squared region'; 
                ! see 'SUBROUTINE Coef_Lwave'.
P0=ILm(i,4)
P1=ILm(i,3)
P2=ILp(i,3)
P3=ILp(i,4)
k=4
t=0.0
DO j=1,3
k=k-1
t=t+1./6.
ILm(i,k)=(1.-t)**3*P0+3.*t*(1.-t)**2*P1+3.*t**2*(1.-t)*P2+t**3*P3
END DO
k=4
t=1.0
DO j=1,3
k=k-1
t=t-1./6.
ILp(i,k)=(1.-t)**3*P0+3.*t*(1.-t)**2*P1+3.*t**2*(1.-t)*P2+t**3*P3
END DO
END DO
!=========================
! ILp,ILm: x=0 axis. 
!
DO k=1,nqz
P0=ILm(4,k)
P1=ILm(3,k)
P2=ILm(3,k)
P3=ILm(4,k)
i=4
t=1.0
DO j=1,3
i=i-1
t=t-1./6.
ILm(i,k)=(1.-t)**3*P0+3.*t*(1.-t)**2*P1+3.*t**2*(1.-t)*P2+t**3*P3
END DO
P0=ILp(4,k)
P1=ILp(3,k)
P2=ILp(3,k)
P3=ILp(4,k)
i=4
t=1.0
DO j=1,3
i=i-1
t=t-1./6.
ILp(i,k)=(1.-t)**3*P0+3.*t*(1.-t)**2*P1+3.*t**2*(1.-t)*P2+t**3*P3
END DO
END DO
!==============================
! L spectrum 'squared region'.
! Arbitrary weights.
!
ILp1=ILp
DO i=isqr,nqx
DO k=ksqr+1,1,-1
ILp1(i,k)=(0.5*ILp1(i-1,k)+1.5*ILp1(i,k+1))/2.
END DO
END DO
!
ILm1=ILm
DO i=isqr,nqx
DO k=ksqr+1,1,-1
ILm1(i,k)=(0.5*ILm1(i-1,k)+1.5*ILm1(i,k+1))/2.
END DO
END DO
!
ILp=ILp1
ILm=ILm1
!==============================
! ILp,ILm: z=0 axis, 'squared region'.
!
DO i=isqr,nqx
P0=ILm(i,4)
P1=ILm(i,3)
P2=ILp(i,3)
P3=ILp(i,4)
k=4
t=0.0
DO j=1,3
k=k-1
t=t+1./6.
ILm(i,k)=(1.-t)**3*P0+3.*t*(1.-t)**2*P1+3.*t**2*(1.-t)*P2+t**3*P3
END DO
k=4
t=1.0
DO j=1,3
k=k-1
t=t-1./6.
ILp(i,k)=(1.-t)**3*P0+3.*t*(1.-t)**2*P1+3.*t**2*(1.-t)*P2+t**3*P3
END DO
END DO
!
RETURN
END SUBROUTINE Rebuild_L

SUBROUTINE Rebuild_S
! Based on the subroutine developed by J. Pavan, April 2009.

USE Common_Params
USE Common_Arrays
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL*8 :: P0,P1,P2,P3,t
INTEGER :: i,k,j,isqr,ksqr

CALL Locate(VQx,nqx,QxSqr,isqr)
CALL Locate(VQz,nqz,QzSqr,ksqr)
!
! Uses Four-Points Bezier Method.
!==============================
! ISp,ISm: z=0 axis.
!
DO i=1,nqx
P0=ISm(i,4)
P1=ISm(i,3)
P2=ISp(i,3)
P3=ISp(i,4)
k=4
t=0.0
DO j=1,3
k=k-1
t=t+1./6.
ISm(i,k)=(1.-t)**3*P0+3.*t*(1.-t)**2*P1+3.*t**2*(1.-t)*P2+t**3*P3
END DO
k=4
t=1.0
DO j=1,3
k=k-1
t=t-1./6.
ISp(i,k)=(1.-t)**3*P0+3.*t*(1.-t)**2*P1+3.*t**2*(1.-t)*P2+t**3*P3
END DO
END DO
!=========================
! ISp,ISm: x=0 axis. 
!
DO k=1,nqz
P0=ISm(4,k)
P1=ISm(3,k)
P2=ISm(3,k)
P3=ISm(4,k)
i=4
t=1.0
DO j=1,3
i=i-1
t=t-1./6.
ISm(i,k)=(1.-t)**3*P0+3.*t*(1.-t)**2*P1+3.*t**2*(1.-t)*P2+t**3*P3
END DO
P0=ISp(4,k)
P1=ISp(3,k)
P2=ISp(3,k)
P3=ISp(4,k)
i=4
t=1.0
DO j=1,3
i=i-1
t=t-1./6.
ISp(i,k)=(1.-t)**3*P0+3.*t*(1.-t)**2*P1+3.*t**2*(1.-t)*P2+t**3*P3
END DO
END DO
!
RETURN
END SUBROUTINE Rebuild_S

SUBROUTINE Rebuild_Fe
! Based on the subroutine developed by J. Pavan, April 2009.

USE Common_Params
USE Common_Arrays
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL*8 :: P0,P1,P2,P3,t
INTEGER :: i,k,j,nuz2,isqr,ksqr

CALL Locate(VQx,nqx,QxSqr,isqr)
CALL Locate(VQz,nqz,QzSqr,ksqr)
!
! Uses Four-Points Bezier Method.
!==============================
! Fe: z=0 axis.
!
! nuz2= (nuz-1)/2+1
! DO i=1,nux
! 	P0=Fe(i,nuz2-3)
! 	P1=Fe(i,nuz2-2)
! 	P2=Fe(i,nuz2+2)
! 	P3=Fe(i,nuz2+3)
! 	k=nuz2-3
! 	t=0.0
! 	DO j=1,5
! 		k=k+1
! 		t=t+1./6.
! 		Fe(i,k)=(1.-t)**3*P0+3.*t*(1.-t)**2*P1+3.*t**2*(1.-t)*P2+t**3*P3
! 	END DO
! END DO
!
! Fe: x=0 axis.
!
DO k=1,nuz
	P0=Fe(4,k)
	P1=Fe(3,k)
	P2=Fe(3,k)
	P3=Fe(4,k)
	i=4
	t=0.0
	IF(VUz(k) .lt. 7.5 .and. VUz(k) .gt. 5.5) THEN
		DO j=1,3
			i=i-1
			t=t+1./6.
			Fe(i,k)=(1.-t)**3*P0+3.*t*(1.-t)**2*P1+3.*t**2*(1.-t)*P2+t**3*P3
		END DO
	END IF
END DO
!
RETURN
END SUBROUTINE Rebuild_Fe

SUBROUTINE Bezier(array)
! Based on the subroutine developed by J. Pavan, April 2009.

USE Common_Params
USE Common_Arrays
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL*8 :: P0,P1,P2,P3,t
REAL*8, DIMENSION(nux, nuz) :: array
INTEGER :: i,k,j,nuz2,isqr,ksqr

CALL Locate(VQx,nqx,QxSqr,isqr)
CALL Locate(VQz,nqz,QzSqr,ksqr)
!
! Uses Four-Points Bezier Method.
!==============================
! Fe: z=0 axis.
!
nuz2= (nuz-1)/2+1
DO i=1,nux
	P0=array(i,nuz2-3)
	P1=array(i,nuz2-2)
	P2=array(i,nuz2+2)
	P3=array(i,nuz2+3)
	k=nuz2-3
	t=0.0
	DO j=1,5
		k=k+1
		t=t+1./6.
		array(i,k)=(1.-t)**3*P0+3.*t*(1.-t)**2*P1+3.*t**2*(1.-t)*P2+t**3*P3
	END DO
END DO
!
! Fe: x=0 axis.
!
DO k=1,nuz
	P0=array(4,k)
	P1=array(3,k)
	P3=array(4,k)
	P2=array(3,k)
	i=4
	t=0.0
	DO j=1,3
		i=i-1
		t=t+1./6.
		array(i,k)=(1.-t)**3*P0+3.*t*(1.-t)**2*P1+3.*t**2*(1.-t)*P2+t**3*P3
	END DO
END DO
!
RETURN
END SUBROUTINE Bezier

SUBROUTINE BezierDDzzduz
! Based on the subroutine developed by J. Pavan, April 2009.

USE Common_Params
USE Common_Arrays
USE Math_Constants
USE Phys_Constants
IMPLICIT NONE
REAL*8 :: P0,P1,P2,P3,t
INTEGER :: i,k,j,nuz2,isqr,ksqr

CALL Locate(VQx,nqx,QxSqr,isqr)
CALL Locate(VQz,nqz,QzSqr,ksqr)
!
! Uses Four-Points Bezier Method.
!==============================
! Fe: z=0 axis.
!
! nuz2= (nuz-1)/2+1
! DO i=1,nux
! 	P0=DDzzduz(i,nuz2-3)
! 	P1=DDzzduz(i,nuz2-2)
! 	P2=DDzzduz(i,nuz2+2)
! 	P3=DDzzduz(i,nuz2+3)
! 	k=nuz2-3
! 	t=0.0
! 	DO j=1,5
! 		k=k+1
! 		t=t+1./6.
! 		DDzzduz(i,k)=(1.-t)**3*P0+3.*t*(1.-t)**2*P1+3.*t**2*(1.-t)*P2+t**3*P3
! 	END DO
! END DO
!
! Fe: x=0 axis.
!
DO k=1,nuz
	P0=DDzzduz(4,k)
	P1=DDzzduz(3,k)
	P2=DDzzduz(3,k)
	P3=DDzzduz(4,k)
	i=4
	t=0.0
	DO j=1,3
		i=i-1
		t=t+1./6.
		DDzzduz(i,k)=(1.-t)**3*P0+3.*t*(1.-t)**2*P1+3.*t**2*(1.-t)*P2+t**3*P3
	END DO
END DO
!
RETURN
END SUBROUTINE BezierDDzzduz

!MODULE Zfn_MODULE
!CONTAINS
!
!Função Z para plasma Maxwelliano inserido ao programa akfvenanisobikappapro1 em 29/08/2013

!c
!c
!c--------------------First Z-function code-----------------------------
!c
      complex*16 function zfn(z)
!c
!c     Evaluates the plasma dispersion function (Fried and Conte 
!c     function) of complex argument with a relative error of 1.e-6.
!c
!c     Algorithm: based closely on that described in Piero Barberio-
!c                Corsetti 'Calculation of the Plasma Dispersion 
!c                Function'. 
!c
!c     Precision: Double
!c
!c     Author: R. L. Mace, Plasma Physics Research Institute
!c               University of Natal, Durban
!c
!c     Modificacoes no programa original:
!c        1)Precisao default: tol= 1.0d-6, dlim= 4.0d+00.  Para aumentar
!c          precisao basta colocar tol= 1.0d-14, dlim= 6.0d+00.
!c
      complex*16 z
!c
!c     constants
!c
      !real*8 tol,zero,half,one,dlim,thrhlf,pid4
      real*8 tol,zero,half,one,dlim,pid4
      parameter( tol=1.0d-14, zero=0.d00, half=0.5d00, one=1.0d00 )
      !parameter( dlim=6.0d00, thrhlf=1.5d00, pid4=0.785398163397448d00 )
      parameter( dlim=6.0d00, pid4=0.785398163397448d00 )
      complex*16 czero,chalf,cone,ctwo,irtpi,i2rtpi
      parameter( czero=(0.d00,0.d00), chalf=(0.5d00,0.d00) )
      parameter( cone=(1.d00,0.d00), ctwo=(2.d00,0.d00) )
      parameter( irtpi=(0.d00,1.772453850905516d00) )
      parameter( i2rtpi=(0.d00,3.544907701811032d00) )
!c
!c     local variables
!c
      real*8 x,y,abx,aby,xymax,fn,cn,aslim,yasm
      complex*16 errz,an,anm1,bn,bnm1,anp1,bnp1
      complex*16 z2,zinv,aa,bb,sum,term,pterm
!c
      x=dreal(z)
      y=dimag(z)
      abx=dabs(x)
      aby=dabs(y)
      if (aby.gt.abx) then
        xymax=aby
      else
        xymax=abx
      endif
      fn=zero
!c
!c     based on the magnitude of the real and imaginary parts of z, 
!c     determine which of power series, continued fraction, or 
!c     asymptotic forms to use
!c
      if (aby.gt.one) then
!c
!c       **********************************
!c       employ the continued fraction form
!c       **********************************
!c
        z2=half-z*z
        an=z
        anm1=czero
        bn=z2
        bnm1=cone
        xymax=one/xymax
!c
!c       compute the continued fraction
!c
        zfn=an/bn
        errz=zfn-anm1/bnm1
        do while (dabs(dreal(errz)).gt.tol*dabs(dreal(zfn)) .or.&
        dabs(dimag(errz)).gt.tol*dabs(dimag(zfn)))
           fn=fn+one
           cn=xymax/fn
           aa=-fn*(fn-half)*cn
           bb=(z2+fn+fn)*cn
           anp1=bb*an+aa*anm1
           bnp1=bb*bn+aa*bnm1
           anm1=an*cn
           an=anp1
           bnm1=bn*cn
           bn=bnp1
           zfn=an/bn
           errz=zfn-anm1/bnm1
         end do
!c
!c        add the contribution from the pole if Im(z) .le. 0
!c
         if (y.le.zero) then
           zfn=zfn+i2rtpi*cdexp(-z*z)
         end if
      else if (abx.gt.dlim)  then
!c
!c        ****************************
!c        use the asmyptotic expansion
!c        ****************************
!c
         zinv=cone/z
         z2=chalf*zinv*zinv
         sum=cone
         term=cone
         aslim=x*x+y*y-one
         do while (&
           (dabs(dreal(term)).gt.tol*dabs(dreal(sum)) .or.&
           dabs(dimag(term)).gt.tol*dabs(dimag(sum))) .and.&
           fn.le.aslim &
           )
           fn=fn+one
           term=term*(fn+fn-one)*z2
           sum=sum+term
         end do
         zfn=-zinv*sum
         yasm=pid4/abx
         if (y.lt.-yasm) then
            zfn=zfn+i2rtpi*cdexp(-z*z)
         else if (y.le.yasm) then
            zfn=zfn+irtpi*cdexp(-z*z)
         end if
      else
!c
!c        *************************
!c        use the power series form
!c        *************************
!c
         z2=z*z
         sum=cone
         term=cone
         pterm=cone
         do while (&
               dabs(dreal(term)).gt.tol*dabs(dreal(sum)) .or.&
               dabs(dimag(term)).gt.tol*dabs(dimag(sum))&
                  )
            fn=fn+one
            pterm=pterm*z2/fn
            term=pterm/(fn+fn+one)
            sum=sum+term
         end do
         zfn=(irtpi-ctwo*z*sum)*cdexp(-z2)
         end if
      return
      end function zfn
!
!END MODULE Zfn_MODULE

END MODULE Sub_Prog

!---------------------------------------------------------------------
! Main program:
!---------------------------------------------------------------------

PROGRAM WT_LST
USE Common_Params
USE Common_Arrays
USE Math_Constants
USE Phys_Constants
USE Sub_Prog
IMPLICIT NONE
REAL*8 :: DTau,Tau,Tau1,Tau2,TauAdd
REAL*8 :: Dqx,Dqz,Dux,Duz,DQ
REAL*8 :: Ewave0,Ewave,Epart0,Epart,EppEw0,EppEw
REAL*8 :: EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3
REAL*8 :: Rn,Rp,Rw,Rs
REAL*8 :: Anorm0,Anorm
INTEGER :: Iflag
INTEGER :: it,Nitera,i
CHARACTER(LEN=6) :: WaveType
CHARACTER(LEN=3) :: TimeEvol
CHARACTER(LEN=3) :: OneDfiles,TwoDfiles,Onecolumnfiles,NewEffectsfiles
CHARACTER(LEN=3) :: ADcoefficients

OPEN(1,FILE='Start.wt')
READ(1,*) DTau
READ(1,*) TauAdd
READ(1,*) TimeEvol
READ(1,*) OneDfiles
READ(1,*) TwoDfiles
READ(1,*) Onecolumnfiles
READ(1,*) ADcoefficients
READ(1,*) NewEffectsfiles
CLOSE(1)

OPEN(2,FILE='Ratios.wt')
OPEN(3,FILE='Ewave.wt')
OPEN(1,FILE='Ini.wt')
READ(1,*) Iflag
READ(1,*) nux,nuz
READ(1,*) nqx,nqz
READ(1,*) nqx2,nqz2
READ(1,*) nqmod
READ(1,*) nrz
READ(1,*) Ulim
READ(1,*) Ucrit
READ(1,*) Qxi 
READ(1,*) Qxf 
READ(1,*) Qzi 
READ(1,*) Qzf 
READ(1,*) Qi
READ(1,*) Qf 
READ(1,*) InitialLevel 
READ(1,*) Iw0 
READ(1,*) AsympT 
READ(1,*) RatioNke 
READ(1,*) RatioNki 
READ(1,*) Kappae,Alphae
READ(1,*) Kappai,Alphai
READ(1,*) RatioNf 
READ(1,*) Uf 
READ(1,*) Thf 
READ(1,*) RTfparTs
READ(1,*) RTfperpTs
READ(1,*) RatioNb 
READ(1,*) Ub 
READ(1,*) Thb 
READ(1,*) RTbparTs
READ(1,*) RTbperpTs
READ(1,*) RTeTi 
READ(1,*) G
READ(1,*) Ve2C2
READ(1,*) Lemis
READ(1,*) LdecayLS
READ(1,*) LdecayLT
READ(1,*) LdecayST
READ(1,*) LdecayTT
READ(1,*) LscatLL
READ(1,*) LscatLT
READ(1,*) Semis
READ(1,*) SdecayLL
READ(1,*) SdecayLT
READ(1,*) Sscat
READ(1,*) TdecayLL
READ(1,*) TdecayLS
READ(1,*) TdecayTL
READ(1,*) TscatLT
READ(1,*) CollTerm
READ(1,*) CollTermForm
READ(1,*) SpontEmis
READ(1,*) ScatElSpo
!READ(1,*) ScatElInd
READ(1,*) NewEffects1
READ(1,*) NewEffects2
READ(1,*) Gcoll
READ(1,*) Bremss
READ(1,*) RenormFe
READ(1,*) DerivLn
READ(1,*) RebuildL
READ(1,*) RebuildS
READ(1,*) RebuildFe
READ(1,*) Lsaturation

IF(Qxi<Qmin .OR. Qzi<Qmin .OR. Qi<Qmin) THEN
 OPEN(98,FILE='Warning_Main_Qi.wt')
 WRITE(98,*) ' Qxi= ',Qxi,'  Qzi= ',Qzi,'  Qi= ',Qi
 WRITE(98,*) ' Qxi, Qzi, and Qi, can not be too small !!'
 CLOSE(98)
 IF(Qxi<Qmin) THEN
  Qxi= Qmin
 ELSE
 END IF
 IF(Qzi<Qmin) THEN
  Qzi= Qmin
 ELSE
 END IF
 IF(Qi<Qmin) THEN
  Qi= Qmin
 ELSE
 END IF
ELSE
END IF

IF(nqx*nqz>6561) THEN
 OPEN(98,FILE='Warning_Main_Size_RK4_Derivs.wt')
 WRITE(98,*) ' nqx= ',nqx,'  nqz= ',nqz,'  nqx*nqz= ',nqx*nqz
 WRITE(98,*) ' nqx*nqz can not be larger than NMAX in RK4 and DERIVS !!'
 CLOSE(98)
 STOP
ELSE
END IF
IF (nqz2>nqz) THEN
 nph= nqz2
ELSE
 nph= nqz 
END IF

IF(RatioNki >= EpsMin) THEN
 OPEN(98,FILE='Warning_Main_RatioNki.wt')
 WRITE(98,*) ' RatioNki ',RatioNki
 WRITE(98,*) ' RatioNki must be zero, in the present version !!'
 CLOSE(98)
 STOP
ELSE
END IF
CALL Allocate_Arrays
CALL Space_Profiles ! Just to generate profiles, which may be useful for
                    ! eventual extension to the case of inhomogeneous medium.
CALL Definitions

SELECT CASE(Iflag)

 CASE(0)
  CLOSE(1)
  CALL Init_Wave(Tau,Dqx,Dqz,Dux,Duz,Anorm0,Epart0,Ewave0,EppEw0,&
     EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3)
  CALL Output("Fe0")
!  CALL Output("Fi ")
  CALL Output("IL0")
  CALL Output("IS0")
  CALL Output("IT0")
  SELECT CASE(NewEffectsfiles)
   CASE("Yes")
    CALL Output("GcL")
    CALL Output("GqL")
    CALL Output("GcS")
    CALL Output("GqS")
    CALL Output("PbL")
    CALL Output("PbS")
   CASE("No ")
   CASE DEFAULT
    OPEN(98,FILE='Warning_WT_LST.wt')
    WRITE(98,*) ' NewEffectsfiles= ',NewEffectsfiles
    WRITE(98,*) ' NewEffectsfiles must be (Yes) or (No ) !! '
    CLOSE(98)
    STOP
  END SELECT

  it= 0
  Rn= 1.     ! Anorm/Anorm0
  Rp= 1.     ! Epart/Epart0
  Rw= 1.     ! Ewave/Ewave0
  Rs= 1.     ! EppEw/EppEw0
!  WRITE(2,2001) it,Tau,Rn,Rp,Rw,Rs
2001 FORMAT(1x,i5,7(1x,e13.5e3))
!  WRITE(3,2001) it,Tau,EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3

 CASE(1)
  CALL Read_Results(Tau,Dqx,Dqz,Dux,Duz,Anorm0,Epart0,Ewave0,EppEw0,&
       Rn,Rp,Rw,Rs,EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3)
  it= 0
!  WRITE(2,2001) it,Tau,Rn,Rp,Rw,Rs
!  WRITE(3,2001) it,Tau,EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3
! One-dimensional array (absolute value of Q):
  Dq= (Qf-Qi)/(nqmod-1)
  DO i= 1,nqmod
   VQQ(i)= Qi+(i-1)*Dq
  END DO
  VQQ(nqmod)= Qf
! Re-initializes the collisional damping and bremsstrahlung coefficients, 
! assumed to be constant along time evolution:
  CALL Coll_Damping
  CALL Bremsstrahlung

 CASE DEFAULT
  OPEN(98,FILE='Warning_Main.wt')
  WRITE(98,*) ' Iflag= ',Iflag
  WRITE(98,*) ' Iflag must be 0 or 1 !!'
  CLOSE(98)
  STOP

END SELECT

SELECT CASE(TimeEvol)

 CASE("Yes")
  ! Start time evolution:

  CALL Res_Cond
  CALL Coef_A
  CALL Coef_Coll
  Auxinterp= 1

  Tau1= Tau
  Tau2= Tau1+TauAdd
  Nitera= NINT((Tau2-Tau1)/DTau)

  !***     Time evolution:    ***
  DO it= 1,Nitera
   Tau= Tau+DTau
   ! Evolution of ILp:
    WaveType= "Lwavep"
    CALL Evol_Iwave(DTau,Tau,WaveType)
    CALL Correct_ILp(it)
   ! Evolution of ILm:
    WaveType= "Lwavem"
    CALL Evol_Iwave(DTau,Tau,WaveType)
    CALL Correct_ILm(it)
   ! Evolution of ISp:
    WaveType= "Swavep"
    CALL Evol_Iwave(DTau,Tau,WaveType)
   ! Evolution of ISm:
    WaveType= "Swavem"
    CALL Evol_Iwave(DTau,Tau,WaveType)
   ! Evolution of ITp:
    WaveType= "Twavep"
    CALL Evol_Iwave(DTau,Tau,WaveType)
   ! Evolution of ITm:
    WaveType= "Twavem"
    CALL Evol_Iwave(DTau,Tau,WaveType)

   CALL Split(Dux,Duz,DTau,it)

    ! IF (MOD(it, 10)==1) THEN
    !  !WRITE (*,*) "entrei na rebuild. it: ", it
    !  CALL Rebuild_Fe
    ! END IF

   Auxinterp= 0
   CALL Fnorm(Anorm)
   CALL Energy(Epart,Ewave,EppEw)
   Rn= Anorm/Anorm0
   Rp= Epart/Epart0
   Rw= Ewave/Ewave0
   Rs= EppEw/EppEw0
   WRITE(2,2001) it,Tau,Rn,Rp,Rw,Rs
   CALL Energy2(EwaveL,EWaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3)
   WRITE(3,2001) it,Tau,EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3
   OPEN(98,FILE='status.wt')
   WRITE(98,*) ' it= ',it,'   Tau= ',Tau
   CLOSE(98)
  END DO

  CLOSE(2)
  CLOSE(3)
  CLOSE(99)

  CALL Save_Results(Tau,Dqx,Dqz,Dux,Duz,Anorm0,Epart0,Ewave0,EppEw0,&
        Rn,Rp,Rw,Rs,EwaveL,EwaveS,EwaveT,EwaveTF,EwaveTH,EwaveT3)

  Iflag= 1

 CASE("No ")

  ! Nothing to be done at this point. Proceed to the next option.

 CASE DEFAULT
  OPEN(98,FILE='Warning_WT_LST.wt')
  WRITE(98,*) ' TimeEvol= ',TimeEvol
  WRITE(98,*) ' TimeEvol must be (Yes) or (No ) !! '
  CLOSE(98)
  STOP

END SELECT

IF (Iflag==1) THEN
 IF (RebuildL=="Yes") THEN
  CALL Rebuild_L
 ELSE
 END IF
 IF (RebuildS=="Yes") THEN
  CALL Rebuild_S
 ELSE
 END IF
 IF (RebuildFe=="Yes") THEN
  CALL Rebuild_Fe
 ELSE
 END IF

 SELECT CASE(TwoDfiles)
  CASE("Yes")
   CALL Output("Fe ")
   CALL Output("IL ")
   CALL Output("IS ")
   CALL Output("IT ")
   CALL Output("DDz")
  CASE("No ")
  CASE DEFAULT
   OPEN(98,FILE='Warning_WT_LST.wt')
   WRITE(98,*) ' TwoDfiles= ',TwoDfiles
   WRITE(98,*) ' TwoDfiles must be (Yes) or (No ) !! '
   CLOSE(98)
   STOP
 END SELECT
 SELECT CASE(OneDfiles)
  CASE("Yes")
   CALL Output("Fe1")
   CALL Output("IL1")
   CALL Output("IS1")
   CALL Output("IT1")
  CASE("No ")
  CASE DEFAULT
   OPEN(98,FILE='Warning_WT_LST.wt')
   WRITE(98,*) ' OneDfiles= ',OneDfiles
   WRITE(98,*) ' OneDfiles must be (Yes) or (No ) !! '
   CLOSE(98)
   STOP
 END SELECT
 SELECT CASE(Onecolumnfiles)
  CASE("Yes")
   CALL Output2("Ui ")
   CALL Output2("Qi ")
   CALL Output2("Q  ")
   CALL Output2("ZTQ")
   CALL Output2("Fe ")
   CALL Output2("IL ")
   CALL Output2("IS ")
   CALL Output2("IT ")
   CALL Output2("ILQ")
   CALL Output2("ISQ")
   CALL Output2("ITQ")
   CALL Output2("Fe1")
   CALL Output2("IL1")
   CALL Output2("IS1")
   CALL Output2("IT1")
  CASE("No ")
  CASE DEFAULT
   OPEN(98,FILE='Warning_WT_LST.wt')
   WRITE(98,*) ' Onecolumnfiles= ',Onecolumnfiles
   WRITE(98,*) ' Onecolumnfiles must be (Yes) or (No ) !! '
   CLOSE(98)
   STOP
 END SELECT
 SELECT CASE(ADcoefficients)
  CASE("Yes")
   CALL Output("Ai ")
   CALL Output("Dij")
  CASE("No ")
  CASE DEFAULT
   OPEN(98,FILE='Warning_WT_LST.wt')
   WRITE(98,*) ' ADCoefficients= ',ADcoefficients
   WRITE(98,*) ' ADcoefficients must be (Yes) or (No ) !! '
   CLOSE(98)
   STOP
 END SELECT
ELSE
 OPEN(98,FILE='Warning_WT_LST.wt')
 WRITE(98,*) ' Iflag= ',Iflag
 WRITE(98,*) ' Iflag must be (1) in order to generate output files!! '
 CLOSE(98)
 STOP
END IF

END PROGRAM WT_LST

! BELOW THIS, I HAVE THE ORIGINAL (FROM 13b) LUX AND LUZ OPERATORS:

 ! DO k= 1,nuz
 !  FORALL(i=2:nux-1)
 !   Auxi(i)= (TANH(SQRT((VUx(i))**2+(VUz(k))**2)))**2
 !   Alpha1(i)= ( ( CAx(i,k)+DDxxdux(i,k) ) * Auxx - CDxx(i,k)*Auxxx ) *Auxi(i)
 !   Beta1(i)= 1. + (- DAxdux(i,k)*DTau/2. + 2.*CDxx(i,k)*Auxxx ) *Auxi(i)
 !   Gamma1(i)= (- ( CAx(i,k)+DDxxdux(i,k) ) * Auxx - CDxx(i,k)*Auxxx ) *Auxi(i)
 !   Psi1(i)= Fe(i,k) + ( DAxdux(i,k)*Fe(i,k)*(DTau/2.) &
 !        + ( CDxz(i,k)*Dfduxduz(i,k) + DDxzdux(i,k)*Dfduz(i,k) ) * (DTau) &
 !        + ( CAx(i,k)+DDxxdux(i,k) ) * Auxx * (Fe(i+1,k)-Fe(i-1,k)) &
 !        + CDxx(i,k)*Auxxx * (Fe(i+1,k)-2.*Fe(i,k)+Fe(i-1,k)) ) *Auxi(i)
 !   Faux1b(i)= Fe(i,k)
 !  END FORALL
 !  i= 2
 !  Psi1(i)= Psi1(i)-Alpha1(i)*Fe(i-1,k)
 !  Alpha1(i)= 0.
 !  Faux1b(i)= Fe(i,k)
 !  i= nux
 !  Auxi(i)= (TANH(SQRT((VUx(i))**2+(VUz(k))**2)))**2
 !  Alpha1(i)= ( CAx(i,k)+DDxxdux(i,k) ) * (DTau/2./Dux) *Auxi(i)
 !  Beta1(i)= 1. + (- DAxdux(i,k)*DTau/2. &
 !       - ( CAx(i,k)+DDxxdux(i,k) ) * (DTau/2./Dux) ) *Auxi(i)
 !  Gamma1(i)= 0.
 !  Psi1(i)= Fe(i,k) + ( DAxdux(i,k)*Fe(i,k)*(DTau/2.) &
 !       + ( CDxz(i,k)*Dfduxduz(i,k) + DDxzdux(i,k)*Dfduz(i,k) ) * (DTau) &
 !       + ( CAx(i,k)+DDxxdux(i,k) ) * (DTau/2./Dux) * (Fe(i,k)-Fe(i-1,k)) &
 !       + CDxx(i,k)*(DTau/Dux/Dux) * (Fe(i,k)-2.*Fe(i-1,k)+Fe(i-2,k)) ) *Auxi(i)
 !  Faux1b(i)= Fe(i,k)