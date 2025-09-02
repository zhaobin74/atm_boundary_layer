#include "unused_dummy.H"

      program surf_layer 
 
      use icepack_kinds
      use icepack_parameters, only: c0, c1, c2, c4, c5, c8, c10
      use icepack_parameters, only: c16, c20, p001, p01, p2, p4, p5, p75, puny
      use icepack_parameters, only: senscoef, latncoef
      use icepack_parameters, only: cp_wv, cp_air, iceruf, zref, qqqice, TTTice, qqqocn, TTTocn
      use icepack_parameters, only: Lsub, Lvap, vonkar, Tffresh, zvir, gravit
      use icepack_parameters, only: pih, dragio, rhoi, rhos, rhow
      use icepack_parameters, only: atmbndy, calc_strair, formdrag
      !use icepack_parameters, only: icepack_chkoptargflag

      implicit none

      integer, parameter :: nx=101

       real, parameter              :: MAPL_GRAV                      = 9.80665  

 real, parameter              :: MAPL_RUNIV                     = 8314.47                        ! J/(Kmole K)


    ! Physical properties
   real, parameter              :: MAPL_H2OMW                     =  18.015                        ! kg/Kmole
   real, parameter              :: MAPL_O3MW                      = 47.9982                        ! kg/Kmole
   real, parameter              :: MAPL_LATENT_HEAT_VAPORIZATION  = 2.4665E6                       ! J/kg @15C @1atm
   real, parameter              :: MAPL_ALHL                      = MAPL_LATENT_HEAT_VAPORIZATION  ! J/kg
   real, parameter              :: MAPL_LATENT_HEAT_FUSION        = 3.3370E5                       ! J/kg @1atm
   real, parameter              :: MAPL_ALHF                      = MAPL_LATENT_HEAT_FUSION        ! J/kg
   real, parameter              :: MAPL_LATENT_HEAT_SUBLIMATION   = MAPL_ALHL+MAPL_ALHF            ! J/kg
   real, parameter              :: MAPL_ALHS                      = MAPL_LATENT_HEAT_SUBLIMATION   ! J/kg

   ! Earth Specific Chemistry and Thermodynamic Constants
   real, parameter              :: MAPL_AIRMW                     =  28.965                        ! kg/Kmole
   real, parameter              :: MAPL_RDRY                      = MAPL_RUNIV/MAPL_AIRMW          ! J/(kg K)
   real, parameter              :: MAPL_CPDRY                     = 3.5*MAPL_RDRY                  ! J/(kg K)
   real, parameter              :: MAPL_CVDRY                     = MAPL_CPDRY-MAPL_RDRY           ! J/(kg K)
   real, parameter              :: MAPL_RVAP                      = MAPL_RUNIV/MAPL_H2OMW          ! J/(kg K)
   real, parameter              :: MAPL_CPVAP                     = 4.*MAPL_RVAP                   ! J/(kg K)
   real, parameter              :: MAPL_CVVAP                     = MAPL_CPVAP-MAPL_RVAP           ! J/(kg K)
   real, parameter              :: MAPL_KAPPA                     = MAPL_RDRY/MAPL_CPDRY           ! (2.0/7.0)

      real, parameter              :: MAPL_EPSILON                   = MAPL_H2OMW/MAPL_AIRMW 

      real, parameter              :: MAPL_VIREPS                    = 1.0/MAPL_EPSILON-1.0

      real, parameter              :: MAPL_KARMAN                    = 0.40
      real, parameter              :: MAPL_RGAS                      = MAPL_RDRY
      real, parameter              :: MAPL_CP                        = MAPL_RGAS/MAPL_KAPPA
      real, parameter              :: MAPL_NUAIR                     = 1.533E-5  

      real*4  :: Z(nx)
      real*4  :: Phim(nx), Phih(nx)
      integer i, j, n, iz, irec


      do i=1,nx
         Z(i) = real(i-1, kind=4)
      enddo   

      call PHI(Z, Phim, Phih, 1, nx)

      do i=1,nx
         print*, Z(i), phim(i)
      enddo
       
contains

      


!=======================================================================

! Atmospheric boundary interface (stability based flux calculations)

! author: Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb
! 2004: Block structure added by William Lipscomb
! 2006: Converted to free source form (F90) by Elizabeth Hunke
! 2013: Form drag routine added (neutral_drag_coeffs) by David Schroeder
! 2014: Adjusted form drag and added high frequency coupling by Andrew Roberts


      !use icepack_parameters, only: c0, c1, c2, c4, c5, c8, c10
      !use icepack_parameters, only: c16, c20, p001, p01, p2, p4, p5, p75, puny
      !use icepack_parameters, only: senscoef, latncoef
      !use icepack_parameters, only: cp_wv, cp_air, iceruf, zref, qqqice, TTTice, qqqocn, TTTocn
      !use icepack_parameters, only: Lsub, Lvap, vonkar, Tffresh, zvir, gravit
      !use icepack_parameters, only: pih, dragio, rhoi, rhos, rhow
      !use icepack_parameters, only: atmbndy, calc_strair, formdrag
      !use icepack_parameters, only: icepack_chkoptargflag
      !use icepack_tracers, only: n_iso
      !use icepack_tracers, only: tr_iso
      !use icepack_warnings, only: warnstr, icepack_warnings_add
      !use icepack_warnings, only: icepack_warnings_setabort, icepack_warnings_aborted



!=======================================================================

!=======================================================================

! Compute coefficients for atm/ice fluxes, stress, and reference
! temperature and humidity. NOTE:
! (1) All fluxes are positive downward,
! (2) Here, tstar = (WT)/U*, and qstar = (WQ)/U*,
! (3a) wind speeds should all be above a minimum speed (eg. 1.0 m/s).
!
! ASSUME:
!  The saturation humidity of air at T(K): qsat(T)  (kg/m**3)
!
! Code originally based on CSM1

      subroutine atmo_boundary_layer (sfctype,            &
                                      calc_strair, formdrag, &
                                      Tsf,      potT,     &
                                      uatm,     vatm,     &
                                      wind,     zlvl,     &
                                      Qa,       rhoa,     &
                                      strx,     stry,     &
                                      Tref,     Qref,     &
                                      delt,     delq,     &
                                      lhcoef,   shcoef,   &
                                      Cdn_atm,            &
                                      Cdn_atm_ratio_n,    &
                                      Qa_iso,   Qref_iso, &
                                      uvel,     vvel,     &
                                      Uref,     zlvs      )

      use icepack_parameters, only: highfreq, natmiter, atmiter_conv

      character (len=3), intent(in) :: &
         sfctype      ! ice or ocean

      logical (kind=log_kind), intent(in) :: &
         calc_strair, &  ! if true, calculate wind stress components
         formdrag        ! if true, calculate form drag

      real (kind=dbl_kind), intent(in) :: &
         Tsf      , & ! surface temperature of ice or ocean
         potT     , & ! air potential temperature  (K)
         uatm     , & ! x-direction wind speed (m/s)
         vatm     , & ! y-direction wind speed (m/s)
         wind     , & ! wind speed (m/s)
         zlvl     , & ! atm level height (m)
         Qa       , & ! specific humidity (kg/kg)
         rhoa         ! air density (kg/m^3)

      real (kind=dbl_kind), intent(inout) :: &
         Cdn_atm      ! neutral drag coefficient

      real (kind=dbl_kind), intent(inout) :: &
         Cdn_atm_ratio_n ! ratio drag coeff / neutral drag coeff

      real (kind=dbl_kind), intent(inout) :: &
         strx     , & ! x surface stress (N)
         stry         ! y surface stress (N)

      real (kind=dbl_kind), intent(inout) :: &
         Tref     , & ! reference height temperature  (K)
         Qref     , & ! reference height specific humidity (kg/kg)
         delt     , & ! potential T difference   (K)
         delq     , & ! humidity difference      (kg/kg)
         shcoef   , & ! transfer coefficient for sensible heat
         lhcoef       ! transfer coefficient for latent heat

      real (kind=dbl_kind), intent(in), dimension(:), optional :: &
         Qa_iso       ! specific isotopic humidity (kg/kg)

      real (kind=dbl_kind), intent(inout), dimension(:), optional :: &
         Qref_iso     ! reference specific isotopic humidity (kg/kg)

      real (kind=dbl_kind), intent(in) :: &
         uvel     , & ! x-direction ice speed (m/s)
         vvel         ! y-direction ice speed (m/s)

      real (kind=dbl_kind), intent(out) :: &
         Uref         ! reference height wind speed (m/s)

      real (kind=dbl_kind), intent(in), optional :: &
         zlvs        ! atm level height (scalar quantities) (m)

      ! local variables

      integer (kind=int_kind) :: &
         k,n         ! iteration index

      real (kind=dbl_kind) :: &
         TsfK  , & ! surface temperature in Kelvin (K)
         psimh , & ! stability function at zlvl   (momentum)
         tau   , & ! stress at zlvl
         fac   , & ! interpolation factor
         al2   , & ! ln(z10   /zTrf)
         psix2 , & ! stability function at zTrf   (heat and water)
         ssq   , & ! sat surface humidity     (kg/kg)
         qqq   , & ! for qsat, dqsfcdt
         TTT   , & ! for qsat, dqsfcdt
         qsat  , & ! the saturation humidity of air (kg/m^3)
         Lheat , & ! Lvap or Lsub, depending on surface type
         umin      ! minimum wind speed (m/s)

      real (kind=dbl_kind) :: &
         ustar , & ! ustar (m/s)
         ustar_prev , & ! ustar_prev (m/s)
         tstar , & ! tstar
         qstar , & ! qstar
         ratio , & ! ratio
         rdn   , & ! sqrt of neutral exchange coefficient (momentum)
         rhn   , & ! sqrt of neutral exchange coefficient (heat)
         ren   , & ! sqrt of neutral exchange coefficient (water)
         rd    , & ! sqrt of exchange coefficient (momentum)
         re    , & ! sqrt of exchange coefficient (water)
         rh    , & ! sqrt of exchange coefficient (heat)
         vmag  , & ! surface wind magnitude   (m/s)
         alzm  , & ! ln(zlvl  /z10)
         alzs  , & ! ln(zlvs  /z10) (if zlvs present)
         thva  , & ! virtual temperature      (K)
         cp    , & ! specific heat of moist air
         holm  , & ! H (at zlvl  ) over L
         hols  , & ! H (at zlvs  ) over L (if zlvs present)
         stable, & ! stability factor
         cpvir , & ! defined as cp_wv/cp_air - 1.
         psixh     ! stability function at zlvl (at zlvs if present) (heat and water)

      real (kind=dbl_kind), parameter :: &
         zTrf  = c2     ! reference height for air temp (m)

      character(len=*),parameter :: subname='(atmo_boundary_layer)'

      al2 = log(zref/zTrf)

      !------------------------------------------------------------
      ! Initialize
      !------------------------------------------------------------

      cpvir = cp_wv/cp_air-c1   ! defined as cp_wv/cp_air - 1.

      if (highfreq) then
       umin  = p5 ! minumum allowable wind-ice speed difference of 0.5 m/s
      else
       umin  = c1 ! minumum allowable wind speed of 1m/s
      endif

      Tref = c0
      Qref = c0
      Uref = c0
      delt = c0
      delq = c0
      shcoef = c0
      lhcoef = c0

      !------------------------------------------------------------
      ! Compute turbulent flux coefficients, wind stress, and
      ! reference temperature and humidity.
      !------------------------------------------------------------

      !------------------------------------------------------------
      ! define variables that depend on surface type
      !------------------------------------------------------------

      if (sfctype(1:3)=='ice') then

         qqq   = qqqice          ! for qsat
         TTT   = TTTice          ! for qsat
         Lheat = Lsub            ! ice to vapor

            if (highfreq) then
               vmag = max(umin, sqrt( (uatm-uvel)**2 + &
                                      (vatm-vvel)**2) )
            else
               vmag = max(umin, wind)
            endif

            if (formdrag .and. Cdn_atm > puny) then
               rdn = sqrt(Cdn_atm)
            else
               rdn  = vonkar/log(zref/iceruf) ! neutral coefficient
               Cdn_atm = rdn * rdn
            endif

      elseif (sfctype(1:3)=='ocn') then

         qqq   = qqqocn
         TTT   = TTTocn
         Lheat = Lvap           ! liquid to vapor
         vmag = max(umin, wind)
         rdn  = sqrt(0.0027_dbl_kind/vmag &
              + .000142_dbl_kind + .0000764_dbl_kind*vmag)

      endif   ! sfctype

      !------------------------------------------------------------
      ! define some more needed variables
      !------------------------------------------------------------

      TsfK = Tsf + Tffresh     ! surface temp (K)
      delt = potT - TsfK       ! pot temp diff (K)
      qsat = qqq * exp(-TTT/TsfK)   ! saturation humidity (kg/m^3)
      ssq  = qsat / rhoa       ! sat surf hum (kg/kg)

      thva = potT * (c1 + zvir * Qa) ! virtual pot temp (K)
      delq = Qa - ssq          ! spec hum dif (kg/kg)
      alzm = log(zlvl/zref)
      if (present(zlvs)) then
         alzs = log(zlvs/zref)
      else
         alzs = alzm
      endif
      cp   = cp_air*(c1 + cpvir*ssq)

      !------------------------------------------------------------
      ! first estimate of Z/L and ustar, tstar and qstar
      !------------------------------------------------------------

      ! neutral coefficients, z/L = 0.0
      rhn = rdn
      ren = rdn

      ! ustar,tstar,qstar
      ustar = rdn * vmag
      tstar = rhn * delt
      qstar = ren * delq

      !------------------------------------------------------------
      ! iterate to converge on Z/L, ustar, tstar and qstar
      !------------------------------------------------------------

      ustar_prev = c2 * ustar

      k = 1
      do while (abs(ustar - ustar_prev)/ustar > atmiter_conv .and. k <= natmiter)
         k = k + 1
         ustar_prev = ustar

         ! compute stability & evaluate all stability functions
         holm = compute_stability_parameter(zlvl , thva , &
                                           ustar, tstar, &
                                           qstar, Qa)
         if (present(zlvs)) then
            hols = compute_stability_parameter(zlvs , thva , &
                                               ustar, tstar, &
                                               qstar, Qa)
         else
            hols = holm
         endif

         call compute_stability_function('momentum', holm, stable, psimh)
         call compute_stability_function('scalar'  , hols, stable, psixh)

         ! shift all coeffs to measurement height and stability
         rd = rdn / (c1+rdn/vonkar*(alzm-psimh))
         rh = rhn / (c1+rhn/vonkar*(alzs-psixh))
         re = ren / (c1+ren/vonkar*(alzs-psixh))

         ! update ustar, tstar, qstar using updated, shifted coeffs
         ustar = rd * vmag
         tstar = rh * delt
         qstar = re * delq

      enddo                     ! end iteration

      if (calc_strair) then

         ! initialize
         strx = c0
         stry = c0

         if (highfreq .and. sfctype(1:3)=='ice') then

            !------------------------------------------------------------
            ! momentum flux for high frequency coupling (RASM/CESM)
            !------------------------------------------------------------
            ! tau = rhoa * rd * rd
            ! strx = tau * |Uatm-U| * (uatm-u)
            ! stry = tau * |Uatm-U| * (vatm-v)
            !------------------------------------------------------------

            tau = rhoa * rd * rd ! not the stress at zlvl

            ! high frequency momentum coupling following Roberts et al. (2014)
            strx = tau * sqrt((uatm-uvel)**2 + (vatm-vvel)**2) * (uatm-uvel)
            stry = tau * sqrt((uatm-uvel)**2 + (vatm-vvel)**2) * (vatm-vvel)

         else

            !------------------------------------------------------------
            ! momentum flux
            !------------------------------------------------------------
            ! tau = rhoa * ustar * ustar
            ! strx = tau * uatm / vmag
            ! stry = tau * vatm / vmag
            !------------------------------------------------------------

            tau = rhoa * ustar * rd ! not the stress at zlvl
            strx = tau * uatm
            stry = tau * vatm

         endif

         Cdn_atm_ratio_n = rd * rd / rdn / rdn

      endif                     ! calc_strair

      !------------------------------------------------------------
      ! coefficients for turbulent flux calculation
      !------------------------------------------------------------
      ! add windless coefficient for sensible heat flux
      ! as in Jordan et al (JGR, 1999)
      !------------------------------------------------------------

      if (trim(atmbndy) == 'mixed') then
         !- Use constant coefficients for sensible and latent heat fluxes
         !    similar to atmo_boundary_const but using vmag instead of wind
         shcoef = senscoef*cp_air*rhoa*vmag
         lhcoef = latncoef*Lheat *rhoa*vmag
      else ! 'similarity'
         !- Monin-Obukhov similarity theory for boundary layer
         shcoef = rhoa * ustar * cp * rh + c1
         lhcoef = rhoa * ustar * Lheat  * re
      endif

      !------------------------------------------------------------
      ! Compute diagnostics: 2m ref T, Q, U
      !------------------------------------------------------------

      hols  = hols*zTrf/zlvl
      psix2 = -c5*hols*stable + (c1-stable)*psi_scalar_unstable(hols)
      fac   = (rh/vonkar) &
            * (alzs + al2 - psixh + psix2)
      Tref  = potT - delt*fac
      Tref  = Tref - p01*zTrf ! pot temp to temp correction
      fac   = (re/vonkar) &
            * (alzs + al2 - psixh + psix2)
      Qref  = Qa - delq*fac

      if (highfreq .and. sfctype(1:3)=='ice') then
         Uref = sqrt((uatm-uvel)**2 + (vatm-vvel)**2) * rd / rdn
      else
         Uref = vmag * rd / rdn
      endif


      end subroutine atmo_boundary_layer

!=======================================================================

! Compute coefficients for atm/ice fluxes, stress
! NOTE: \\
! (1) all fluxes are positive downward,  \\
! (2) reference temperature and humidity are NOT computed

      subroutine atmo_boundary_const (sfctype,  calc_strair, &
                                      uatm,     vatm,     &
                                      wind,     rhoa,     &
                                      strx,     stry,     &
                                      Tsf,      potT,     &
                                      Qa,                 &
                                      delt,     delq,     &
                                      lhcoef,   shcoef    )

      character (len=3), intent(in) :: &
         sfctype      ! ice or ocean

      logical (kind=log_kind), intent(in) :: &
         calc_strair  ! if true, calculate wind stress components

      real (kind=dbl_kind), intent(in) :: &
         Tsf      , & ! surface temperature of ice or ocean
         potT     , & ! air potential temperature  (K)
         Qa       , & ! specific humidity (kg/kg)
         uatm     , & ! x-direction wind speed (m/s)
         vatm     , & ! y-direction wind speed (m/s)
         wind     , & ! wind speed (m/s)
         rhoa         ! air density (kg/m^3)

      real (kind=dbl_kind), intent(inout):: &
         strx     , & ! x surface stress (N)
         stry         ! y surface stress (N)

      real (kind=dbl_kind), intent(out):: &
         delt     , & ! potential T difference   (K)
         delq     , & ! humidity difference      (kg/kg)
         shcoef   , & ! transfer coefficient for sensible heat
         lhcoef       ! transfer coefficient for latent heat

       ! local variables

      real (kind=dbl_kind) :: &
         TsfK, & ! surface temperature in Kelvin (K)
         qsat, & ! the saturation humidity of air (kg/m^3)
         ssq , & ! sat surface humidity     (kg/kg)
         tau, &  ! stress at zlvl
         Lheat   ! Lvap or Lsub, depending on surface type

      character(len=*),parameter :: subname='(atmo_boundary_const)'

      !------------------------------------------------------------
      ! Initialize
      !------------------------------------------------------------

      delt = c0
      delq = c0
      shcoef = c0
      lhcoef = c0

      if (calc_strair) then

         strx = c0
         stry = c0

      !------------------------------------------------------------
      ! momentum flux
      !------------------------------------------------------------
         tau = rhoa * 0.0012_dbl_kind * wind
!AOMIP         tau = rhoa * (1.10_dbl_kind + c4*p01*wind) &
!AOMIP                         * wind * p001
         strx = tau * uatm
         stry = tau * vatm

      endif                     ! calc_strair

      !------------------------------------------------------------
      ! define variables that depend on surface type
      !------------------------------------------------------------

      if (sfctype(1:3)=='ice') then
         Lheat = Lsub           ! ice to vapor
      elseif (sfctype(1:3)=='ocn') then
         Lheat = Lvap           ! liquid to vapor
      endif   ! sfctype

      !------------------------------------------------------------
      ! potential temperature and specific humidity differences
      !------------------------------------------------------------

      TsfK     = Tsf + Tffresh    ! surface temp (K)
      qsat     = qqqocn * exp(-TTTocn/TsfK) ! sat humidity (kg/m^3)
      ssq      = qsat / rhoa      ! sat surf hum (kg/kg)

      delt= potT - TsfK      ! pot temp diff (K)
      delq= Qa - ssq         ! spec hum dif (kg/kg)

      !------------------------------------------------------------
      ! coefficients for turbulent flux calculation
      !------------------------------------------------------------

      shcoef = senscoef*cp_air*rhoa*wind
      lhcoef = latncoef*Lheat *rhoa*wind

      end subroutine atmo_boundary_const

!=======================================================================

! Neutral drag coefficients for ocean and atmosphere also compute the
! intermediate necessary variables ridge height, distance, floe size
! based upon Tsamados et al. (2014), JPO, DOI: 10.1175/JPO-D-13-0215.1.
! Places where the code varies from the paper are commented.
!
! authors: Michel Tsamados, CPOM
!          David Schroeder, CPOM
!
! changes: Andrew Roberts, NPS (RASM/CESM coupling and documentation)

      subroutine neutral_drag_coeffs (apnd,     hpnd,     &
                                      ipnd,               &
                                      alvl,     vlvl,     &
                                      aice,     vice,     &
                                      vsno,     aicen,    &
                                      vicen, &
                                      Cdn_ocn,  Cdn_ocn_skin,    &
                                      Cdn_ocn_floe, Cdn_ocn_keel,&
                                      Cdn_atm,  Cdn_atm_skin,    &
                                      Cdn_atm_floe, Cdn_atm_pond,&
                                      Cdn_atm_rdg, hfreebd,      &
                                      hdraft,   hridge,          &
                                      distrdg,  hkeel,           &
                                      dkeel,    lfloe,           &
                                      dfloe,    ncat)

      !use icepack_tracers, only: tr_pond

      integer (kind=int_kind), intent(in) :: &
         ncat

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         apnd     ,& ! melt pond fraction of sea ice
         hpnd     ,& ! mean melt pond depth over sea ice
         ipnd     ,& ! mean ice pond depth over sea ice in cat n
         alvl     ,& ! level ice area fraction (of grid cell ?)
         vlvl        ! level ice mean thickness

      real (kind=dbl_kind), intent(in) :: &
         aice     , & ! concentration of ice
         vice     , & ! volume per unit area of ice
         vsno         ! volume per unit area of snow

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         aicen    , & ! concentration of ice
         vicen        ! volume per unit area of ice (m)

      real (kind=dbl_kind), intent(out) :: &
         hfreebd      , & ! freeboard (m)
         hdraft       , & ! draught of ice + snow column (Stoessel1993)
         hridge       , & ! ridge height
         distrdg      , & ! distance between ridges
         hkeel        , & ! keel depth
         dkeel        , & ! distance between keels
         lfloe        , & ! floe length (m)
         dfloe        , & ! distance between floes
         Cdn_ocn      , & ! ocean-ice neutral drag coefficient
         Cdn_ocn_skin , & ! drag coefficient due to skin drag
         Cdn_ocn_floe , & ! drag coefficient due to floe edges
         Cdn_ocn_keel , & ! drag coefficient due to keels
         Cdn_atm      , & ! ice-atmosphere drag coefficient
         Cdn_atm_skin , & ! drag coefficient due to skin drag
         Cdn_atm_floe , & ! drag coefficient due to floe edges
         Cdn_atm_pond , & ! drag coefficient due to ponds
         Cdn_atm_rdg      ! drag coefficient due to ridges

      real (kind=dbl_kind), parameter :: &
                                      ! [,] = range of values that can be tested
         csw       = 0.002_dbl_kind ,&! ice-ocn drag coefficient [0.0005,0.005]
         csa       = 0.0005_dbl_kind,&! ice-air drag coefficient [0.0001,0.001]
         mrdg      = c20            ,&! screening effect see Lu2011 [5,50]
         mrdgo     = c10            ,&! screening effect see Lu2011 [5,50]
         beta      = p5             ,&! power exponent appearing in astar and
                                      ! L=Lmin(A*/(A*-A))**beta [0,1]
         Lmin      = c8             ,&! min length of floe (m) [5,100]
         Lmax      = 300._dbl_kind  ,&! max length of floe (m) [30,3000]
         cfa       = p2             ,&! Eq. 12 ratio of local from drag over
                                      ! geometrical parameter [0,1]
         cfw       = p2             ,&! Eq. 15 ratio of local from drag over
                                      ! geometrical parameter [0,1]
         cpa       = p2             ,&! Eq. 16 ratio of local form drag over
                                      ! geometrical parameter [0,1]
         cra       = p2             ,&! Eq. 10 local form drag coefficient [0,1]
         crw       = p2             ,&! Eq. 11 local form drag coefficient [0,1]
         sl        = 22._dbl_kind   ,&! Sheltering parameter Lupkes2012 [10,30]
         lpmin     = 2.26_dbl_kind  ,&! min pond length (m) see Eq. 17 [1,10]
         lpmax     = 24.63_dbl_kind ,&! max pond length (m) see Eq. 17 [10,100]
         tanar     = p4             ,&! 0.25 sail slope = 14 deg [0.4,1]
         tanak     = p4             ,&! 0.58 keel slope = 30 deg [0.4,1]
         phir      = 0.8_dbl_kind   ,&! porosity of ridges [0.4,1]
         phik      = 0.8_dbl_kind   ,&! porosity of keels  [0.4,1]
         hkoverhr  = c4             ,&! hkeel/hridge ratio [4,8]
         dkoverdr  = c1             ,&! dkeel/distrdg ratio [1,5]
         sHGB      = 0.18_dbl_kind  ,&! Lupkes2012 Eq. 28, Hanssen1988,
                                      ! Steele1989 suggest instead 0.18
         alpha2    = c0             ,&! weight functions for area of
         beta2     = p75              ! ridged ice [0,1]

       integer (kind=int_kind) :: &
         n            ! category index

      real (kind=dbl_kind) :: &
         astar,     & ! new constant for form drag
         ctecaf,    & ! constante
         ctecwf,    & ! constante
         sca,       & ! wind attenuation function
         scw,       & ! ocean attenuation function
         lp,        & ! pond length (m)
         ctecar,    &
         ctecwk,    &
         ai, aii,   & ! ice area and its inverse
         ocnrufi,   & ! inverse ocean roughness
         icerufi,   & ! inverse ice roughness
         tmp1         ! temporary

      real (kind=dbl_kind) :: &
         apond    , & ! melt pond fraction of grid cell
         ardg     , & ! ridged ice area fraction of grid cell
         vrdg         ! ridged ice mean thickness

      real (kind=dbl_kind), parameter :: &
         ocnruf   = 0.000327_dbl_kind ! ocean surface roughness (m)

      real (kind=dbl_kind), parameter :: &
         camax    = 0.02_dbl_kind , & ! Maximum for atmospheric drag
         cwmax    = 0.06_dbl_kind     ! Maximum for ocean drag

      logical :: tr_pond

      character(len=*),parameter :: subname='(neutral_drag_coeffs)'

      astar = c1/(c1-(Lmin/Lmax)**(c1/beta))

      !-----------------------------------------------------------------
      ! Initialize across entire grid
      !-----------------------------------------------------------------
    
      tr_pond = .false.

      ocnrufi  = c1/ocnruf    ! inverse ocean roughness
      icerufi  = c1/iceruf    ! inverse ice roughness
      hfreebd=c0
      hdraft =c0
      hridge =c0
      distrdg=c0
      hkeel  =c0
      dkeel  =c0
      lfloe  =c0
      dfloe  =c0
      Cdn_ocn=dragio
      Cdn_ocn_skin=c0
      Cdn_ocn_floe=c0
      Cdn_ocn_keel=c0
      Cdn_atm = (vonkar/log(zref/iceruf)) * (vonkar/log(zref/iceruf))
      Cdn_atm_skin=c0
      Cdn_atm_floe=c0
      Cdn_atm_pond=c0
      Cdn_atm_rdg =c0

      if (aice > p001) then

         Cdn_atm_skin = csa
         Cdn_ocn_skin = csw

         ai  = aice
         aii = c1/ai

      !------------------------------------------------------------
      ! Compute average quantities
      !------------------------------------------------------------

         ! ponds
         apond = c0
         if (tr_pond) then
            do n = 1,ncat
               ! area of pond per unit area of grid cell
               apond = apond+apnd(n)*aicen(n)
            enddo
         endif

         ! draft and freeboard (see Eq. 27)
         hdraft = (rhoi*vice+rhos*vsno)*aii/rhow ! without ponds
         hfreebd = (vice+vsno)*aii-hdraft

         ! Do not allow draft larger than ice thickness (see Eq. 28)
         if (hdraft >= vice*aii) then
            ! replace excess snow with ice so hi~=hdraft
            hfreebd = (hdraft*ai*(c1-rhoi/rhow) + &
                      (vsno-(vice-hdraft*ai)*rhoi/rhos) * &
                      (c1-rhos/rhow))*aii ! Stoessel1993
         endif

         ! floe size parameterization see Eq. 13
         lfloe = Lmin * (astar / (astar - ai))**beta

         ! distance between floes parameterization see Eq. 14
         dfloe = lfloe * (c1/sqrt(ai) - c1)

         ! Relate ridge height and distance between ridges to
         ! ridged ice area fraction and ridged ice mean thickness
         ! Assumes total volume of ridged ice is split into ridges and keels.
         ! Then assume total ridges volume over total area of ridges =
         ! volume of one average ridge / area of one average ridge
         ! Same for keels.

         ardg=c0
         vrdg=c0
         do n=1,ncat
            ! ridged ice area fraction over grid cell
            ardg=ardg+(c1-alvl(n))*aicen(n)
            ! total ridged ice volume per unit grid cell area
            vrdg=vrdg+(c1-vlvl(n))*vicen(n)
         enddo

         ! hridge, hkeel, distrdg and dkeel estimates from CICE for
         ! simple triangular geometry
         if (ardg > p001) then
            ! see Eq. 25 and Eq. 26
            hridge = vrdg/ardg*c2 &
                   * (alpha2+beta2*hkoverhr/dkoverdr*tanar/tanak) &
                   / (phir*c1+phik*tanar/tanak*hkoverhr**c2/dkoverdr)
            distrdg = c2*hridge*ai/ardg &
                    * (alpha2/tanar+beta2/tanak*hkoverhr/dkoverdr)
            hkeel = hkoverhr * hridge
            dkeel = dkoverdr * distrdg

          ! Use the height of ridges relative to the mean freeboard of
          ! the pack.  Therefore skin drag and ridge drag differ in
          ! this code as compared to  Tsamados et al. (2014) equations
          ! 10 and 18, which reference both to sea level.
          tmp1 = max(c0,hridge - hfreebd)

      !------------------------------------------------------------
      ! Skin drag (atmo)
      !------------------------------------------------------------

          Cdn_atm_skin = csa*(c1 - mrdg*tmp1/distrdg)
          Cdn_atm_skin = max(min(Cdn_atm_skin,camax),c0)

      !------------------------------------------------------------
      ! Ridge effect (atmo)
      !------------------------------------------------------------

          if (tmp1 > puny) then
            sca = c1 - exp(-sHGB*distrdg/tmp1) ! see Eq. 9
            ctecar = cra*p5
            Cdn_atm_rdg = ctecar*tmp1/distrdg*sca* &
                       (log(tmp1*icerufi)/log(zref*icerufi))**c2
            Cdn_atm_rdg = min(Cdn_atm_rdg,camax)
          endif

          ! Use the depth of keels relative to the mean draft of
          ! the pack.  Therefore skin drag and keel drag differ in
          ! this code as compared to  Tsamados et al. (2014) equations
          ! 11 and 19, which reference both to  sea level. In some
          ! circumstances, hkeel can be less than hdraft because hkoverhr
          ! is constant, and max(c0,...) temporarily addresses this.
          tmp1 = max(c0,hkeel - hdraft)

      !------------------------------------------------------------
      ! Skin drag bottom ice (ocean)
      !------------------------------------------------------------

          Cdn_ocn_skin = csw * (c1 - mrdgo*tmp1/dkeel)
          Cdn_ocn_skin = max(min(Cdn_ocn_skin,cwmax), c0)

      !------------------------------------------------------------
      ! Keel effect (ocean)
      !------------------------------------------------------------

          if (tmp1 > puny) then
            scw = c1 - exp(-sHGB*dkeel/tmp1)
            ctecwk = crw*p5
            Cdn_ocn_keel = ctecwk*tmp1/dkeel*scw* &
                        (log(tmp1*icerufi)/log(zref*icerufi))**c2
            Cdn_ocn_keel = max(min(Cdn_ocn_keel,cwmax),c0)
          endif

         endif ! ardg > 0.001

      !------------------------------------------------------------
      ! Floe edge drag effect (atmo)
      !------------------------------------------------------------

        if (hfreebd > puny) then
          sca = c1 - exp(-sl*beta*(c1-ai))
          ctecaf = cfa*p5*(log(hfreebd*ocnrufi)/log(zref*ocnrufi))**c2*sca
          Cdn_atm_floe = ctecaf * hfreebd / lfloe
          Cdn_atm_floe = max(min(Cdn_atm_floe,camax),c0)
        endif

      !------------------------------------------------------------
      ! Pond edge effect (atmo)
      !------------------------------------------------------------

        if (hfreebd > puny) then
          sca = (apond)**(c1/(zref*beta))
          lp  = lpmin*(1-apond)+lpmax*apond
          Cdn_atm_pond = cpa*p5*sca*apond*hfreebd/lp &
                   * (log(hfreebd*ocnrufi)/log(zref*ocnrufi))**c2
          Cdn_atm_pond = min(Cdn_atm_pond,camax)
        endif

      !------------------------------------------------------------
      ! Floe edge drag effect (ocean)
      !------------------------------------------------------------

        if (hdraft > puny) then
          scw = c1 - exp(-sl*beta*(c1-ai))
          ctecwf = cfw*p5*(log(hdraft*ocnrufi)/log(zref*ocnrufi))**c2*scw
          Cdn_ocn_floe = ctecwf * hdraft / lfloe
          Cdn_ocn_floe = max(min(Cdn_ocn_floe,cwmax),c0)
        endif

      !------------------------------------------------------------
      ! Total drag coefficient (atmo)
      !------------------------------------------------------------

         Cdn_atm = Cdn_atm_skin + Cdn_atm_floe + Cdn_atm_pond + Cdn_atm_rdg
         Cdn_atm = min(Cdn_atm,camax)

      !------------------------------------------------------------
      ! Total drag coefficient (ocean)
      !------------------------------------------------------------

         Cdn_ocn = Cdn_ocn_skin + Cdn_ocn_floe + Cdn_ocn_keel
         Cdn_ocn = min(Cdn_ocn,cwmax)

      endif

      end subroutine neutral_drag_coeffs

!=======================================================================
!autodocument_start icepack_atm_boundary
!

      subroutine icepack_atm_boundary(sfctype,                   &
                                     Tsf,         potT,          &
                                     uatm,        vatm,          &
                                     wind,        zlvl,          &
                                     Qa,          rhoa,          &
                                     strx,        stry,          &
                                     Tref,        Qref,          &
                                     delt,        delq,          &
                                     lhcoef,      shcoef,        &
                                     Cdn_atm,                    &
                                     Cdn_atm_ratio_n,            &
                                     Qa_iso,      Qref_iso,      &
                                     uvel,        vvel,          &
                                     Uref,        zlvs)

      character (len=3), intent(in) :: &
         sfctype      ! ice or ocean

      real (kind=dbl_kind), intent(in) :: &
         Tsf      , & ! surface temperature of ice or ocean
         potT     , & ! air potential temperature  (K)
         uatm     , & ! x-direction wind speed (m/s)
         vatm     , & ! y-direction wind speed (m/s)
         wind     , & ! wind speed (m/s)
         zlvl     , & ! atm level height for momentum (and scalars if zlvs is not present) (m)
         Qa       , & ! specific humidity (kg/kg)
         rhoa         ! air density (kg/m^3)

      real (kind=dbl_kind), intent(inout) :: &
         Cdn_atm  , &    ! neutral drag coefficient
         Cdn_atm_ratio_n ! ratio drag coeff / neutral drag coeff

      real (kind=dbl_kind), intent(inout) :: &
         strx     , & ! x surface stress (N)
         stry         ! y surface stress (N)

      real (kind=dbl_kind), intent(inout) :: &
         Tref     , & ! reference height temperature  (K)
         Qref     , & ! reference height specific humidity (kg/kg)
         delt     , & ! potential T difference   (K)
         delq     , & ! humidity difference      (kg/kg)
         shcoef   , & ! transfer coefficient for sensible heat
         lhcoef       ! transfer coefficient for latent heat

      real (kind=dbl_kind), intent(in), dimension(:), optional :: &
         Qa_iso       ! specific isotopic humidity (kg/kg)

      real (kind=dbl_kind), intent(inout), dimension(:), optional :: &
         Qref_iso     ! reference specific isotopic humidity (kg/kg)

      real (kind=dbl_kind), intent(in), optional :: &
         uvel     , & ! x-direction ice speed (m/s)
         vvel     , & ! y-direction ice speed (m/s)
         zlvs         ! atm level height for scalars (if different than zlvl) (m)

      real (kind=dbl_kind), intent(out), optional :: &
         Uref         ! reference height wind speed (m/s)

!autodocument_end

      ! local variables

      real (kind=dbl_kind) :: &
         l_uvel, l_vvel, l_Uref

      logical (kind=log_kind), save :: &
         first_call_ice = .true.   ! first call flag

      character(len=*),parameter :: subname='(icepack_atm_boundary)'

      !------------------------------------------------------------
      ! Check optional arguments
      ! Need separate first_call flags for 'ice' and 'ocn' sfctype
      !------------------------------------------------------------


      l_uvel = c0
      l_vvel = c0
      l_Uref = c0
      if (present(uvel)) then
         l_uvel = uvel
      endif
      if (present(vvel)) then
         l_vvel = vvel
      endif

      Cdn_atm_ratio_n = c1

      if (trim(atmbndy) == 'constant') then
         call atmo_boundary_const (sfctype,  calc_strair, &
                                   uatm,     vatm,     &
                                   wind,     rhoa,     &
                                   strx,     stry,     &
                                   Tsf,      potT,     &
                                   Qa,                 &
                                   delt,     delq,     &
                                   lhcoef,   shcoef    )
      else
         call atmo_boundary_layer (sfctype,                 &
                                   calc_strair, formdrag,   &
                                   Tsf,      potT,          &
                                   uatm,     vatm,          &
                                   wind,     zlvl,          &
                                   Qa,       rhoa,          &
                                   strx,     stry,          &
                                   Tref,     Qref,          &
                                   delt,     delq,          &
                                   lhcoef,   shcoef,        &
                                   Cdn_atm,                 &
                                   Cdn_atm_ratio_n,         &
                                   Qa_iso=Qa_iso,           &
                                   Qref_iso=Qref_iso,       &
                                   uvel=l_uvel, vvel=l_vvel,&
                                   Uref=l_Uref, zlvs=zlvs   )
      endif ! atmbndy

      if (present(Uref)) then
         Uref = l_Uref
      endif

      end subroutine icepack_atm_boundary

!=======================================================================

      function compute_stability_parameter(zlvl , thva , &
                                           ustar, tstar, &
                                           qstar, Qa)    &
                                           result(hol)

      real (kind=dbl_kind), intent(in) :: &
         zlvl     , & ! atm level height (m)
         thva     , & ! virtual temperature      (K)
         ustar    , & ! turbulent scale for momentum
         tstar    , & ! turbulent scale for temperature
         qstar    , & ! turbulent scale for humidity
         Qa           ! specific humidity (kg/kg)

      real (kind=dbl_kind) :: &
         hol          ! H (at zlvl) over L

      character(len=*),parameter :: subname='(compute_stability_parameter)'

      hol = vonkar * gravit * zlvl &
               * (tstar/thva &
               + qstar/(c1/zvir+Qa)) &
               / ustar**2
      hol    = sign( min(abs(hol),c10), hol)

      end function compute_stability_parameter

!=======================================================================

      subroutine compute_stability_function(qty, hol, stable, psi)

      character (len=*), intent(in) :: &
         qty          ! 'momentum' or 'scalar'

      real (kind=dbl_kind), intent(in) :: &
         hol          ! H over L

      real (kind=dbl_kind), intent(out) :: &
         psi          , & ! stability function at hol
         stable           ! unit step function at hol

      ! local variables

      real (kind=dbl_kind) :: &
         psi_stable   , & ! stable stability funcion at hol
         psi_unstable     ! unstable stability funcion at hol

      character(len=*),parameter :: subname='(compute_stability_function)'

      stable = p5 + sign(p5 , hol)

      psi_stable = -(0.7_dbl_kind*hol &
                 + 0.75_dbl_kind*(hol-14.3_dbl_kind) &
                 * exp(-0.35_dbl_kind*hol) + 10.7_dbl_kind)

      if(trim(qty) == 'momentum') then
         psi_unstable = psi_momentum_unstable(hol)
      elseif(trim(qty) == 'scalar') then
         psi_unstable = psi_scalar_unstable(hol)
      else
         call icepack_warnings_add(subname//' incorrect qty: ' // qty)
         call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
      endif

      psi = psi_stable*stable + (c1 - stable)*psi_unstable

   end subroutine compute_stability_function

!------------------------------------------------------------
! Define functions
!------------------------------------------------------------

!=======================================================================

      real(kind=dbl_kind) function psi_momentum_unstable(hol)

      real(kind=dbl_kind), intent(in) :: hol

      real(kind=dbl_kind) :: xd

      xd = capital_X(hol)

      psi_momentum_unstable = log((c1+xd*(c2+xd))*(c1+xd*xd)/c8) &
                              - c2*atan(xd) + pih

      end function psi_momentum_unstable

!=======================================================================

      real(kind=dbl_kind) function psi_scalar_unstable(hol)

      real(kind=dbl_kind), intent(in) :: hol

      real(kind=dbl_kind) :: xd

      xd = capital_X(hol)

      psi_scalar_unstable =  c2 * log((c1 + xd*xd)/c2)

      end function psi_scalar_unstable

      real(kind=dbl_kind) function capital_X(hol)

      real(kind=dbl_kind), intent(in) :: hol

      capital_X =  sqrt(max(sqrt(abs(c1 - c16*hol)) , c1))

      end function capital_X

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: helfsurface
! !INTERFACE:
       SUBROUTINE helfsurface(VUS,VVS,VT1,VT2,VSH1,VSH2,VP,VPE, &
        VZ0,LAI,IVWATER,VHS,N,IRUN, &
        VRHO,VKH,VKM,VUSTAR,VXX,VYY,VCU,VCT,VRIB,VZETA,VWS, &
        t2m,q2m,u2m,v2m,t10m,q10m,u10m,v10m,u50m,v50m,CHOOSEZ0,WMCHARNOCK)
!**********************************************************************
!  SUBROUTINE helfsurface - COMPUTES SURFACE TRANSFER COEFFICIENTS
!
!   ARGUMENTS ::
!
!     INPUT:
!     ------
!    US            -         U - COMPONENT OF SURFACE WIND
!    VS            -         V - COMPONENT OF SURFACE WIND
!    THV1          -         VIRTUAL POTENTIAL TEMPERATURE AT NLAY
!    THV2          -         VIRTUAL POTENTIAL TEMPERATURE AT GROUND
!    TH1           -         POTENTIAL TEMPERATURE AT NLAY
!    TH2           -         POTENTIAL TEMPERATURE AT GROUND
!    SH1           -         SPECIFIC HUMIDITY AT NLAY
!    SH2           -         SPECIFIC HUMIDITY AT GROUND
!    PK            -         EVEN LEVEL PRESSURE ** KAPPA AT LEVEL NLAY
!    PKE           -         EDGE LEVEL PRESSURE ** KAPPA AT GROUND
!    PE            -         SURFACE PRESSURE
!    Z0            -         SURFACE ROUGHNESS
!    WATER         -         ARRAY WITH '1' OVER OCEANS
!    HS            -         DEPTH OF SURFACE LAYER
!    N             -         NUMBER OF helfsurface ITERATIONS
!    CHOOSEZ0      -         INTEGER FLAG: 0 - L&P Z0, no high wind limit
!                                          1 - Edson Z0 for mom. and heat, high wind limit
!                                          2 - L&P Z0, high wind limit
!                                          3 - Edson Z0 for mom. only, high wind limit
!                                          4 - wave model Charnock coefficient
!     OUTPUT:
!     -------
!    RHO           -         DENSITY AT SURFACE
!    KH            -         HEAT TRANSFER COEFFICIENT (CT*USTAR)
!    KM            -         MOMENTUM TRANSFER COEFFICIENT (CU*USTAR)
!    USTAR         -         FRICTION VELOCITY
!    XX            -         PHIM(ZETA) - DIMENSIONLESS WIND SHEAR
!    YY            -         PHIH(ZETA) - DIMENSIONLESS TEMP GRADIENT
!    CU            -         MOMENTUM TRANSPORT COEFFICIENT
!    CT            -         HEAT TRANSPORT COEFFICIENT
!
!**********************************************************************
      implicit none

! Argument List Declarations
      integer n,irun,CHOOSEZ0
      real VUS(:),VVS(:),VT1(:),VT2(:),VSH1(:),VSH2(:)
      real VPE(:),VP(:),VZ0(:),LAI(:),VHS(:)
      integer IVWATER(:)
      real VRHO(:)
      real VKM(:),VKH(:),VUSTAR(:),VXX(:)
      real VYY(:),VCU(:),VCT(:),VRIB(:)
      real VZETA(:),VWS(:)
      real, intent(OUT) :: t2m(:),q2m(:),u2m(:),v2m(:)
      real, intent(OUT) :: t10m(:),q10m(:),u10m(:),v10m(:)
      real, intent(OUT) :: u50m(:),v50m(:)
      LOGICAL LWATER
      integer IVBITRIB(irun)
      real, optional, intent(in) :: WMCHARNOCK(:)

! Local Variables
      real VHZ(irun),VPSIM(irun),VAPSIM(irun),VPSIG(irun),VPSIHG(irun)
      real VTEMP(irun),VDZETA(irun),VDZ0(irun),VDPSIM(irun)
      real VDPSIH(irun),VZH(irun),VXX0(irun),VYY0(irun)
      real VAPSIHG(irun),VRIB1(irun)
      real VPSIH(irun),VPSIH2(irun),VH0(irun)
      real VX0PSIM(irun),VG(irun),VG0(irun),VR1MG0(irun)
      real VZ2(irun),VDZSEA(irun),VAZ0(irun),VXNUM1(irun)
      real VPSIGB2(irun),VDX(irun),VDXPSIM(irun),VDY(irun)
      real VXNUM2(irun),VDEN(irun),VAWS1(irun),VXNUM3(irun)
      real VXNUM(irun),VDZETA1(irun),VDZETA2(irun)
      real VZCOEF2(irun),VZCOEF1(irun),VTEMPLIN(irun)
      real VDPSIMC(irun),VDPSIHC(irun),VAHS(irun)
      real VTHV1(IRUN),VTHV2(IRUN),VTH1(IRUN),VTH2(IRUN),VPKE(IRUN),VPK(IRUN)

      real vz0h(irun),vh0h(irun),dummy1(irun),dummy2(irun),dummy3(irun),dummy4(irun),dummy5(irun)

! Local Variables
      real USTMX3,USTZ0S,Z0MIN,H0BYZ0,USTH0S,H0VEG,Z0VEGM,PRFAC,Z0MAX
      real XPFAC,DIFSQT
      PARAMETER ( USTMX3 =   0.0632456)
      PARAMETER ( USTZ0S =   0.2030325E-5)
      PARAMETER ( Z0MIN  =  USTZ0S/USTMX3)
      PARAMETER ( Z0MAX  =  USTZ0S/USTMX3)
      PARAMETER ( H0BYZ0 =    30.0    )
      PARAMETER ( USTH0S =  H0BYZ0*USTZ0S )
      PARAMETER ( Z0VEGM =   0.005    )
      PARAMETER ( H0VEG  =  H0BYZ0*Z0VEGM )  !! This prevents discontinuity
      PARAMETER ( PRFAC  = 0.595864   )
      PARAMETER ( XPFAC  = .55        )  
      PARAMETER ( DIFSQT  = 3.872983E-3)

      real psihdiag(irun),psimdiag(irun)
      real rvk,vk2,bmdl(irun)
      integer iwater,itype
      integer i,iter

      real VCH(irun)

!
      if (present(WMCHARNOCK)) then
         VCH = WMCHARNOCK
      else
         VCH = 0.018
      end if

!
      UNUSED_DUMMY(LAI)
      rvk = 1./MAPL_KARMAN
      vk2 = MAPL_KARMAN*MAPL_KARMAN
      DO I = 1,IRUN
      if( ivwater(i) .eq. 3 ) then 
       BMDL(i)    = 0.
!scale BMDL(i)    = (MAPL_KARMAN * XPFAC * PRFAC / DIFSQT) * exp(-lai(i)*2.)
      else
       BMDL(i)    = (MAPL_KARMAN * XPFAC * PRFAC / DIFSQT)
      endif
      enddo

!     INITIALIZATION 

      DO I = 1,IRUN
       VAHS(I) = 1. / VHS(I)
       VPKE(I) = VPE(I) ** MAPL_KAPPA
       VPK(I) = VP(I) ** MAPL_KAPPA
       VTH1(I) = VT1(I)/VPK(I)
       VTH2(I) = VT2(I)/VPKE(I)
       VTHV1(I) = VTH1(I)*( 1.0 + MAPL_VIREPS*VSH1(I))
       VTHV2(I) = VTH2(I)*( 1.0 + MAPL_VIREPS*VSH2(I))
      ENDDO

!     DETERMINE SURFACE WIND MAGNITUDE AND BULK RICHARDSON NUMBER
!
      DO I = 1,IRUN
       VWS(I) = max(VUS(I)*VUS(I) + VVS(I)*VVS(I),1.e-4)
       VRIB(I) = MAPL_CP*(VPKE(I)-VPK(I))*(VTHV1(I)-VTHV2(I)) / VWS(I)
       VWS(I) = SQRT( VWS(I) )
      ENDDO    

!  INITIAL GUESS FOR ROUGHNESS LENGTH Z0 OVER WATER
!
      IWATER = 0
      DO 9002 I = 1,IRUN
       IF (IVWATER(I).EQ.1)  IWATER = IWATER + 1
 9002 CONTINUE
      LWATER = .FALSE.
      IF(IWATER.GE.1)LWATER = .TRUE.
!
      IF(LWATER)THEN
       DO 9004 I = 1,IRUN
        IF (IVWATER(I).EQ.1) VZ0(I) = 0.0003
 9004  CONTINUE
      ENDIF
      do i = 1,irun
       vh0(i) = h0byz0 * vz0(i)
       if(vz0(i).ge.z0vegm)vh0(i) = h0veg
      enddo
       DO I = 1,IRUN
        VZ0H(I) = 0.001
       ENDDO

!     CU AND PSIHG FOR NEUTRALLY STRATIFIED FLOW
!
      DO 9006 I = 1,IRUN
       VHZ(I) = (VHS(I) / VZ0(I) + 1.)
       VPSIM(I) = LOG( VHZ(I) )
       VAPSIM(I) = 1. / VPSIM(I)
       VCU(I) = MAPL_KARMAN * VAPSIM(I)
       VUSTAR(I) = VCU(I) * VWS(I)
!
       VPSIG(I) = BMDL(i)*sqrt(max(VH0(I)*VUSTAR(I)-USTH0S,0.))
       VPSIHG(I) = VPSIM(I) + VPSIG(I)
 9006 CONTINUE

!
!     LINEAR CORRECTION FOR ERROR IN ROUGHNESS LENGTH Z0
!
      IF(LWATER)THEN
       DO 9008 I = 1,IRUN
        VTEMP(I) = 0.
 9008  CONTINUE
       CALL LINADJ(VRIB,VRIB,VWS,VWS,VZ0,VUSTAR,IVWATER,VAPSIM, &
        VTEMP,VTEMP,VTEMP,VTEMP,VTEMP,VTEMP,VTEMP,1,.TRUE.,IRUN,VDZETA, &
        VDZ0,VDPSIM,VDPSIH,IVBITRIB, &
        VX0PSIM,VG,VG0,VR1MG0,VZ2,VDZSEA,VAZ0,VXNUM1,VPSIGB2,VDX, &
        VDXPSIM,VDY,VXNUM2,VDEN,VAWS1,VXNUM3,VXNUM,VDZETA1,VDZETA2, &
        VZCOEF2,VZCOEF1,VTEMPLIN,VDPSIMC,VDPSIHC,MAPL_KARMAN,bmdl,CHOOSEZ0,VCH)
       DO 9010 I = 1,IRUN
        IF ( IVWATER(I).EQ.1 ) THEN
         VCU(I) = VCU(I) * (1. - VDPSIM(I)*VAPSIM(I))
         VZ0(I) = VZ0(I) + VDZ0(I)
         ENDIF 
         IF ( IVWATER(I).EQ.1) THEN
         IF ( VZ0(I) .LE. Z0MIN ) VZ0(I) = Z0MIN 
         vh0(i) = h0byz0 * vz0(i)
         VPSIG(I) = VH0(I) * VCU(I) * VWS(I) - USTH0S
         if(VPSIG(I).lt.0.)  VPSIG(I) = 0.
         VPSIG(I) = SQRT( VPSIG(I) )
         VPSIG(I) = BMDL(i) * VPSIG(I)
         VPSIHG(I) = VPSIM(I) + VDPSIH(I) + VPSIG(I)
        ENDIF  
 9010  CONTINUE
!
      ENDIF
!
!  INITIAL GUESS FOR STABILITY PARAMETER ZETA
!
      DO 9012 I = 1,IRUN
       VZETA(I) = VK2 * VRIB(I) / (VCU(I) * VCU(I) * VPSIHG(I))
 9012 CONTINUE
!
!  RECOMPUTE CU, ESTIMATE PSIHG AND UPDATE ZETA AND Z0
!
      DO 9014 I = 1,IRUN
!      VZH(I) = VZ0(I) * VAHS(I)
       VZH(I) = VZ0(I) / (VHS(I) + VZ0(I))
 9014 CONTINUE
      CALL PSI (VZETA,VZH,VPSIM,VTEMP,IRUN,VXX,VXX0,VYY,VYY0,2)
      DO 9016 I = 1,IRUN
       VCU(I) = MAPL_KARMAN / VPSIM(I)
       VPSIG(I) = VH0(I) * VCU(I) * VWS(I) - USTH0S
       if(VPSIG(I).lt.0.)  VPSIG(I) = 0.
       VPSIG(I) = SQRT(VPSIG(I))
       VPSIG(I) = BMDL(i) * VPSIG(I)
       VPSIHG(I) = VPSIM(I) + VPSIG(I)
       VZETA(I) = VK2 * VRIB(I) / (VCU(I) * VCU(I) * VPSIHG(I))
 9016 CONTINUE
!
      IF(LWATER)THEN
       DO 9018 I = 1,IRUN
        IF (IVWATER(I).EQ.1) VUSTAR(I) = VCU(I) * VWS(I)
 9018  CONTINUE
       CALL ZCSUB ( VUSTAR,VCH,VHZ,IVWATER,.FALSE.,IRUN,VTEMP,CHOOSEZ0)
       CALL ZCSUB ( VUSTAR,VCH,VHZ,IVWATER,.FALSE.,IRUN,vz0h,2)
       DO 9020 I = 1,IRUN
        IF (IVWATER(I).EQ.1 ) then
         VZ0(I) = VTEMP(I)
         IF ( VZ0(I) .LE. Z0MIN ) VZ0(I) = Z0MIN
         IF ( VZ0H(I) .LE. Z0MIN ) VZ0H(I) = Z0MIN
         vh0(i) = h0byz0 * vz0(i)
         vh0h(i) = h0byz0 * vz0h(i)
        endif
 9020  CONTINUE
      ENDIF
!
!  ITERATIVE LOOP - N ITERATIONS
!     COMPUTE CU AND CT
!
      DO 200 ITER = 1,N

       DO 9026 I = 1,IRUN
!       VZH(I) = VZ0(I) * VAHS(I)
        VZH(I) = VZ0(I) / (VHS(I) + VZ0(I))
 9026  CONTINUE
       CALL PSI (VZETA,VZH,VPSIM,VPSIH,IRUN,VXX,VXX0,VYY,VYY0,1)
       DO I = 1,IRUN
!       VZH(I) = VZ0H(I) * VAHS(I)
        VZH(I) = VZ0H(I) / (VHS(I) + VZ0H(I))
       ENDDO
       if( choosez0.eq.3 .AND. Lwater ) CALL PSI (VZETA,VZH,dummy1,VPSIH,IRUN,dummy2,dummy3,dummy4,dummy5,3)
       DO 9028 I = 1,IRUN
        VCU(I) = MAPL_KARMAN / VPSIM(I)
        VUSTAR(I) = VCU(I) * VWS(I)
!
        VPSIG(I) = VH0(I) * VUSTAR(I) - USTH0S
        if(VPSIG(I).lt.0.)  VPSIG(I) = 0.
        VPSIG(I) = SQRT(VPSIG(I))
        VPSIG(I) = BMDL(i) * VPSIG(I)
        VPSIHG(I) = VPSIH(I) + VPSIG(I)
!
!  LINEAR CORRECTIONS FOR CU, CT, ZETA, AND Z0
!
        VAPSIM(I) = VCU(I) * RVK
        VAPSIHG(I) = 1. / VPSIHG(I)
        VRIB1(I) = VAPSIM(I) * VAPSIM(I) * VPSIHG(I) * VZETA(I)
 9028  CONTINUE
!
       ITYPE = 3
       IF(ITER.EQ.N) ITYPE = 5
!
       CALL LINADJ(VRIB1,VRIB,VWS, &
        VWS,VZ0,VUSTAR,IVWATER, &
        VAPSIM,VAPSIHG,VPSIH, &
        VPSIG,VXX,VXX0, &
        VYY,VYY0,ITYPE,LWATER,IRUN,VDZETA, &
        VDZ0,VDPSIM,VDPSIH, &
        IVBITRIB, &
       VX0PSIM,VG,VG0,VR1MG0,VZ2,VDZSEA,VAZ0,VXNUM1,VPSIGB2,VDX, &
       VDXPSIM,VDY,VXNUM2,VDEN,VAWS1,VXNUM3,VXNUM,VDZETA1,VDZETA2, &
       VZCOEF2,VZCOEF1,VTEMPLIN,VDPSIMC,VDPSIHC,MAPL_KARMAN,bmdl,CHOOSEZ0,VCH)
!
!  UPDATES OF ZETA, Z0, CU AND CT
!
       DO 9032 I = 1,IRUN
        VZETA(I) = VZETA(I) * ( 1. + VDZETA(I) )
        IF (IVBITRIB(I).EQ.1 ) VZETA(I) = VPSIM(I) * VPSIM(I) * VRIB(I) * VAPSIHG(I)
 9032  CONTINUE
!
       IF ( LWATER ) THEN
        DO 9034 I = 1,IRUN
         IF (IVWATER(I).EQ.1 ) then
          VZ0(I) = VZ0(I) * ( 1. + VDZ0(I) )
          VZ0H(I) = VZ0H(I) * ( 1. + VDZ0(I) )
          IF (VZ0(I) .LE. Z0MIN ) VZ0(I) = Z0MIN
          IF (VZ0H(I) .LE. Z0MIN ) VZ0H(I) = Z0MIN
          vh0(i) = h0byz0 * vz0(i)
          vh0h(i) = h0byz0 * vz0h(i)
         endif
 9034   CONTINUE
       ENDIF
!
       IF ( ITER .EQ. N ) THEN
        DO 9036 I = 1,IRUN
         VPSIM(I) = VPSIM(I) + VDPSIM(I)
         VCU(I) = MAPL_KARMAN / VPSIM(I)
         VUSTAR(I) = VCU(I) * VWS(I)
!
         VPSIG(I) = VH0(I) * VUSTAR(I) - USTH0S
         if(VPSIG(I).lt.0.)  VPSIG(I) = 0.
         VPSIG(I) = SQRT(VPSIG(I))
         VPSIG(I) = BMDL(i) * VPSIG(I)
         VPSIHG(I) = VPSIH(I) + VDPSIH(I) + VPSIG(I)
         VCT(I) = MAPL_KARMAN / VPSIHG(I)
 9036   CONTINUE
       ENDIF

!
!  SAVE VALUES OF RIB AND WS
!
        DO 9038 I = 1,IRUN
         VRIB1(I) = VRIB(I)
 9038   CONTINUE
!
 200  CONTINUE
!
!  CALCULATE RHO-SURFACE ( KG / M**3 )
!
       DO I = 1,IRUN
        VTEMP(I) =  10. * VAHS(I) * VZETA(I)
!       VZH(I) = VZ0(I) * 0.1
        VZH(I) = VZ0(I) / (10. + VZ0(I))
       ENDDO
       CALL PSI (VTEMP,VZH,VHZ,VPSIH2,IRUN,VHZ,VHZ,VHZ,VHZ,3)
       DO I = 1,IRUN
        VTEMP(I) = min(( VPSIH2(I) + VPSIG(I) ) / VPSIHG(I),1.)
        VRHO(I) = VPKE(I)*( VTH2(I) + VTEMP(I) * (VTH1(I)-VTH2(I)) )
        VRHO(I) = VPE(I)*100. / ( MAPL_RGAS * VRHO(I) )
       ENDDO
!
! interpolate uvtq to 2, 10 and 50 meters for diagnostic output
!  use psih and psim which represent non-dim change from ground
!                 to specified level
! and multiply theta by surface p**kappa to get temperatures
!
        do i = 1,irun
         vtemp(i) = 2. * vahs(i) * vzeta(i)
!        vzh(i) = min(vz0(i),2.) * 0.5
         VZH(I) = min(VZ0(I),2.) / (2. + min(VZ0(I),2.))
        enddo
        call psi(vtemp,vzh,psimdiag,psihdiag,irun,vhz,vhz,vhz,vhz,1)
        do i = 1,irun
         vtemp(i) = min(( psihdiag(i) + vpsig(i) ) / vpsihg(i),1.)
         t2m(i) = ( (vth2(i) + vtemp(i)* (vth1(i)-vth2(i))) ) * vpke(i)
         q2m(i) = (vsh2(i) + vtemp(i)* (vsh1(i)-vsh2(i)))
         u2m(i) = (psimdiag(i)/vpsim(i) * vus(i))
         v2m(i) = (psimdiag(i)/vpsim(i) * vvs(i))
        enddo

        do i = 1,irun
         vtemp(i) = 10. * vahs(i) * vzeta(i)
!        vzh(i) = vz0(i) * 0.1
         VZH(I) = VZ0(I) / (10. + VZ0(I))
        enddo
        call psi(vtemp,vzh,psimdiag,psihdiag,irun,vhz,vhz,vhz,vhz,1)
        do i = 1,irun
         vtemp(i) = min(( psihdiag(i) + vpsig(i) ) / vpsihg(i),1.)
         t10m(i) = ( (vth2(i) + vtemp(i)* (vth1(i)-vth2(i))) ) * vpke(i)
         q10m(i) = (vsh2(i) + vtemp(i)* (vsh1(i)-vsh2(i)))
         u10m(i) = (psimdiag(i)/vpsim(i) * vus(i))
         v10m(i) = (psimdiag(i)/vpsim(i) * vvs(i))
        enddo

        do i = 1,irun
         vtemp(i) = 50. * vahs(i) * vzeta(i)
!        vzh(i) = vz0(i) * 0.02
         VZH(I) = VZ0(I) / (50. + VZ0(I))
        enddo
        call psi(vtemp,vzh,psimdiag,psihdiag,irun,vhz,vhz,vhz,vhz,1)
        do i = 1,irun
         u50m(i) = (psimdiag(i)/vpsim(i) * vus(i))
         v50m(i) = (psimdiag(i)/vpsim(i) * vvs(i))
        enddo
!
!  EVALUATE TURBULENT TRANSFER COEFFICIENTS
!

      DO 9044 I = 1,IRUN
!!     VKH(I) = VUSTAR(I) * VCT(I)
!!     VKM(I) = VUSTAR(I) * VCU(I)
       VKH(I) = VUSTAR(I) * VCT(I) * VRHO(I)
       VKM(I) = VUSTAR(I) * VCU(I) * VRHO(I)
 9044 CONTINUE

      DO I = 1,IRUN
       VRIB(I) = MAPL_CP*(VPKE(I)-VPK(I))*(VTHV1(I)-VTHV2(I)) /    &
                max(VUS(I)*VUS(I) + VVS(I)*VVS(I),1.e-1)
      ENDDO    

end subroutine helfsurface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: phi
! !INTERFACE:
      SUBROUTINE PHI(Z,PHIM,PHIH,IFLAG,N)
!**********************************************************************
!
!  FUNCTION PHI - SOLVES KEYPS EQUATIONS
!               - CALLED FROM PSI
!
!  DESCRIPTION OF PARAMETERS
!     Z     -  INPUTED VALUE OF MONIN- OBUKHOV STABILITY PARAMETER ZETA
!               TIMES APPROPRIATE CONSTANT
!     PHIM  -  OUTPUTED SOLUTION OF KEYPS EQUATION FOR MOMENTUM
!     PHIH  -  OUTPUTED SOLUTION OF KEYPS EQUATION FOR SCALARS
!     IFLAG -  FLAG TO DETERMINE IF X IS NEEDED (IFLAG=2), Y IS NEEDED
!                  (IFLAG=3), OR BOTH (IFLAG=1)
!     N     -  LENGTH OF VECTOR TO BE SOLVED
!
!**********************************************************************
      implicit none

! Argument List Declarations
      integer n,iflag
      real PHIM(:),PHIH(:),Z(:)

! Local Variables
      integer I1(N),I2(N)
      real ZSTAR(N),E1(N),E2(N),TEMP1(N)
!
      real PHIM0(385),ZLINM1(75),ZLINM2(75),ZLINM3(36)
      real ZLOGM1(74),ZLOGM2(75),ZLOGM3(50)
      real PHIH0(385),ZLINH1(75),ZLINH2(75),ZLINH3(36)
      real ZLOGH1(74),ZLOGH2(75),ZLOGH3(50)
      EQUIVALENCE (PHIM0(1),ZLINM1(1)),(PHIM0(76),ZLINM2(1))
      EQUIVALENCE (PHIM0(151),ZLINM3(1))
      EQUIVALENCE (PHIM0(187),ZLOGM1(1)),(PHIM0(261),ZLOGM2(1))
      EQUIVALENCE (PHIM0(336),ZLOGM3(1))
      EQUIVALENCE (PHIH0(1),ZLINH1(1)),(PHIH0(76),ZLINH2(1))
      EQUIVALENCE (PHIH0(151),ZLINH3(1))
      EQUIVALENCE (PHIH0(187),ZLOGH1(1)),(PHIH0(261),ZLOGH2(1))
      EQUIVALENCE (PHIH0(336),ZLOGH3(1))
!
       DATA ZLOGM1/ &
                   0.697894,0.678839,0.659598,0.640260, &
        0.620910,0.601628,0.582486,0.563550,0.544877, &
        0.526519,0.508516,0.490903,0.473708,0.456951, &
        0.440649,0.424812,0.409446,0.394553,0.380133, &
        0.366182,0.352695,0.339664,0.327082,0.314938, &
        0.303222,0.291923,0.281029,0.270528,0.260409, &
        0.250659,0.241267,0.232221,0.223509,0.215119, &
        0.207041,0.199264,0.191776,0.184568,0.177628, &
        0.170949,0.164519,0.158331,0.152374,0.146641, &
        0.141123,0.135813,0.130702,0.125783,0.121048, &
        0.116492,0.112107,0.107887,0.103826,0.0999177, &
        0.0961563,0.0925364,0.0890528,0.0857003,0.0824739, &
        0.0793690,0.0763810,0.0735054,0.0707380,0.0680749, &
        0.0655120,0.0630455,0.0606720,0.0583877,0.0561895, &
        0.0540740,0.0520382,0.0500790,0.0481936,0.0463791/
       DATA ZLOGM2/ &
        0.0446330,0.0429526,0.0413355,0.0397792,0.0382816, &
        0.0368403,0.0354533,0.0341185,0.0328340,0.0315978, &
        0.0304081,0.0292633,0.0281616,0.0271013,0.0260809, &
        0.0250990,0.0241540,0.0232447,0.0223695,0.0215273, &
        0.0207168,0.0199369,0.0191862,0.0184639,0.0177687, &
        0.0170998,0.0164560,0.0158364,0.0152402,0.0146664, &
        0.0141142,0.0135828,0.0130714,0.0125793,0.0121057, &
        0.0116499,0.0112113,0.0107892,0.0103830,0.999210E-2, &
        0.961590E-2,0.925387E-2,0.890547E-2,0.857018E-2,0.824752E-2, &
        0.793701E-2,0.763818E-2,0.735061E-2,0.707386E-2,0.680754E-2, &
        0.655124E-2,0.630459E-2,0.606722E-2,0.583880E-2,0.561897E-2, &
        0.540742E-2,0.520383E-2,0.500791E-2,0.481937E-2,0.463792E-2, &
        0.446331E-2,0.429527E-2,0.413355E-2,0.397793E-2,0.382816E-2, &
        0.368403E-2,0.354533E-2,0.341185E-2,0.328340E-2,0.315978E-2, &
        0.304082E-2,0.292633E-2,0.281616E-2,0.271013E-2,0.260809E-2/
       DATA ZLOGM3/ &
        0.250990E-2,0.241541E-2,0.232447E-2,0.223695E-2,0.215273E-2, &
        0.207168E-2,0.199369E-2,0.191862E-2,0.184639E-2,0.177687E-2, &
        0.170998E-2,0.164560E-2,0.158364E-2,0.152402E-2,0.146664E-2, &
        0.141142E-2,0.135828E-2,0.130714E-2,0.125793E-2,0.121057E-2, &
        0.116499E-2,0.112113E-2,0.107892E-2,0.103830E-2,0.999210E-3, &
        0.961590E-3,0.925387E-3,0.890547E-3,0.857018E-3,0.824752E-3, &
        0.793701E-3,0.763818E-3,0.735061E-3,0.707386E-3,0.680754E-3, &
        0.655124E-3,0.630459E-3,0.606722E-3,0.583880E-3,0.561897E-3, &
        0.540742E-3,0.520383E-3,0.500791E-3,0.481937E-3,0.463792E-3, &
        0.446331E-3,0.429527E-3,0.413355E-3,0.397793E-3,0.382816E-3/
       DATA ZLOGH1/ &
                   0.640529,0.623728,0.606937,0.590199, &
        0.573552,0.557032,0.540672,0.524504,0.508553, &
        0.492843,0.477397,0.462232,0.447365,0.432809, &
        0.418574,0.404670,0.391103,0.377878,0.364999, &
        0.352468,0.340284,0.328447,0.316954,0.305804, &
        0.294992,0.284514,0.274364,0.264538,0.255028, &
        0.245829,0.236933,0.228335,0.220026,0.211999, &
        0.204247,0.196762,0.189537,0.182564,0.175837, &
        0.169347,0.163088,0.157051,0.151231,0.145620, &
        0.140211,0.134998,0.129974,0.125133,0.120469, &
        0.115975,0.111645,0.107475,0.103458,0.995895E-1, &
        0.958635E-1,0.922753E-1,0.888199E-1,0.854925E-1,0.822886E-1, &
        0.792037E-1,0.762336E-1,0.733739E-1,0.706208E-1,0.679704E-1, &
        0.654188E-1,0.629625E-1,0.605979E-1,0.583217E-1,0.561306E-1, &
        0.540215E-1,0.519914E-1,0.500373E-1,0.481564E-1,0.463460E-1/
       DATA ZLOGH2/ &
        0.446034E-1,0.429263E-1,0.413120E-1,0.397583E-1,0.382629E-1, &
        0.368237E-1,0.354385E-1,0.341053E-1,0.328222E-1,0.315873E-1, &
        0.303988E-1,0.292550E-1,0.281541E-1,0.270947E-1,0.260750E-1, &
        0.250937E-1,0.241494E-1,0.232405E-1,0.223658E-1,0.215240E-1, &
        0.207139E-1,0.199342E-1,0.191839E-1,0.184618E-1,0.177669E-1, &
        0.170981E-1,0.164545E-1,0.158351E-1,0.152390E-1,0.146653E-1, &
        0.141133E-1,0.135820E-1,0.130707E-1,0.125786E-1,0.121051E-1, &
        0.116494E-1,0.112108E-1,0.107888E-1,0.103826E-1,0.999177E-2, &
        0.961561E-2,0.925360E-2,0.890523E-2,0.856997E-2,0.824733E-2, &
        0.793684E-2,0.763803E-2,0.735048E-2,0.707375E-2,0.680743E-2, &
        0.655114E-2,0.630450E-2,0.606715E-2,0.583873E-2,0.561891E-2, &
        0.540737E-2,0.520379E-2,0.500787E-2,0.481933E-2,0.463789E-2, &
        0.446328E-2,0.429524E-2,0.413353E-2,0.397790E-2,0.382814E-2, &
        0.368401E-2,0.354532E-2,0.341184E-2,0.328338E-2,0.315977E-2, &
        0.304081E-2,0.292632E-2,0.281615E-2,0.271012E-2,0.260809E-2/
       DATA ZLOGH3/ &
        0.250990E-2,0.241540E-2,0.232446E-2,0.223695E-2,0.215273E-2, &
        0.207168E-2,0.199368E-2,0.191862E-2,0.184639E-2,0.177687E-2, &
        0.170997E-2,0.164559E-2,0.158364E-2,0.152402E-2,0.146664E-2, &
        0.141142E-2,0.135828E-2,0.130714E-2,0.125793E-2,0.121057E-2, &
        0.116499E-2,0.112113E-2,0.107892E-2,0.103830E-2,0.999209E-3, &
        0.961590E-3,0.925387E-3,0.890546E-3,0.857018E-3,0.824752E-3, &
        0.793700E-3,0.763818E-3,0.735061E-3,0.707386E-3,0.680754E-3, &
        0.655124E-3,0.630459E-3,0.606722E-3,0.583880E-3,0.561897E-3, &
        0.540742E-3,0.520383E-3,0.500791E-3,0.481937E-3,0.463792E-3, &
        0.446331E-3,0.429527E-3,0.413355E-3,0.397793E-3,0.382816E-3/
 
       DATA ZLINM1/ &
        0.964508,0.962277,0.960062,0.957863,0.955680, &
        0.953512,0.951359,0.949222,0.947100,0.944992, &
        0.942899,0.940821,0.938758,0.936709,0.934673, &
        0.932652,0.930645,0.928652,0.926672,0.924706, &
        0.922753,0.920813,0.918886,0.916973,0.915072, &
        0.913184,0.911308,0.909445,0.907594,0.905756, &
        0.903930,0.902115,0.900313,0.898522,0.896743, &
        0.894975,0.893219,0.891475,0.889741,0.888019, &
        0.886307,0.884607,0.882917,0.881238,0.879569, &
        0.877911,0.876264,0.874626,0.872999,0.871382, &
        0.869775,0.868178,0.866591,0.865013,0.863445, &
        0.861887,0.860338,0.858798,0.857268,0.855747, &
        0.854235,0.852732,0.851238,0.849753,0.848277, &
        0.846809,0.845350,0.843900,0.842458,0.841025, &
        0.839599,0.838182,0.836774,0.835373,0.833980/
       DATA ZLINM2/ &
        0.832596,0.831219,0.829850,0.828489,0.827136, &
        0.825790,0.824451,0.823121,0.821797,0.820481, &
        0.819173,0.817871,0.816577,0.815289,0.814009, &
        0.812736,0.811470,0.810210,0.808958,0.807712, &
        0.806473,0.805240,0.804015,0.802795,0.801582, &
        0.800376,0.799176,0.797982,0.796794,0.795613, &
        0.794438,0.793269,0.792106,0.790949,0.789798, &
        0.788652,0.787513,0.786380,0.785252,0.784130, &
        0.783014,0.781903,0.780798,0.779698,0.778604, &
        0.777516,0.776432,0.775354,0.774282,0.773215, &
        0.772153,0.771096,0.770044,0.768998,0.767956, &
        0.766920,0.765888,0.764862,0.763840,0.762824, &
        0.761812,0.760805,0.759803,0.758805,0.757813, &
        0.756824,0.755841,0.754862,0.753888,0.752918, &
        0.751953,0.750992,0.750035,0.749083,0.748136/
       DATA ZLINM3/ &
        0.747192,0.746253,0.745318,0.744388,0.743462, &
        0.742539,0.741621,0.740707,0.739798,0.738892, &
        0.737990,0.737092,0.736198,0.735308,0.734423, &
        0.733540,0.732662,0.731788,0.730917,0.730050, &
        0.729187,0.728328,0.727472,0.726620,0.725772, &
        0.724927,0.724086,0.723248,0.722414,0.721584, &
        0.720757,0.719933,0.719113,0.718296,0.717483, &
        0.716673/
       DATA ZLINH1/ &
        0.936397,0.932809,0.929287,0.925827,0.922429, &
        0.919089,0.915806,0.912579,0.909405,0.906284, &
        0.903212,0.900189,0.897214,0.894284,0.891399, &
        0.888558,0.885759,0.883001,0.880283,0.877603, &
        0.874962,0.872357,0.869788,0.867255,0.864755, &
        0.862288,0.859854,0.857452,0.855081,0.852739, &
        0.850427,0.848144,0.845889,0.843662,0.841461, &
        0.839287,0.837138,0.835014,0.832915,0.830841, &
        0.828789,0.826761,0.824755,0.822772,0.820810, &
        0.818869,0.816949,0.815050,0.813170,0.811310, &
        0.809470,0.807648,0.805845,0.804060,0.802293, &
        0.800543,0.798811,0.797095,0.795396,0.793714, &
        0.792047,0.790396,0.788761,0.787141,0.785535, &
        0.783945,0.782369,0.780807,0.779259,0.777724, &
        0.776204,0.774696,0.773202,0.771720,0.770251/
       DATA ZLINH2/ &
        0.768795,0.767351,0.765919,0.764499,0.763091, &
        0.761694,0.760309,0.758935,0.757571,0.756219, &
        0.754878,0.753547,0.752226,0.750916,0.749616, &
        0.748326,0.747045,0.745775,0.744514,0.743262, &
        0.742020,0.740787,0.739563,0.738348,0.737141, &
        0.735944,0.734755,0.733574,0.732402,0.731238, &
        0.730083,0.728935,0.727795,0.726664,0.725539, &
        0.724423,0.723314,0.722213,0.721119,0.720032, &
        0.718952,0.717880,0.716815,0.715756,0.714704, &
        0.713660,0.712621,0.711590,0.710565,0.709547, &
        0.708534,0.707529,0.706529,0.705536,0.704549, &
        0.703567,0.702592,0.701623,0.700660,0.699702, &
        0.698750,0.697804,0.696863,0.695928,0.694998, &
        0.694074,0.693155,0.692241,0.691333,0.690430, &
        0.689532,0.688639,0.687751,0.686868,0.685990/
       DATA ZLINH3/ &
        0.685117,0.684249,0.683386,0.682527,0.681673, &
        0.680824,0.679979,0.679139,0.678303,0.677472, &
        0.676645,0.675823,0.675005,0.674191,0.673381, &
        0.672576,0.671775,0.670978,0.670185,0.669396, &
        0.668611,0.667830,0.667054,0.666281,0.665512, &
        0.664746,0.663985,0.663227,0.662473,0.661723, &
        0.660977,0.660234,0.659495,0.658759,0.658027, &
        0.657298/

        integer i
!
      DO 9002 I = 1,N
       ZSTAR(I)    = 100. * Z(I) - 14.
 9002 CONTINUE
!
      DO 9004 I = 1,N
       TEMP1(I) = Z(I)*0.5
       IF( Z(I) .LE. 2. )TEMP1(I) = 1.
       TEMP1(I) = LOG10(TEMP1(I))
       TEMP1(I) = (TEMP1(I) + 9.3) * 20.
       IF( Z(I) .GT. 2. ) ZSTAR(I) = TEMP1(I)
       IF( Z(I).GT.1.78e10 ) ZSTAR(I) = 384.9999
 9004  CONTINUE
!
 60    CONTINUE
!
      DO 9006 I = 1,N
       I1(I) = ZSTAR(I)
       I2(I) = I1(I) + 1
       TEMP1(I) = ZSTAR(I) - I1(I)
!
 9006  CONTINUE
!
      IF( IFLAG .GT. 2 ) GO TO 100
       DO 9008 I = 1,N
       if( z(i).ge.0.15 ) then
       E1(I) = PHIM0( I1(I) )
       E2(I) = PHIM0( I2(I) )
       PHIM(I)  = TEMP1(I) * ( E2(I)-E1(I) )
       PHIM(I)  = PHIM(I) +   E1(I)
       endif
 9008  CONTINUE

  100 CONTINUE
!
      IF( IFLAG .EQ. 2 ) GO TO 200
       DO 9010 I = 1,N
       if( z(i).ge.0.15 ) then
       E1(I) = PHIH0( I1(I) )
       E2(I) = PHIH0( I2(I) )
       PHIH(I)  = TEMP1(I) * ( E2(I)-E1(I) )
       PHIH(I)  = PHIH(I) +   E1(I)
       endif
 9010  CONTINUE

  200 CONTINUE
!
       DO 9012 I = 1,N
       ZSTAR(I) = -Z(I)
 9012  CONTINUE
!
      IF( IFLAG .GT. 2 ) GO TO 300
       DO 9014 I = 1,N
       IF( Z(I) .LT. 0.15 ) PHIM(I) = 1. + ZSTAR(I) &
           *(0.25+ZSTAR(I)*(0.09375+ZSTAR(I)* &
           (0.03125+0.00732422 * ZSTAR(I))))
 9014  CONTINUE
!
  300 CONTINUE
      IF( IFLAG .EQ. 2 ) GO TO 500
       DO 9016 I = 1,N
       IF( Z(I) .LT. 0.15 ) THEN
       PHIH(I) =1.+ Z(I) * (0.5+ZSTAR(I)*(0.375+ZSTAR(I)* &
           (0.5+ZSTAR(I)*(0.8203125+ZSTAR(I)* &
           (1.5+2.93262*ZSTAR(I))))))
       PHIH(I) = 1. / PHIH(I)
      ENDIF
 9016  CONTINUE
!
  500 CONTINUE

end subroutine phi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: psi
! !INTERFACE:
      SUBROUTINE PSI(VZZ,VZH,VPSIM,VPSIH,IRUN,VX,VXS,VY,VYS,IFLAG)
!**********************************************************************
!
!  SUBROUTINE PSI - DETERMINES DIMENSIONLESS WIND AND
!                    SCALAR PROFILES IN SURFACE LAYER
!                 - CALLED FROM helfsurface
!
!  DESCRIPTION OF PARAMETERS
!     ZZ   -  INPUTED VALUE OF MONIN- OBUKHOV STABILITY PARAMETER ZETA
!     ZH   -  INPUTED VALUE OF Z0 DIVIDED BY SFC LAYER HEIGHT
!     PSIM -  OUTPUTED VALUE OF DIMENSIONLESS WIND
!     PSIH -  OUTPUTED VALUE OF DIMENSIONLESS SCALAR
!     X    -  OUTPUTED VALUE OF PHIM(ZETA)
!     XS   -  OUTPUTED VALUE OF PHIM(ZETA0)
!     Y    -  OUTPUTED VALUE OF PHIH(ZETA)
!     YS   -  OUTPUTED VALUE OF PHIH(ZETA0)
!     IFLAG-  FLAG TO DETERMINE IF CU IS NEEDED (IFLAG=2),
!                  IF CT IS NEEDED (IFLAG=3), OR BOTH (IFLAG=1)
!  SUBPROGRAMS NEEDED
!     PHI  -  COMPUTES SIMILARITY FUNCTION FOR MOMENTUM AND SCALARS
!
!**********************************************************************
      implicit none

! Argument List Declarations
      integer irun,iflag
      real VZZ(:),VZH(:),VPSIM(:),VPSIH(:), &
           VX(:),VXS(:),VY(:),VYS(:)
 
! Local Variables
      real ZWM,RZWM,Z0M,ZCM,RZCM,CM1,CM2,CM6,CM7,CM8ARG,YCM
      PARAMETER ( ZWM     =    1.    )
      PARAMETER ( RZWM    =  1./ZWM  )
      PARAMETER ( Z0M     =    0.2    )
      PARAMETER ( ZCM     =    42.    )
      PARAMETER ( RZCM    =  1./ZCM  )
      PARAMETER ( CM1     =  1./126. )
      PARAMETER ( CM2     =  1./(6.*CM1)  )
      PARAMETER ( CM6     =  6. / ( 1. + 6.*CM1 )  )
      PARAMETER ( CM7     =  CM2 + ZWM  )
      PARAMETER ( CM8ARG  =  CM7*ZCM*RZWM / (CM2+ZCM)  )
      PARAMETER ( YCM     =  6. / ( 1. + 6.*CM1*ZCM )  )

      integer INTSTB(irun),INTZ0(irun)
      real ZZ0(irun),Z(irun),Z2(irun),Z1(irun),Z0(irun)
      real X0(irun),X1(irun),Y0(irun),Y1(irun)
      real PSI2(irun),TEMP(irun)
      real HZ(irun),ARG0(irun),ARG1(irun),DX(irun)
      real X0NUM(irun),X1NUM(irun),X0DEN(irun)
      real X1DEN(irun),Y1DEN(irun),Z2ZWM(irun)
      real cm3,cm4,cm5,cm8
      integer ibit,indx
      integer i
!
      CM3 =   sqrt( 0.2/CM1-0.01 )
      CM4 =   1./CM3
      CM5 =  (10.-CM1) / (10.*CM1*CM3)
      CM8 =   6. * LOG(CM8ARG)
!
      DO 9000 I = 1,IRUN
       VPSIM(I) = 0.
       VPSIH(I) = 0.
       VX(I) = 0.
       VXS(I) = 0.
       VY(I) = 0.
       VYS(I) = 0.
       ZZ0(I) = VZH(I)*VZZ(I)
 9000 CONTINUE
      IBIT = 0
      DO 9122 I = 1,IRUN
       IF(VZZ(I).LE.-1.e-7)IBIT = IBIT + 1
 9122 CONTINUE
      DO 9022 I = 1,IRUN
       IF(VZZ(I).LE.-1.e-7)THEN
        INTSTB(I) = 1
       ELSE
        INTSTB(I) = 0
       ENDIF
 9022 CONTINUE
!
! ****************************************
! *****    UNSTABLE SURFACE LAYER    *****
! ****************************************
!
      IF(IBIT.LE.0)  GO TO 100
!
      indx = 0
      DO 9002 I = 1,IRUN
       IF (INTSTB(I).EQ.1)THEN
        indx = indx + 1
        Z(indx) = VZZ(I)
        Z0(indx) = ZZ0(I)
       ENDIF
 9002 CONTINUE
!
      DO 9004 I = 1,IBIT
       Z(I) = -18. * Z(I)
       Z0(I) = -18. * Z0(I)
 9004 CONTINUE
 
      CALL PHI( Z,X1,Y1,IFLAG,IBIT )
      CALL PHI( Z0,X0,Y0,IFLAG,IBIT )
 
! ****************************
! *****    COMPUTE PSIM  *****
! ****************************
!
      IF(IFLAG.GE.3) GO TO 75
!
      DO 9006 I = 1,IBIT
       ARG1(I) = 1. - X1(I)
       IF ( Z(I) .LT. 0.013 ) ARG1(I) = Z(I) * ( 0.25 -  0.09375 * Z(I) )
!
       ARG0(I)  = 1. - X0(I)
       IF ( Z0(I) .LT. 0.013 ) ARG0(I) = Z0(I) * ( 0.25 -  0.09375 * Z0(I) )
!
       ARG1(I) = ARG1(I) * ( 1.+X0(I) )
       ARG0(I) = ARG0(I) * ( 1.+X1(I) )
       DX(I) = X1(I) - X0(I)
       ARG1(I) = ARG1(I) / ARG0(I)
       ARG0(I) = -DX(I) / ( 1. + X1(I)*X0(I) )
       ARG0(I) = ATAN( ARG0(I) )
       ARG1(I) = LOG( ARG1(I) )
       PSI2(I) = 2. * ARG0(I) + ARG1(I)
       PSI2(I) = PSI2(I) + DX(I)
 9006 CONTINUE
!
      indx = 0
      DO 9008 I = 1,IRUN
       IF( INTSTB(I).EQ.1 ) THEN
        indx = indx + 1
        VPSIM(I) = PSI2(indx)
        VX(I) = X1(indx)
        VXS(I) = X0(indx)
       ENDIF
 9008 CONTINUE
!
! ****************************
! *****    COMPUTE PSIH  *****
! ****************************
!
      IF(IFLAG.EQ.2) GO TO 100
!
  75  CONTINUE
      DO 9010 I = 1,IBIT
       ARG1(I) = 1. - Y1(I)
       IF( Z(I) .LT. 0.0065 ) ARG1(I) = Z(I) * ( 0.5 -  0.625 * Z(I) )
!
       ARG0(I)  = 1. - Y0(I)
       IF( Z0(I) .LT. 0.0065 ) ARG0(I) = Z0(I) * ( 0.5 -  0.625 * Z0(I) )
!
       ARG1(I) = ARG1(I) * ( 1. + Y0(I) )
       ARG0(I) = ARG0(I) * ( 1. + Y1(I) )
       ARG1(I) = ARG1(I) / ARG0(I)
       PSI2(I) = LOG( ARG1(I) )
       PSI2(I) = PSI2(I) - Y1(I) + Y0(I)
 9010 CONTINUE
!
      indx = 0
      DO 9012 I = 1,IRUN
       IF( INTSTB(I).EQ.1 ) THEN
       indx = indx + 1
       VPSIH(I) = PSI2(indx)
       VY(I) = Y1(indx)
       VYS(I) = Y0(indx)
       ENDIF
 9012 CONTINUE
!
! **************************************
! *****    STABLE SURFACE LAYER    *****
! **************************************
!
  100 CONTINUE
      IBIT = 0
      DO 9114 I = 1,IRUN
       IF(VZZ(I).GT.-1.e-7)THEN
        IBIT = IBIT + 1
       ENDIF
 9114 CONTINUE
      DO 9014 I = 1,IRUN
       IF(VZZ(I).GT.-1.e-7)THEN
        INTSTB(I) = 1
       ELSE
        INTSTB(I) = 0
       ENDIF
 9014 CONTINUE
      IF(IBIT.LE.0)  GO TO 300
      indx = 0
      DO 9016 I = 1,IRUN
       IF (INTSTB(I).EQ.1)THEN
        indx = indx + 1
        Z(indx) = VZZ(I)
        Z0(indx) = ZZ0(I)
        ARG1(indx) = VZH(I)
       ENDIF
 9016 CONTINUE

      DO 9018 I = 1,IBIT
       HZ(I) = 1. / ARG1(I)
       Z1(I) = Z(I)
       Z2(I) = ZWM
!
       IF ( Z(I) .GT. ZWM ) THEN
        Z1(I) = ZWM
        Z2(I) = Z(I)
       ENDIF
!
       IF ( Z0(I) .GT. Z0M ) THEN
        Z0(I) = Z0M
        INTZ0(I) = 1
       ELSE
        INTZ0(I) = 0
       ENDIF
!
       X1NUM(I) = 1. + 5. * Z1(I)
       X0NUM(I) = 1. + 5. * Z0(I)
       X1DEN(I) = 1. / (1. + CM1 * (X1NUM(I) * Z1(I)) )
       X0DEN(I) = 1. + CM1 * (X0NUM(I) * Z0(I))
!
       IF ( (INTZ0(I).EQ.1) .OR. (Z(I).GT.ZWM) ) &
            HZ(I) = Z1(I) / Z0(I)
       ARG1(I) = HZ(I)*HZ(I)*X0DEN(I)*X1DEN(I)
       ARG1(I) = LOG( ARG1(I) )
       ARG1(I) = 0.5 * ARG1(I)
       ARG0(I) = (Z1(I) + 0.1) * (Z0(I) + 0.1)
       ARG0(I) = CM3 + ARG0(I) * CM4
       ARG0(I) = ( Z1(I) - Z0(I) ) / ARG0(I)
       ARG0(I) = ATAN( ARG0(I) )
       TEMP(I) = ARG1(I) + CM5 * ARG0(I)
!
       X0(I) = X0NUM(I) / X0DEN(I)
       IF ( INTZ0(I).EQ.1 ) X0(I) = 0.
       Z2ZWM(I) = Z2(I) * RZWM
 9018 CONTINUE
!
! ****************************
! *****    COMPUTE PSIM  *****
! ****************************
!
      IF( IFLAG.GE.3 ) GO TO 225
!
      DO 9020 I = 1,IBIT
       X1(I) = X1NUM(I) * X1DEN(I)
       ARG1(I) = LOG( Z2ZWM(I) )
       PSI2(I) = TEMP(I) + CM6 * ARG1(I)
 9020 CONTINUE
!
      indx = 0
      DO 9030 I = 1,IRUN
       IF( INTSTB(I).EQ.1 ) THEN
       indx = indx + 1
       VPSIM(I) = PSI2(indx)
       VX(I) = X1(indx)
       VXS(I) = X0(indx)
       ENDIF
 9030 CONTINUE
!
! ****************************
! *****    COMPUTE PSIH  *****
! ****************************
!
       IF(IFLAG.EQ.2)GO TO 300
!
  225 CONTINUE
      DO 9024 I = 1,IBIT
       Y1DEN(I) = 1. + CM1 * ( X1NUM(I) * Z(I) )
       Y1(I) = X1NUM(I) / Y1DEN(I)
       ARG1(I) = CM7 * Z2ZWM(I) / ( CM2 + Z2(I) )
       ARG0(I) = 6.
       IF ( Z2(I) .GT. ZCM ) THEN
        Y1(I) = YCM
        ARG1(I) = Z2(I) * RZCM
        ARG0(I) = YCM
        TEMP(I) = TEMP(I) + CM8
       ENDIF
       ARG1(I) = LOG( ARG1(I) )
       PSI2(I) = TEMP(I) + ARG0(I) * ARG1(I)
 9024 CONTINUE
!
      indx = 0
      DO 9026 I = 1,IRUN
       IF( INTSTB(I).EQ.1 ) THEN
       indx = indx + 1
       VPSIH(I) = PSI2(indx)
       VY(I) = Y1(indx)
       VYS(I) = X0(indx)
       ENDIF
 9026 CONTINUE
!
  300 CONTINUE
!
end subroutine psi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: linadj
! !INTERFACE:
      SUBROUTINE LINADJ ( VRIB1,VRIB2,VWS1,VWS2,VZ1,VUSTAR,IWATER, &
       VAPSIM, VAPSIHG,VPSIH,VPSIG,VX,VX0,VY,VY0,ITYPE,LWATER,IRUN, &
       VDZETA,VDZ0,VDPSIM,VDPSIH,INTRIB, &
       VX0PSIM,VG,VG0,VR1MG0,VZ2,VDZSEA,VAZ0,VXNUM1,VPSIGB2,VDX, &
       VDXPSIM,VDY,VXNUM2,VDEN,VAWS1,VXNUM3,VXNUM,VDZETA1,VDZETA2, &
       VZCOEF2,VZCOEF1,VTEMPLIN,VDPSIMC,VDPSIHC,vk,bmdl,CHOOSEZ0,VCHARNOCK)
!
!**********************************************************************
!
!  ARGUMENTS ::
!
!     INPUT:
!     ------
!    RIB1          -         BULK RICHARDSON NUMBER OF INPUT STATE
!    RIB2          -         DESIRED BULK RICH NUMBER OF OUTPUT STATE
!    WS1           -         SURFACE WIND SPEED OF INPUT STATE
!    WS2           -         DESIRED SURFACE WIND SPEED OF OUTPUT STATE
!    Z1            -         INPUT VALUE OF ROUGHNESS HEIGHT
!    USTAR         -         INPUT VALUE OF CU * WS
!    WATER         -         BIT ARRAY - '1' WHERE OCEAN
!    APSIM         -         (1/PSIM)
!    APSIHG        -         ( 1 / (PSIH+PSIG) )
!    PSIH          -         NON-DIM TEMP GRADIENT
!    PSIG          -         PSIH FOR THE MOLECULAR LAYER
!    X             -         PHIM(ZETA) - DERIVATIVE OF PSIM
!    X0            -         PHIM(ZETA0)
!    Y             -         PHIH(ZETA) - DERIVATIVE OF PSIH
!    Y0            -         PHIH(ZETA0)
!    ITYPE         -         INTEGER FLAG :
!                               1    = NEUTRAL ADJUSTMENT
!                               3, 5 = ADJUSTMENT INSIDE LOOP
!                               5    = ADJUST CU AND CT
!    LWATER        -         LOGICAL - .TRUE. IF THERE ARE WATER POINTS
!    CHOOSEZ0      -         INTEGER FLAG: 0 - L&P Z0, no high wind limit
!                                          1 - Edson Z0 for mom. and heat, high wind limit
!                                          2 - L&P Z0, high wind limit
!                                          3 - Edson Z0 for mom. only, high wind limit
!                                          4 - wave model Charnock coefficient
!
!     OUTPUT:
!     -------
!    DZETA         -         D LOG ZETA
!    DZ0           -         D Z0 (ITYPE 1) OR D LOG Z0 (ITYPE 2-5)
!    DPSIM         -         D PSIM
!    DPSIH         -         D PSIH
!    BITRIB        -         BIT ARRAY - '1' WHERE RIB1 = 0
!
!**********************************************************************
      implicit none

! Argument List Declarations
      integer irun,itype,CHOOSEZ0
      real VRIB1(:),VRIB2(:)
      real VWS1(:),VWS2(:),VZ1(:),VUSTAR(:)
      integer IWATER(:)
      real VAPSIM(:),VAPSIHG(:)
      real VPSIH(:),VPSIG(:),VX(:)
      real VX0(:),VY(:),VY0(:)
      LOGICAL LWATER
      real VDZETA(:),VDZ0(:),VDPSIM(:)
      real VDPSIH(:)
      integer INTRIB(:)
      real VX0PSIM(:),VG(:),VG0(:),VR1MG0(:)
      real VZ2(:),VDZSEA(:),VAZ0(:),VXNUM1(:)
      real VPSIGB2(:),VDX(:),VDXPSIM(:),VDY(:)
      real VXNUM2(:),VDEN(:),VAWS1(:),VXNUM3(:)
      real VXNUM(:),VDZETA1(:),VDZETA2(:)
      real VZCOEF2(:),VZCOEF1(:),VTEMPLIN(:)
      real VDPSIMC(:),VDPSIHC(:),bmdl(:)
      real VCHARNOCK(:)

! Local Variables
      real xx0max,prfac,xpfac,difsqt,ustz0s,h0byz0,usth0s
      PARAMETER ( XX0MAX  =   1.49821 )
      PARAMETER ( PRFAC  = 0.595864   )
      PARAMETER ( XPFAC  = .55        )  
      PARAMETER ( DIFSQT  = 3.872983E-3)
      PARAMETER ( USTZ0S =   0.2030325E-5)
      PARAMETER ( H0BYZ0 =    30.0    )
      PARAMETER ( USTH0S =  H0BYZ0*USTZ0S )

      integer VINT1(irun),VINT2(irun)
      real vk,b2uhs(irun)
      integer i
!
      UNUSED_DUMMY(VWS2)
      UNUSED_DUMMY(vk)

      do i = 1,irun
      B2UHS(i)   = BMDL(i) * BMDL(i) * USTH0S
      enddo

!   COMPUTE X0/PSIM, 1/Z0, G, G0, 1/(1-G0),
!     DEL LOG Z0, D LOG ZO / D USTAR
!
      IF ( (ITYPE.EQ.1) .AND. LWATER ) THEN
       DO 9000 I = 1,IRUN
        IF (IWATER(I).EQ.1) VX0PSIM(I) = VAPSIM(I)
 9000  CONTINUE
      ENDIF
      IF ( ITYPE .GE. 3 ) THEN
       DO 9002 I = 1,IRUN
        VX0PSIM(I) = VX0(I) * VAPSIM(I)
 9002  CONTINUE
      ENDIF
!
       DO 9004 I = 1,IRUN
        VDZ0(I) = 0.
        VG(I) = 0.
        VG0(I) = 0.
        VR1MG0(I) = 1.
 9004  CONTINUE
!
       IF ( LWATER ) THEN
        CALL ZCSUB ( VUSTAR,VCHARNOCK,VDZSEA,IWATER,.TRUE.,IRUN,VZ2,CHOOSEZ0)

        VDZSEA = min( VDZSEA, 0.2*VZ1/VAPSIM ) ! To prevent Divide by Zero as VG0 => 1.0
!
        DO 9006 I = 1,IRUN
         IF ( IWATER(I).EQ.1) THEN
          VAZ0(I) = 1. / VZ1(I)
          VG(I) = VDZSEA(I) * VAZ0(I)
          VG0(I) = VX0PSIM(I) * VG(I)
          VR1MG0(I) = 1. / ( 1. - VG0(I) )
          VDZ0(I) = ( VZ2(I) - VZ1(I) ) * VR1MG0(I)
         ENDIF
 9006   CONTINUE
       ENDIF
!
      IF ( LWATER .AND. (ITYPE.GE.3) ) THEN
       DO 9008 I = 1,IRUN
        IF (IWATER(I).EQ.1) VDZ0(I) = VDZ0(I) * VAZ0(I)
 9008  CONTINUE
      ENDIF
!
!   COMPUTE NUM1,NUM2,NUM3, DEN
!
      IF (ITYPE.GE.3) THEN
       DO 9010 I = 1,IRUN
        VXNUM1(I) = 0.
        IF (VRIB1(I).EQ.0.) THEN
         INTRIB(I) = 1
        ELSE
         INTRIB(I) = 0
        ENDIF
        IF ( INTRIB(I).EQ.0 ) VXNUM1(I) = 1. / VRIB1(I)
        VPSIGB2(I) = 0.
        if(vpsig(i).gt.0.)VPSIGB2(I) = &
              0.5 * ( vpsig(i)*vpsig(i) + b2uhs(i) ) / vpsig(i)
        VDX(I) = VX(I) - VX0(I)
        VDXPSIM(I) = VDX(I) * VAPSIM(I)
        VDY(I) = VY(I) - VY0(I)
        VXNUM3(I) = - VPSIGB2(I)
!
        IF ( LWATER ) THEN
         IF (IWATER(I).EQ.1) THEN
          VDXPSIM(I) = VDXPSIM(I) * VR1MG0(I)
          VXNUM3(I) = VXNUM3(I) + VG(I) * ( VY0(I) - VPSIGB2(I) )
          VXNUM2(I) = VY0(I) - VPSIGB2(I) - VX0PSIM(I) * VPSIGB2(I)
          VXNUM2(I) = (VXNUM2(I) * VAPSIHG(I)) - 2. * VX0PSIM(I)
          VXNUM2(I) = VXNUM2(I) * VDZ0(I)
         ENDIF
        ENDIF
!
        VDEN(I) = VDY(I) + VDXPSIM(I) * VXNUM3(I)
        VDEN(I) = ( 1. + VDEN(I) * VAPSIHG(I) ) - 2. * VDXPSIM(I)
 9010  CONTINUE
      ENDIF
!
      IF (ITYPE.EQ.5) THEN
       DO 9012 I = 1,IRUN
        VAWS1(I) = VR1MG0(I) / VWS1(I)
        VXNUM3(I) = VXNUM3(I) * VAPSIHG(I)
!
        IF ( LWATER ) THEN
         IF(IWATER(I).EQ.1) THEN
          VXNUM3(I) = VXNUM3(I) - 2. * VG0(I)
          VXNUM3(I) = VAWS1(I) * VXNUM3(I)
         ENDIF
        ENDIF
 9012  CONTINUE
      ENDIF
!
!   COMPUTE D LOG ZETA
!
      IF (ITYPE.GE.3) THEN
       DO 9014 I = 1,IRUN
        VXNUM(I) = VRIB2(I) - VRIB1(I)
        IF( (VX0(I).GT.XX0MAX).AND.(VXNUM(I).GE.0.) )VXNUM(I) = 0.
        VXNUM(I) = VXNUM1(I) * VXNUM(I)
 9014  CONTINUE
!
       DO 9018 I = 1,IRUN
        VDZETA1(I) = VXNUM(I)
        IF(LWATER.AND.(IWATER(I).EQ.1)) VXNUM(I) = VXNUM(I) + VXNUM2(I)
        IF ( VDEN(I) .LT.0.1 ) VDEN(I) = 0.1
 9018  CONTINUE
!
       DO 9020 I = 1,IRUN
        VDZETA(I) = VXNUM(I) / VDEN(I)
 9020  CONTINUE
       DO 9022 I = 1,IRUN
        IF((VRIB2(I).EQ.0.).OR.(VDZETA(I).LE.-1.))VDZETA(I) = VDZETA1(I)
 9022  CONTINUE
      ENDIF
!
!   COMPUTE D LOG Z0
!
      IF ( LWATER .AND. (ITYPE.GE.3) )THEN
       DO 9026 I = 1,IRUN
        IF( IWATER(I).EQ.1 ) THEN
         VZCOEF2(I) = VG(I) * VDXPSIM(I)
         VDZ0(I) = VDZ0(I) - VZCOEF2(I) * VDZETA(I)
        ENDIF
 9026  CONTINUE
      ENDIF
!
      IF ( LWATER .AND. (ITYPE.EQ.5) ) THEN
       DO 9028 I = 1,IRUN
        IF(IWATER(I).EQ.1) VZCOEF1(I) = VG(I) * VAWS1(I)
 9028  CONTINUE
      ENDIF
!
!   CALCULATE D PSIM AND D PSIH
!
      IF ( (ITYPE.EQ.1) .AND. LWATER ) THEN
       DO 9032 I = 1,IRUN
        IF (IWATER(I).EQ.1) THEN
         VDPSIM(I) = - VDZ0(I) * VAZ0(I)
         VDPSIH(I) = VDPSIM(I)
        ENDIF
 9032  CONTINUE
      ENDIF
!
      IF (ITYPE.GE.3) THEN
       DO 9034 I = 1,IRUN
        VDPSIM(I) = VDX(I) * VDZETA(I)
        VDPSIH(I) = VDY(I) * VDZETA(I)
        IF ( LWATER ) THEN
         IF (IWATER(I).EQ.1 ) THEN
          VDPSIM(I) = VDPSIM(I) - VX0(I) * VDZ0(I)
          VDPSIH(I) = VDPSIH(I) - VY0(I) * VDZ0(I)
         ENDIF
        ENDIF
 9034  CONTINUE
      ENDIF
!
!   PREVENT OVERCORRECTION OF PSIM OR PSIH FOR UNSTABLE CASE
!
      IF (ITYPE.GE.4) THEN
       DO 9036 I = 1,IRUN
        VDPSIMC(I) = -0.9 - VDPSIM(I) * VAPSIM(I)
        VDPSIHC(I) = -0.9 *  VPSIH(I) - VDPSIH(I)
        IF ( VDPSIMC(I).GT.0.  ) THEN
         VINT1(I) = 1
        ELSE
         VINT1(I) = 0
        ENDIF
        IF ( VDPSIHC(I).GT.0.  ) THEN
         VINT2(I) = 1
        ELSE
         VINT2(I) = 0
        ENDIF
        VDZETA1(I) = 0.
        IF(VINT1(I).EQ.1) VDZETA1(I) = VDPSIMC(I) / VDXPSIM(I)
        IF((VINT1(I).EQ.1).OR.(VINT2(I).EQ.1)) VTEMPLIN(I) = &
              VDY(I) + VY0(I) * VG(I) * VDXPSIM(I)
!AMM    IF (VINT2(I).EQ.1 .and. VTEMPLIN(I).GT.tiny(1.0)) then
        IF (VINT2(I).EQ.1) then
             VDZETA2(I) =  VDPSIHC(I) / VTEMPLIN(I)
        IF ( VDZETA2(I).LT.VDZETA1(I) ) VDZETA1(I) = VDZETA2(I)
        endif
        IF((VINT1(I).EQ.1).OR.(VINT2(I).EQ.1)) THEN
         VDZETA(I) = VDZETA1(I) + VDZETA(I)
         VDPSIM(I) = VDPSIM(I) + VDX(I) * VR1MG0(I) * VDZETA1(I)
         VDPSIH(I) = VDPSIH(I) + VTEMPLIN(I) * VDZETA1(I)
         IF ( IWATER(I).EQ.1 ) &
           VDZ0(I) = VDZ0(I) - VG(I) * VDXPSIM(I) * VDZETA1(I)
        ENDIF
 9036  CONTINUE
      ENDIF
!
end subroutine linadj
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: zcsub
! !INTERFACE:
      SUBROUTINE ZCSUB (VUSTAR,VCHARNOCK,VDZSEA,IWATER,LDZSEA,IRUN,VZSEA,CHOOSEZ0)
!**********************************************************************
!  FUNCTION ZSEA
!  PURPOSE
!     COMPUTES Z0 AS A FUNCTION OF USTAR OVER WATER SURFACES
!  USAGE
!     CALLED BY helfsurface
!  DESCRIPTION OF PARAMETERS
!     USTAR    -  INPUTED VALUE OF SURFACE-STRESS VELOCITY
!     DZSEA    -  OUTPUTED VALUE OF DERIVATIVE  D(ZSEA)/D(USTAR)
!     WATER    -  INPUTED BIT VECTOR TO DETERMINE WATER POINTS
!     LDZSEA   -  LOGICAL FLAG TO DETERMINE IF DZSEA SHOULD BE COMPUTED
!     ZSEA     -  OUTPUTED VALUE OF ROUGHNESS LENGTH
!     CHOOSEZ0 -  INTEGER FLAG: 0 - L&P Z0, no high wind limit
!                               1 - Edson Z0 for mom. and heat, high wind limit
!                               2 - L&P Z0, high wind limit
!                               3 - Edson Z0 for mom. only, high wind limit
!                               4 - wave-model Charnock coefficient
!  SUBPROGRAMS NEEDED
!     NONE
!  RECORD OF MODIFICATIONS
!   Molod 6/8/2011 - Implement new choozez0 options (expand from 0,1 choice)
!  REMARKS:
!        COMPUTE ROUGHNESS LENGTH FOR OCEAN POINTS
!          BASED ON FUNCTIONS OF LARGE AND POND
!          AND OF KONDO --- DESIGNED FOR K = .4
! *********************************************************************
      implicit none 

! Argument List Delcarations
      integer irun, CHOOSEZ0
      real VZSEA(:),VUSTAR(:),VDZSEA(:),VCHARNOCK(:)
      integer IWATER(:)
      LOGICAL LDZSEA

! Local Variables
      real USTMX1_OLD,USTMX2_OLD
      real USTMX1_NEW,USTMX2_NEW
      real USTMX1,USTMX2,USTMX3

      PARAMETER ( USTMX1_NEW =   0.80 )
      PARAMETER ( USTMX2_NEW =   0.80 )
      PARAMETER ( USTMX1_OLD =   1.1  )
      PARAMETER ( USTMX2_OLD =   0.381844 )
      PARAMETER ( USTMX3     =   0.0632456)

      real AA(IRUN,5),TEMP(IRUN)
      integer INT2(IRUN),INT3(IRUN),INT4(IRUN)
      integer i,k
      real ustloc(irun)

      real AA1(5),AA2(5),AA3(5),AA4(5)
      real AA2_NEW(5),AA3_NEW(5),AA4_NEW(5)
      real AA2_OLD(5),AA3_OLD(5),AA4_OLD(5)

      DATA AA1/.2030325E-5,0.0,0.0,0.0,0.0/

      DATA AA2_NEW/-1.102451E-08,0.1593E-04,0.1E-03,2.918E-03, &
               0.695649E-04/
      DATA AA3_NEW/-1.102451E-08,0.12E-04,0.1E-03,2.918E-03, &
               1.5649E-04/
      DATA AA4_NEW/0.085E-03,1.5E-03,-0.210E-03,0.215E-02, &
               -0.0/

      DATA AA2_OLD/-0.402451E-08,0.239597E-04,0.117484E-03,0.191918E-03, &
               0.395649E-04/
      DATA AA3_OLD/-0.237910E-04,0.228221E-03,-0.860810E-03,0.176543E-02, &
               0.784260E-04/
      DATA AA4_OLD/-0.343228E-04,0.552305E-03,-0.167541E-02,0.250208E-02, &
               -0.153259E-03/

      CHARNOCK: if ( CHOOSEZ0 == 4 ) then
          ustloc = max(1e-6, vustar)
          VZSEA = (0.11*MAPL_NUAIR)/ustloc + (VCHARNOCK/MAPL_GRAV)*ustloc**2
          
          DERIVATIVE: if ( LDZSEA ) then
              VDZSEA = -(0.11*MAPL_NUAIR)/ustloc**2 + (VCHARNOCK/MAPL_GRAV)*2*ustloc
          end if DERIVATIVE

          return
      end if CHARNOCK

      if( CHOOSEZ0.eq.0 .OR. CHOOSEZ0.eq.2) then
          USTMX1 = USTMX1_OLD
          USTMX2 = USTMX2_OLD
             AA2 =    AA2_OLD
             AA3 =    AA3_OLD
             AA4 =    AA4_OLD
      else
          USTMX1 = USTMX1_NEW
          USTMX2 = USTMX2_NEW
             AA2 =    AA2_NEW
             AA3 =    AA3_NEW
             AA4 =    AA4_NEW
      endif
!
!**********************************************************************
!*****              LOWER CUTOFF CONDITION FOR USTAR                ***
!**********************************************************************
!
      DO 9000 I = 1,IRUN
       IF(VUSTAR(I) .LT. 1.e-6)THEN
        INT3(I) = 1
       ELSE
        INT3(I) = 0
       ENDIF
 9000 CONTINUE
      DO 9002 I = 1,IRUN
       IF(INT3(I).EQ.1) VUSTAR(I) = 1.e-6
 9002 CONTINUE

      ustloc = vustar
!
!***********************************
!*****  LOAD THE ARRAY A(I,K)  *****
!***********************************
!
      DO 9004 I = 1,IRUN
       IF( (ustloc(I) .GT. USTMX1) .AND. (IWATER(I).EQ.1) ) THEN
        if( CHOOSEZ0.gt.0 ) ustloc(i) = ustmx1
        INT4(I) = 1
       ELSE
        INT4(I) = 0
       ENDIF
 9004 CONTINUE
      DO 9006 I = 1,IRUN
       IF(ustloc(I) .GT. USTMX2) THEN
        INT3(I) = 1
       ELSE
        INT3(I) = 0
       ENDIF
 9006 CONTINUE
      DO 9008 I = 1,IRUN
       IF(ustloc(I) .GE. USTMX3) THEN
        INT2(I) = 1
       ELSE
        INT2(I) = 0
       ENDIF
 9008 CONTINUE
!
      DO 100 K=1,5
       DO 9010 I = 1,IRUN
        AA(I,K) = AA1(K)
        IF( INT2(I).EQ.1 )  AA(I,K) = AA2(K)
        IF( INT3(I).EQ.1 )  AA(I,K) = AA3(K)
        IF( INT4(I).EQ.1 )  AA(I,K) = AA4(K)
 9010  CONTINUE
  100 CONTINUE
!
!********************************************************
!*****  EVALUATE THE ENHANCED POLYNOMIAL FOR ZSEA  *****
!********************************************************
!
      DO 9012 I = 1,IRUN
       VDZSEA(I)  =  ( AA(I,4) + AA(I,5) * ustloc(I) ) * ustloc(I)
       VZSEA(I)  =  AA(I,2) + ( AA(I,3) + VDZSEA(I) ) * ustloc(I)
       TEMP(I) = AA(I,1) / ustloc(I)
       VZSEA(I)  =  VZSEA(I) + TEMP(I)
 9012 CONTINUE
!
!**********************************************************************
!*****        EVALUATE THE DERIVATIVE DZSEA IF LDZSEA IS TRUE       ***
!**********************************************************************
!
      IF( LDZSEA ) THEN
       DO 9014 I = 1,IRUN
        VDZSEA(I)  =  3. * VDZSEA(I) -(AA(I,4)*ustloc(I) - AA(I,3))
        VDZSEA(I)  =  VDZSEA(I) * ustloc(I) - TEMP(I)
 9014  CONTINUE
      ENDIF
!
end subroutine zcsub


!=======================================================================


!=======================================================================
      end Program surf_layer 


!=======================================================================
