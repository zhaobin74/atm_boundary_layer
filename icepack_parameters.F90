!=========================================================================
!
! flags for the column package
!
! authors: Elizabeth C. Hunke, LANL

      module icepack_parameters

      use icepack_kinds

      implicit none
      private


      !-----------------------------------------------------------------
      ! control options
      !-----------------------------------------------------------------

      character (char_len), public :: &
         argcheck = 'first'      !  optional argument checks, 'never','first','always'

      !-----------------------------------------------------------------
      ! parameter constants
      !-----------------------------------------------------------------

      real (kind=dbl_kind), parameter, public :: &
         c0   = 0.0_dbl_kind, &
         c1   = 1.0_dbl_kind, &
         c1p5 = 1.5_dbl_kind, &
         c2   = 2.0_dbl_kind, &
         c3   = 3.0_dbl_kind, &
         c4   = 4.0_dbl_kind, &
         c5   = 5.0_dbl_kind, &
         c6   = 6.0_dbl_kind, &
         c8   = 8.0_dbl_kind, &
         c10  = 10.0_dbl_kind, &
         c15  = 15.0_dbl_kind, &
         c16  = 16.0_dbl_kind, &
         c20  = 20.0_dbl_kind, &
         c25  = 25.0_dbl_kind, &
         c100 = 100.0_dbl_kind, &
         c180 = 180.0_dbl_kind, &
         c1000= 1000.0_dbl_kind, &
         p001 = 0.001_dbl_kind, &
         p01  = 0.01_dbl_kind, &
         p1   = 0.1_dbl_kind, &
         p2   = 0.2_dbl_kind, &
         p4   = 0.4_dbl_kind, &
         p5   = 0.5_dbl_kind, &
         p6   = 0.6_dbl_kind, &
         p05  = 0.05_dbl_kind, &
         p15  = 0.15_dbl_kind, &
         p25  = 0.25_dbl_kind, &
         p75  = 0.75_dbl_kind, &
         p333 = c1/c3, &
         p666 = c2/c3, &
         spval_const= -1.0e36_dbl_kind

      real (kind=dbl_kind), public :: &
         secday = 86400.0_dbl_kind ,&! seconds in calendar day
         puny   = 1.0e-11_dbl_kind, &
         bignum = 1.0e+30_dbl_kind, &
         pi     = 3.14159265358979323846_dbl_kind

      !-----------------------------------------------------------------
      ! derived physical constants
      !    Lfresh = Lsub-Lvap     ,&! latent heat of melting of fresh ice (J/kg)
      !    cprho  = cp_ocn*rhow   ,&! for ocean mixed layer (J kg / K m^3)
      !    Cp     = 0.5_dbl_kind*gravit*(rhow-rhoi)*rhoi/rhow ,&! proport const for PE
      !-----------------------------------------------------------------

      real (kind=dbl_kind), public :: &
         pih        = spval_const     ,&! 0.5 * pi
         piq        = spval_const     ,&! 0.25 * pi
         pi2        = spval_const     ,&! 2 * pi
         rad_to_deg = spval_const     ,&! conversion factor, radians to degrees
         Lfresh     = spval_const     ,&! latent heat of melting of fresh ice (J/kg)
         cprho      = spval_const     ,&! for ocean mixed layer (J kg / K m^3)
         Cp         = spval_const       ! proport const for PE

      !-----------------------------------------------------------------
      ! Densities
      !-----------------------------------------------------------------

      real (kind=dbl_kind), public :: &
         rhos      = 330.0_dbl_kind   ,&! density of snow (kg/m^3)
         rhoi      = 917.0_dbl_kind   ,&! density of ice (kg/m^3)
         rhosi     = 940.0_dbl_kind   ,&! average sea ice density
                                        ! Cox and Weeks, 1982: 919-974 kg/m^2
         rhow      = 1026.0_dbl_kind  ,&! density of seawater (kg/m^3)
         rhofresh  = 1000.0_dbl_kind    ! density of fresh water (kg/m^3)

!-----------------------------------------------------------------------
! Parameters for thermodynamics
!-----------------------------------------------------------------------

      real (kind=dbl_kind), public :: &
         hfrazilmin = 0.05_dbl_kind   ,&! min thickness of new frazil ice (m)
         cp_ice    = 2106._dbl_kind   ,&! specific heat of fresh ice (J/kg/K)
         cp_ocn    = 4218._dbl_kind   ,&! specific heat of ocn    (J/kg/K)
                                        ! freshwater value needed for enthalpy
         depressT  = 0.054_dbl_kind   ,&! Tf:brine salinity ratio (C/ppt)
         viscosity_dyn = 1.79e-3_dbl_kind, & ! dynamic viscosity of brine (kg/m/s)
         tscale_pnd_drain = c10       ,&! mushy macroscopic drainage timescale (days)
         Tocnfrz   = -1.8_dbl_kind    ,&! freezing temp of seawater (C),
                                        ! used as Tsfcn for open water
         Tffresh   = 273.15_dbl_kind  ,&! freezing temp of fresh ice (K)
         Lsub      = 2.835e6_dbl_kind ,&! latent heat, sublimation freshwater (J/kg)
         Lvap      = 2.501e6_dbl_kind ,&! latent heat, vaporization freshwater (J/kg)
         Timelt    = 0.0_dbl_kind     ,&! melting temperature, ice top surface  (C)
         Tsmelt    = 0.0_dbl_kind     ,&! melting temperature, snow top surface (C)
         ice_ref_salinity =4._dbl_kind,&! (ppt)
                                        ! kice is not used for mushy thermo
         kice      = 2.03_dbl_kind    ,&! thermal conductivity of fresh ice(W/m/deg)
         ksno      = 0.30_dbl_kind    ,&! thermal conductivity of snow  (W/m/deg)
         hs_min    = 1.e-4_dbl_kind   ,&! min snow thickness for computing zTsn (m)
         snowpatch = 0.02_dbl_kind    ,&! parameter for fractional snow area (m)
         saltmax   = 3.2_dbl_kind     ,&! max salinity at ice base for BL99 (ppt)
                                        ! phi_init, dSin0_frazil are for mushy thermo
         phi_init  = 0.75_dbl_kind    ,&! initial liquid fraction of frazil
         min_salin = p1               ,&! threshold for brine pocket treatment
         salt_loss = 0.4_dbl_kind     ,&! fraction of salt retained in zsalinity
         Tliquidus_max = c0           ,&! maximum liquidus temperature of mush (C)
         dSin0_frazil = c3            ,&! bulk salinity reduction of newly formed frazil
         dts_b     = 50._dbl_kind     ,&! zsalinity timestep
         ustar_min = 0.005_dbl_kind   ,&! minimum friction velocity for ocean heat flux (m/s)
         hi_min    = p01              ,&! minimum ice thickness allowed (m) for thermo
         ! mushy thermo
         a_rapid_mode      =  0.5e-3_dbl_kind,&! channel radius for rapid drainage mode (m)
         Rac_rapid_mode    =    10.0_dbl_kind,&! critical Rayleigh number
         aspect_rapid_mode =     1.0_dbl_kind,&! aspect ratio (larger is wider)
         dSdt_slow_mode    = -1.5e-7_dbl_kind,&! slow mode drainage strength (m s-1 K-1)
         phi_c_slow_mode   =    0.05_dbl_kind,&! critical liquid fraction porosity cutoff
         phi_i_mushy       =    0.85_dbl_kind  ! liquid fraction of congelation ice

      integer (kind=int_kind), public :: &
         ktherm = 1      ! type of thermodynamics
                         ! -1 none
                         ! 1 = Bitz and Lipscomb 1999
                         ! 2 = mushy layer theory

      character (char_len), public :: &
         conduct = 'bubbly', &          ! 'MU71' or 'bubbly'
         fbot_xfer_type = 'constant', & ! transfer coefficient type for ice-ocean heat flux
         cpl_frazil = 'fresh_ice_correction' ! type of coupling for frazil ice

      logical (kind=log_kind), public :: &
         calc_Tsfc     = .true. ,&! if true, calculate surface temperature
                                  ! if false, Tsfc is computed elsewhere and
                                  ! atmos-ice fluxes are provided to CICE
         update_ocn_f = .false. ,&! include fresh water and salt fluxes for frazil
         solve_zsal   = .false. ,&! if true, update salinity profile from solve_S_dt
         modal_aero   = .false. ,&! if true, use modal aerosal optical properties
                                  ! only for use with tr_aero or tr_zaero
         conserv_check = .false.  ! if true, do conservations checks and abort

      character(len=char_len), public :: &
         tfrz_option  = 'mushy'   ! form of ocean freezing temperature
                                  ! 'minus1p8' = -1.8 C
                                  ! 'constant' = Tocnfrz
                                  ! 'linear_salt' = -depressT * sss
                                  ! 'mushy' conforms with ktherm=2

      character(len=char_len), public :: &
         saltflux_option  = 'constant'! Salt flux computation
                                      ! 'constant' reference value of ice_ref_salinity
                                      ! 'prognostic' prognostic salt flux

!-----------------------------------------------------------------------
! Parameters for radiation
!-----------------------------------------------------------------------

      real (kind=dbl_kind), public :: &
         ! (Briegleb JGR 97 11475-11485  July 1992)
         emissivity = 0.985_dbl_kind,&! emissivity of snow and ice
         albocn     = 0.06_dbl_kind ,&! ocean albedo
         vonkar     = 0.4_dbl_kind  ,&! von Karman constant
         stefan_boltzmann = 567.0e-10_dbl_kind,&!  W/m^2/K^4
         ! (Ebert, Schramm and Curry JGR 100 15965-15975 Aug 1995)
         kappav     = 1.4_dbl_kind  ,&! vis extnctn coef in ice, wvlngth<700nm (1/m)
         hi_ssl     = 0.050_dbl_kind,&! ice surface scattering layer thickness (m)
         hs_ssl     = 0.040_dbl_kind,&! snow surface scattering layer thickness (m)
         ! baseline albedos for ccsm3 shortwave, set in namelist
         albicev    = 0.78_dbl_kind ,&! visible ice albedo for h > ahmax
         albicei    = 0.36_dbl_kind ,&! near-ir ice albedo for h > ahmax
         albsnowv   = 0.98_dbl_kind ,&! cold snow albedo, visible
         albsnowi   = 0.70_dbl_kind ,&! cold snow albedo, near IR
         ahmax      = 0.3_dbl_kind  ,&! thickness above which ice albedo is constant (m)
         ! dEdd tuning parameters, set in namelist
         R_ice      = c0   ,&! sea ice tuning parameter; +1 > 1sig increase in albedo
         R_pnd      = c0   ,&! ponded ice tuning parameter; +1 > 1sig increase in albedo
         R_snw      = c1p5 ,&! snow tuning parameter; +1 > ~.01 change in broadband albedo
         dT_mlt     = c1p5 ,&! change in temp for non-melt to melt snow grain
                             ! radius change (C)
         rsnw_mlt   = 1500._dbl_kind,&! maximum melting snow grain radius (10^-6 m)
         kalg       = 0.60_dbl_kind   ! algae absorption coefficient for 0.5 m thick layer
                                      ! 0.5 m path of 75 mg Chl a / m2
      ! weights for albedos
      ! 4 Jan 2007 BPB  Following are appropriate for complete cloud
      ! in a summer polar atmosphere with 1.5m bare sea ice surface:
      ! .636/.364 vis/nir with only 0.5% direct for each band.
      real (kind=dbl_kind), public :: &                 ! currently used only
         awtvdr = 0.00318_dbl_kind, &! visible, direct  ! for history and
         awtidr = 0.00182_dbl_kind, &! near IR, direct  ! diagnostics
         awtvdf = 0.63282_dbl_kind, &! visible, diffuse
         awtidf = 0.36218_dbl_kind   ! near IR, diffuse

      character (len=char_len), public :: &
         shortwave   = 'dEdd', & ! shortwave method, 'ccsm3' or 'dEdd' or 'dEdd_snicar_ad'
         albedo_type = 'ccsm3'   ! albedo parameterization, 'ccsm3' or 'constant'
                                 ! shortwave='dEdd' overrides this parameter

      ! Parameters for shortwave redistribution
      logical (kind=log_kind), public :: &
         sw_redist     = .false.

      real (kind=dbl_kind), public :: &
         sw_frac      = 0.9_dbl_kind    , & ! Fraction of internal shortwave moved to surface
         sw_dtemp     = 0.02_dbl_kind       ! temperature difference from melting

      ! Parameters for dEdd_snicar_ad
      character (len=char_len), public :: &
         snw_ssp_table = 'test'   ! lookup table: 'snicar' or 'test'

!-----------------------------------------------------------------------
! Parameters for dynamics, including ridging and strength
!-----------------------------------------------------------------------

      integer (kind=int_kind), public :: & ! defined in namelist
         kstrength   = 1, & ! 0 for simple Hibler (1979) formulation
                            ! 1 for Rothrock (1975) pressure formulation
         krdg_partic = 1, & ! 0 for Thorndike et al. (1975) formulation
                            ! 1 for exponential participation function
         krdg_redist = 1    ! 0 for Hibler (1980) formulation
                            ! 1 for exponential redistribution function

      real (kind=dbl_kind), public :: &
         Cf       = 17._dbl_kind     ,&! ratio of ridging work to PE change in ridging
         Pstar    = 2.75e4_dbl_kind  ,&! constant in Hibler strength formula
                                       ! (kstrength = 0)
         Cstar    = 20._dbl_kind     ,&! constant in Hibler strength formula
                                       ! (kstrength = 0)
         dragio   = 0.00536_dbl_kind ,&! ice-ocn drag coefficient
         thickness_ocn_layer1 = 2.0_dbl_kind,&! thickness of first ocean level (m)
         iceruf_ocn = 0.03_dbl_kind  ,&! under-ice roughness (m)
         gravit   = 9.80616_dbl_kind ,&! gravitational acceleration (m/s^2)
         mu_rdg = 3.0_dbl_kind ! e-folding scale of ridged ice, krdg_partic=1 (m^0.5)
                                       ! (krdg_redist = 1)

      logical (kind=log_kind), public :: &
         calc_dragio     = .false.     ! if true, calculate dragio from iceruf_ocn and thickness_ocn_layer1

!-----------------------------------------------------------------------
! Parameters for atmosphere
!-----------------------------------------------------------------------

      real (kind=dbl_kind), public :: &
         cp_air = 1005.0_dbl_kind    ,&! specific heat of air (J/kg/K)
         cp_wv  = 1.81e3_dbl_kind    ,&! specific heat of water vapor (J/kg/K)
         zvir   = 0.606_dbl_kind     ,&! rh2o/rair - 1.0
         zref   = 10._dbl_kind       ,&! reference height for stability (m)
         iceruf = 0.0005_dbl_kind    ,&! ice surface roughness (m)
         qqqice = 11637800._dbl_kind ,&! for qsat over ice
         TTTice = 5897.8_dbl_kind    ,&! for qsat over ice
         qqqocn = 627572.4_dbl_kind  ,&! for qsat over ocn
         TTTocn = 5107.4_dbl_kind    ,&! for qsat over ocn
         senscoef= 0.0012_dbl_kind   ,&! Sensible heat flux coefficient for constant-based boundary layer
         latncoef= 0.0015_dbl_kind     ! Latent heat flux coefficient for constant-based boundary layer

      character (len=char_len), public :: &
         atmbndy = 'similarity'        ! atmo boundary method, 'similarity', 'constant' or 'mixed'

      logical (kind=log_kind), public :: &
         calc_strair     = .true.  , & ! if true, calculate wind stress
         formdrag        = .false. , & ! if true, calculate form drag
         highfreq        = .false.     ! if true, calculate high frequency coupling

      integer (kind=int_kind), public :: &
         natmiter        = 5 ! number of iterations for atm boundary layer calcs

      ! Flux convergence tolerance
      real (kind=dbl_kind), public :: atmiter_conv = c0

!-----------------------------------------------------------------------
! Parameters for the ice thickness distribution
!-----------------------------------------------------------------------

      integer (kind=int_kind), public :: &
         kitd      = 1 ,&! type of itd conversions
                         !   0 = delta function
                         !   1 = linear remap
         kcatbound = 1   !   0 = old category boundary formula
                         !   1 = new formula giving round numbers
                         !   2 = WMO standard
                         !   3 = asymptotic formula

!-----------------------------------------------------------------------
! Parameters for the floe size distribution
!-----------------------------------------------------------------------

      integer (kind=int_kind), public :: &
         nfreq = 25                   ! number of frequencies

      real (kind=dbl_kind), public :: &
         floeshape = 0.66_dbl_kind    ! constant from Steele (unitless)

      real (kind=dbl_kind), public :: &
         floediam  = 300.0_dbl_kind   ! effective floe diameter for lateral melt (m)

      logical (kind=log_kind), public :: &
         wave_spec = .false.          ! if true, use wave forcing

      character (len=char_len), public :: &
         wave_spec_type = 'constant'  ! 'none', 'constant', or 'random'

!-----------------------------------------------------------------------
! Parameters for melt ponds
!-----------------------------------------------------------------------

      real (kind=dbl_kind), public :: &
         hs0       = 0.03_dbl_kind    ! snow depth for transition to bare sea ice (m)

      ! level-ice ponds
      character (len=char_len), public :: &
         frzpnd    = 'cesm'           ! pond refreezing parameterization

      real (kind=dbl_kind), public :: &
         dpscale   = 0.001_dbl_kind,& ! alter e-folding time scale for flushing (ktherm=1)
         rfracmin  = 0.15_dbl_kind, & ! minimum retained fraction of meltwater
         rfracmax  = 0.85_dbl_kind, & ! maximum retained fraction of meltwater
         pndaspect = 0.8_dbl_kind, &  ! ratio of pond depth to area fraction
         hs1       = 0.03_dbl_kind    ! snow depth for transition to bare pond ice (m)

      ! topo ponds
      real (kind=dbl_kind), public :: &
         hp1       = 0.01_dbl_kind    ! critical pond lid thickness for topo ponds

!-----------------------------------------------------------------------
! Parameters for snow redistribution, metamorphosis
!-----------------------------------------------------------------------

      character (len=char_len), public :: &
         snwredist       = 'none', &     ! type of snow redistribution
         snw_aging_table = 'test'        ! lookup table: 'snicar' or 'test' or 'file'

      logical (kind=log_kind), public :: &
         use_smliq_pnd = .false.     , & ! use liquid in snow for ponds
         snwgrain      = .false.         ! snow metamorphosis

      real (kind=dbl_kind), public :: &
         rsnw_fall  = 54.526_dbl_kind, & ! radius of new snow (10^-6 m)
         rsnw_tmax  = 1500.0_dbl_kind, & ! maximum snow radius (10^-6 m)
         rhosnew    =  100.0_dbl_kind, & ! new snow density (kg/m^3)
         rhosmin    =  100.0_dbl_kind, & ! minimum snow density (kg/m^3)
         rhosmax    =  450.0_dbl_kind, & ! maximum snow density (kg/m^3)
         windmin    =   10.0_dbl_kind, & ! minimum wind speed to compact snow (m/s)
         drhosdwind =   27.3_dbl_kind, & ! wind compaction factor for snow (kg s/m^4)
         snwlvlfac  =    0.3_dbl_kind    ! fractional increase in snow
                                         ! depth for bulk redistribution
      ! indices for aging lookup table
      integer (kind=int_kind), public :: &
         isnw_T,    & ! maximum temperature index
         isnw_Tgrd, & ! maximum temperature gradient index
         isnw_rhos    ! maximum snow density index

      ! dry snow aging parameters
      real (kind=dbl_kind), dimension(:), allocatable, public :: &
         snowage_rhos,  & ! snowage table dimension data for rhos (kg/m^3)
         snowage_Tgrd,  & ! snowage table dimension data for temp gradient (deg K/m)
         snowage_T        ! snowage table dimension data for temperature (deg K)
      real (kind=dbl_kind), dimension(:,:,:), allocatable, public :: &
         snowage_tau,   & ! snowage table 3D data for tau (10^-6 m)
         snowage_kappa, & ! snowage table 3D data for kappa (10^-6 m)
         snowage_drdt0    ! snowage table 3D data for drdt0 (10^-6 m/hr)

!-----------------------------------------------------------------------
! Parameters for biogeochemistry
!-----------------------------------------------------------------------

      character(char_len), public :: &
      ! skl biology parameters
         bgc_flux_type = 'Jin2006'  ! type of ocean-ice piston velocity (or 'constant')

      logical (kind=log_kind), public :: &
         z_tracers  = .false.,    & ! if .true., bgc or aerosol tracers are vertically resolved
         scale_bgc  = .false.,    & ! if .true., initialize bgc tracers proportionally with salinity
         solve_zbgc = .false.,    & ! if .true., solve vertical biochemistry portion of code
         dEdd_algae = .false.,    & ! if .true., algal absorption of shortwave is computed in the
         skl_bgc    = .false.       ! if true, solve skeletal biochemistry

      real (kind=dbl_kind), public :: &
         phi_snow     = p5              , & ! snow porosity
         grid_o       = c5              , & ! for bottom flux
         initbio_frac = c1              , & ! fraction of ocean trcr concentration in bio trcrs
         l_sk         = 7.0_dbl_kind    , & ! characteristic diffusive scale (m)
         grid_oS      = c5              , & ! for bottom flux
         l_skS        = 7.0_dbl_kind    , & ! characteristic skeletal layer thickness (m) (zsalinity)
         algal_vel    = 1.11e-8_dbl_kind, & ! 0.5 cm/d(m/s) Lavoie 2005  1.5 cm/day
         R_dFe2dust   = 0.035_dbl_kind  , & !  g/g (3.5% content) Tagliabue 2009
         dustFe_sol   = 0.005_dbl_kind  , & ! solubility fraction
         frazil_scav  = c1              , & ! fraction or multiple of bgc concentrated in frazil ice
         sk_l         = 0.03_dbl_kind   , & ! skeletal layer thickness (m)
         min_bgc      = 0.01_dbl_kind   , & ! fraction of ocean bgc concentration in surface melt
         T_max        = c0              , & ! maximum temperature (C)
         fsal         = c1              , & ! Salinity limitation (1)
         op_dep_min   = p1              , & ! light attenuates for optical depths exceeding min
         fr_graze_s   = p5              , & ! fraction of grazing spilled or slopped
         fr_graze_e   = p5              , & ! fraction of assimilation excreted
         fr_mort2min  = p5              , & ! fractionation of mortality to Am
         fr_dFe       = 0.3_dbl_kind    , & ! fraction of remineralized nitrogen
                                            ! (in units of algal iron)
         k_nitrif     = c0              , & ! nitrification rate (1/day)
         t_iron_conv  = 3065.0_dbl_kind , & ! desorption loss pFe to dFe (day)
         max_loss     = 0.9_dbl_kind    , & ! restrict uptake to % of remaining value
         max_dfe_doc1 = 0.2_dbl_kind    , & ! max ratio of dFe to saccharides in the ice
                                            ! (nM Fe/muM C)
         fr_resp      = 0.05_dbl_kind   , & ! fraction of algal growth lost due to respiration
         fr_resp_s    = 0.75_dbl_kind   , & ! DMSPd fraction of respiration loss as DMSPd
         y_sk_DMS     = p5              , & ! fraction conversion given high yield
         t_sk_conv    = 3.0_dbl_kind    , & ! Stefels conversion time (d)
         t_sk_ox      = 10.0_dbl_kind       ! DMS oxidation time (d)

!=======================================================================

      contains

!=======================================================================


!=======================================================================


    end module icepack_parameters

!=======================================================================
