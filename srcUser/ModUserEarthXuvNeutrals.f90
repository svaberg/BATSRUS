module ModUser

  use BATL_lib, ONLY: test_start, test_stop

  use ModUserEmpty, ONLY: &
       user_set_boundary_cells, user_set_face_boundary, &
       user_set_cell_boundary, user_initial_perturbation, &
       user_set_ics_empty => user_set_ics, &
       user_specify_region, user_amr_criteria, &
       user_calc_sources_impl, &
       user_init_point_implicit, user_get_b0, user_update_states, &
       user_calc_timestep, user_normalization, user_io_units, &
       user_set_resistivity, user_material_properties, i_type_block_user

  include 'user_module.h'

  character(len=*), parameter :: NameUserFile = &
       'srcUser/ModUserEarthXuvNeutrals.f90'
  character(len=*), parameter :: NameUserModule = &
       'EARTH XUV NEUTRALS (MULTISPECIES Hp+H)'

  ! Physics scope for this module (current status, single-fluid only):
  ! Implemented:
  ! 1. XUV optical-depth sweep from +x and volumetric XUV heating diagnostic.
  ! 2. Heating source added to bulk MHD pressure/energy update.
  ! 3. Diagnostics from species state (Hp,H): fneutral, nH, nHp, ne.
  ! 4. Optional user IC profile for Hp/H from planet-centered radius.
  ! 5. Hp/H boundary split is synchronized from the same user neutral fractions.
  ! 6. Optional local H<->Hp chemistry source (photoionization + recombination).
  ! 7. Optional per-cell source limiter that preserves Hp/H positivity and
  !    chemistry-source mass conservation.
  ! 8. Optional fixed explicit chemistry subcycling for improved stiff handling.
  !
  ! Not implemented (roughly easiest -> hardest):
  ! 1. Implicit/adaptive chemistry integration controls beyond fixed subcycling.
  ! 2. Verification suite for limiting cases (no-XUV, strong-XUV, near-pure-ion/neutral starts).
  ! 3. Charge-exchange source terms and associated momentum/energy coupling.
  ! 4. Electron thermodynamics/chemistry beyond ne=nHp diagnostics.
  ! 5. Ion-fraction-dependent electromagnetic coupling model
  !    (including any explicit ion-neutral drag/slip representation).

  ! XUV controls
  logical :: UseXuv = .false.
  logical :: UseXuvBodyShadow = .true.
  logical :: UseHeatingSource = .true.
  real    :: XuvFluxSi = 4.64e-3
  real    :: XuvSigmaSi = 1.0e-22
  integer :: XuvMaxIter = 12

  ! Species handling:
  ! Use species state (Hp,H) everywhere. Optional one-time initialization from
  ! total rho is useful for bootstrapping from MHD-like starts.
  logical :: UseInitSpeciesFromRho = .false.
  real    :: NeutralFracInit = 1.0
  real    :: NeutralFracWind = 0.10
  real    :: NeutralFracPlanet = 0.30
  real    :: NeutralPlanetScale = 40.0
  real    :: NeutralFracFloor = 1.0e-6
  real    :: NeutralFracMax = 1.0

  ! Optional local chemistry (single-fluid species conversion only).
  logical :: UseNeutralChem = .false.
  logical :: UseChemTauAttenuation = .true.
  logical :: UseChemEnergyBookkeeping = .false.
  logical :: UseChemSourceLimiter = .true.
  logical :: UseChemSubcycling = .false.
  real    :: PhotoIonCoeffSi = 1.0e-7
  real    :: RecombCoeffSi = 2.6e-19
  real    :: ChemEnergyPerIonSi = 2.179872e-18
  integer :: nChemSubstep = 4

  logical :: IsXuvNeutralsReady = .false.
  integer :: nStepXuvNeutrals = -1

  real, allocatable :: TauXuv_GB(:,:,:,:)
  real, allocatable :: XuvHeat_GB(:,:,:,:)
  real, allocatable :: NeutralFrac_GB(:,:,:,:)
  real, allocatable :: Nh_GB(:,:,:,:)
  real, allocatable :: Nhp_GB(:,:,:,:)
  real, allocatable :: Ne_GB(:,:,:,:)

contains

  subroutine user_read_inputs

    use BATL_lib, ONLY: iProc, lVerbose
    use ModReadParam, ONLY: read_line, read_command, read_var
    use ModIO, ONLY: write_prefix, iUnitOut

    character(len=100) :: NameCommand

    do
       if(.not.read_line()) EXIT
       if(.not.read_command(NameCommand)) CYCLE

       select case(NameCommand)
       case('#EARTHXUV')
          call read_var('UseXuv', UseXuv)
          if(UseXuv) then
             call read_var('XuvFluxSi', XuvFluxSi)
             call read_var('XuvSigmaSi', XuvSigmaSi)
          end if
       case('#EARTHXUVSWITCH')
          call read_var('UseXuvBodyShadow', UseXuvBodyShadow)
          call read_var('UseHeatingSource', UseHeatingSource)
       case('#EARTHXUVNEUTRALSWITCH')
          call read_var('UseInitSpeciesFromRho', UseInitSpeciesFromRho)
       case('#EARTHXUVNEUTRALCHEMSWITCH')
          call read_var('UseNeutralChem', UseNeutralChem)
          call read_var('UseChemTauAttenuation', UseChemTauAttenuation)
          call read_var('UseChemEnergyBookkeeping', UseChemEnergyBookkeeping)
          call read_var('UseChemSourceLimiter', UseChemSourceLimiter)
          call read_var('UseChemSubcycling', UseChemSubcycling)
       case('#EARTHXUVNEUTRALS')
          call read_var('NeutralFracWind', NeutralFracWind)
          call read_var('NeutralFracPlanet', NeutralFracPlanet)
          call read_var('NeutralPlanetScale', NeutralPlanetScale)
          call read_var('NeutralFracFloor', NeutralFracFloor)
          call read_var('NeutralFracMax', NeutralFracMax)
       case('#EARTHXUVNEUTRALCHEM')
          call read_var('PhotoIonCoeffSi', PhotoIonCoeffSi)
          call read_var('RecombCoeffSi', RecombCoeffSi)
          call read_var('ChemEnergyPerIonSi', ChemEnergyPerIonSi)
          call read_var('nChemSubstep', nChemSubstep)
       case('#USERINPUTEND')
          if(iProc == 0 .and. lVerbose > 0) then
             call write_prefix
             write(iUnitOut,*) 'User read_input EARTHXUVNEUTRALS ends'
          end if
          EXIT
       case default
          cycle
       end select
    end do

    NeutralFracFloor = min(max(NeutralFracFloor, 0.0), 1.0)
    NeutralFracMax = min(max(NeutralFracMax, NeutralFracFloor), 1.0)
    NeutralFracWind = min(max(NeutralFracWind, NeutralFracFloor), NeutralFracMax)
    NeutralFracPlanet = min(max(NeutralFracPlanet, NeutralFracFloor), NeutralFracMax)
    NeutralPlanetScale = max(NeutralPlanetScale, 1.0e-6)
    NeutralFracInit = NeutralFracWind
    PhotoIonCoeffSi = max(PhotoIonCoeffSi, 0.0)
    RecombCoeffSi = max(RecombCoeffSi, 0.0)
    ChemEnergyPerIonSi = max(ChemEnergyPerIonSi, 0.0)
    nChemSubstep = max(1, nChemSubstep)

    call sync_neutral_fraction_controls

  end subroutine user_read_inputs

  subroutine user_init_session

    ! Keep BC fractions synchronized for each session/restart.
    call sync_neutral_fraction_controls

  end subroutine user_init_session

  subroutine sync_neutral_fraction_controls

    use BATL_lib, ONLY: iProc, lVerbose
    use ModAdvance, ONLY: nSpecies
    use ModPhysics, ONLY: LowDensityRatio, BodyNDim_I, BodyNSpeciesDim_I
    use ModIO, ONLY: write_prefix, iUnitOut

    real :: FracWind, FracPlanet, BodyTotalN

    FracWind = min(max(NeutralFracWind, NeutralFracFloor), NeutralFracMax)
    FracPlanet = min(max(NeutralFracPlanet, NeutralFracFloor), NeutralFracMax)

    ! Outer species split fallback follows the wind neutral fraction.
    LowDensityRatio = FracWind

    if(nSpecies < 2) RETURN

    ! Inner body species BC are always derived from the same neutral fraction
    ! used by the IC profile, so manual PARAM species values cannot drift.
    BodyTotalN = max(0.0, BodyNDim_I(1))
    if(BodyTotalN <= 0.0) then
       BodyTotalN = max(0.0, BodyNSpeciesDim_I(1)) + max(0.0, BodyNSpeciesDim_I(2))
    end if
    if(BodyTotalN <= 0.0) RETURN

    BodyNSpeciesDim_I = 0.0
    BodyNSpeciesDim_I(1) = (1.0 - FracPlanet)*BodyTotalN
    BodyNSpeciesDim_I(2) = FracPlanet*BodyTotalN
    BodyNDim_I(1) = BodyNSpeciesDim_I(1) + BodyNSpeciesDim_I(2)

    if(iProc == 0 .and. lVerbose > 0) then
      call write_prefix
      write(iUnitOut,*) 'EARTHXUVNEUTRALS BC sync: FracWind,FracPlanet=', &
           FracWind, FracPlanet
      call write_prefix
      write(iUnitOut,*) 'EARTHXUVNEUTRALS BC sync: BodyNDim Hp,H (/ccm)=', &
           BodyNSpeciesDim_I(1), BodyNSpeciesDim_I(2)
    end if

  end subroutine sync_neutral_fraction_controls


  subroutine user_action(NameAction)

    character(len=*), intent(in) :: NameAction

    select case(trim(NameAction))
    case('initial condition done')
       call timing_start('xuvneu_user_action')
       call update_xuv_neutrals
       call timing_stop('xuvneu_user_action')
    case('timestep done')
       call timing_start('xuvneu_user_action')
       call update_xuv_neutrals
       call timing_stop('xuvneu_user_action')
    case('load balance done')
       IsXuvNeutralsReady = .false.
       nStepXuvNeutrals = -1
    case default
       ! No action needed.
    end select

  end subroutine user_action

  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use BATL_lib, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK

    integer,          intent(in)    :: iBlock
    character(len=*), intent(in)    :: NameVar
    logical,          intent(in)    :: IsDimensional
    real,             intent(inout) :: PlotVar_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    real,             intent(out)   :: PlotVarBody
    logical,          intent(out)   :: UsePlotVarBody
    character(len=*), intent(inout) :: NameTecVar
    character(len=*), intent(inout) :: NameTecUnit
    character(len=*), intent(inout) :: NameIdlUnit
    logical,          intent(out)   :: IsFound

    select case(NameVar)
    case('xuvtau','tauxuv')
       if(allocated(TauXuv_GB)) then
          PlotVar_G = TauXuv_GB(:,:,:,iBlock)
       else
          PlotVar_G = 0.0
       end if
       PlotVarBody = 0.0
       UsePlotVarBody = .false.
       NameTecVar = 'XUVTAU'
       NameTecUnit = '[none]'
       NameIdlUnit = 'none'
       IsFound = .true.
    case('xuvheat')
       if(allocated(XuvHeat_GB)) then
          PlotVar_G = XuvHeat_GB(:,:,:,iBlock)
       else
          PlotVar_G = 0.0
       end if
       PlotVarBody = 0.0
       UsePlotVarBody = .false.
       NameTecVar = 'XUVHEAT'
       NameTecUnit = '[W/m^3]'
       NameIdlUnit = 'none'
       IsFound = .true.
    case('fneutral')
       if(allocated(NeutralFrac_GB)) then
          PlotVar_G = NeutralFrac_GB(:,:,:,iBlock)
       else
          PlotVar_G = NeutralFracInit
       end if
       PlotVarBody = 0.0
       UsePlotVarBody = .false.
       NameTecVar = 'FNEUTRAL'
       NameTecUnit = '[none]'
       NameIdlUnit = 'none'
       IsFound = .true.
    case('nh')
       if(allocated(Nh_GB)) then
          PlotVar_G = Nh_GB(:,:,:,iBlock)
       else
          PlotVar_G = 0.0
       end if
       PlotVarBody = 0.0
       UsePlotVarBody = .false.
       NameTecVar = 'NH'
       NameTecUnit = '[1/m^3]'
       NameIdlUnit = 'none'
       IsFound = .true.
    case('nhp')
       if(allocated(Nhp_GB)) then
          PlotVar_G = Nhp_GB(:,:,:,iBlock)
       else
          PlotVar_G = 0.0
       end if
       PlotVarBody = 0.0
       UsePlotVarBody = .false.
       NameTecVar = 'NHP'
       NameTecUnit = '[1/m^3]'
       NameIdlUnit = 'none'
       IsFound = .true.
    case('ne')
       if(allocated(Ne_GB)) then
          PlotVar_G = Ne_GB(:,:,:,iBlock)
       else
          PlotVar_G = 0.0
       end if
       PlotVarBody = 0.0
       UsePlotVarBody = .false.
       NameTecVar = 'NE'
       NameTecUnit = '[1/m^3]'
       NameIdlUnit = 'none'
       IsFound = .true.
    case default
       IsFound = .false.
       PlotVarBody = 0.0
       UsePlotVarBody = .false.
    end select

  end subroutine user_set_plot_var

  subroutine user_get_log_var(VarValue, NameVar, Radius)

    real, intent(out)            :: VarValue
    character(len=*), intent(in) :: NameVar
    real, intent(in), optional   :: Radius

    select case(trim(NameVar))
    case('usexuv')
       VarValue = merge(1.0, 0.0, UseXuv)
    case('xuvfluxsi')
       VarValue = XuvFluxSi
    case('xuvsigmasi')
       VarValue = XuvSigmaSi
    case('usexuvbshd')
       VarValue = merge(1.0, 0.0, UseXuvBodyShadow)
    case('useheatsrc')
       VarValue = merge(1.0, 0.0, UseHeatingSource)
    case('useinitspeciesrho','initsprho','useinitspe')
       VarValue = merge(1.0, 0.0, UseInitSpeciesFromRho)
    case('neutralfracinit','nfrinit')
       VarValue = NeutralFracInit
    case('neutralfracwind','nfrwind')
       VarValue = NeutralFracWind
    case('neutralfracplanet','nfrplan')
       VarValue = NeutralFracPlanet
    case('neutralplanetscale','nplscale','neutralpla')
       VarValue = NeutralPlanetScale
    case('neutralfracfloor','nfrfloor')
       VarValue = NeutralFracFloor
    case('useneutralchem','usechem','useneutral')
       VarValue = merge(1.0, 0.0, UseNeutralChem)
    case('usechemtauatt','chmtauatt','usechemtau')
       VarValue = merge(1.0, 0.0, UseChemTauAttenuation)
    case('usechemenergy','chmenergy','usechemene')
       VarValue = merge(1.0, 0.0, UseChemEnergyBookkeeping)
    case('usechemlimiter','chmlimit','usechemlim')
       VarValue = merge(1.0, 0.0, UseChemSourceLimiter)
    case('usechemsub','chemsub')
       VarValue = merge(1.0, 0.0, UseChemSubcycling)
    case('photoioncoeffsi','phioncoef','photoionco')
       VarValue = PhotoIonCoeffSi
    case('recombcoeffsi','recombcoef')
       VarValue = RecombCoeffSi
    case('chemenergyperionsi','chemener','chemenergy')
       VarValue = ChemEnergyPerIonSi
    case('nchemsubstep','nchemsub','nchemsubst')
       VarValue = real(nChemSubstep)
    case default
       VarValue = -7777.0
    end select

  end subroutine user_get_log_var

  subroutine user_calc_sources_expl(iBlock)

    use BATL_lib, ONLY: nI, nJ, nK, Used_GB, Xyz_DGB
    use ModAdvance, ONLY: Source_VC, State_VGB, nSpecies
    use ModConservative, ONLY: UseNonConservative, nConservCrit, IsConserv_CB
    use ModMain, ONLY: Dt, UseBody
    use ModPhysics, ONLY: Si2No_V, No2Si_V, UnitEnergyDens_, UnitT_, UnitRho_, &
         UnitX_, GammaMinus1, rBody
    use ModVarIndexes, ONLY: p_, Energy_, RhoHp_, RhoH_
    use ModConst, ONLY: cProtonMass
    use ModBatsrusUtility, ONLY: stop_mpi

    integer, intent(in) :: iBlock
    integer :: i, j, k, iSub, nSubLocal
    real :: HeatRateNo
    real :: DtNo, DtSi, DtSubSi, DnIonSi, DnRecombSi, DnNetSi
    real :: RhoHpSi, RhoHSi, NHSi, NHpSi, NESi
    real :: NHWorkSi, NHpWorkSi
    real :: RhoSourceNo, PhotoRateSi, ChemPowerSi, ChemPowerNo
    real :: RhoHpNo, RhoHNo, SourceBaseHpNo, SourceBaseHNo, SourceMinNo, SourceMaxNo
    real :: TauLocal, XSi, YSi, ZSi, RBodySi, RPerp2Si, XOCCSi
    logical :: IsShadowed
    logical :: DoTest
    character(len=*), parameter :: NameSub = 'user_calc_sources_expl'
    real, parameter :: Small = 1.0e-30

    call test_start(NameSub, DoTest)

    if(UseHeatingSource .and. UseXuv .and. allocated(XuvHeat_GB)) then
       call timing_start('xuvneu_heat_src')
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(.not.Used_GB(i,j,k,iBlock)) CYCLE

          ! Volumetric XUV heating Q [W/m^3] from Beer-Lambert attenuation.
          HeatRateNo = XuvHeat_GB(i,j,k,iBlock) * &
               Si2No_V(UnitEnergyDens_) / Si2No_V(UnitT_)

          ! Non-conservative cells advance pressure directly: dp/dt=(gamma-1)Q.
          if(UseNonConservative .and. &
               .not.(nConservCrit > 0 .and. IsConserv_CB(i,j,k,iBlock))) then
             Source_VC(p_,i,j,k) = Source_VC(p_,i,j,k) + GammaMinus1*HeatRateNo
          else
             Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) + HeatRateNo
          end if
       end do; end do; end do
       call timing_stop('xuvneu_heat_src')
    end if

    if(UseNeutralChem) then
       if(nSpecies < 2) call stop_mpi('UseNeutralChem requires nSpecies>=2 (Hp+H)')

       call timing_start('xuvneu_chem_src')
       DtNo = max(Dt, Small)
       DtSi = max(Dt*No2Si_V(UnitT_), Small)
       RBodySi = rBody*No2Si_V(UnitX_)

       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(.not.Used_GB(i,j,k,iBlock)) CYCLE

          RhoHpSi = max(0.0, State_VGB(RhoHp_,i,j,k,iBlock)*No2Si_V(UnitRho_))
          RhoHSi  = max(0.0, State_VGB(RhoH_, i,j,k,iBlock)*No2Si_V(UnitRho_))
          RhoHpNo = max(0.0, State_VGB(RhoHp_,i,j,k,iBlock))
          RhoHNo  = max(0.0, State_VGB(RhoH_, i,j,k,iBlock))
          NHpSi = RhoHpSi/max(cProtonMass, Small)
          NHSi  = RhoHSi /max(cProtonMass, Small)
          NESi  = NHpSi

          ! Local photoionization rate [1/s], optionally attenuated by XUV optical depth.
          PhotoRateSi = PhotoIonCoeffSi
          if(UseChemTauAttenuation .and. allocated(TauXuv_GB)) then
             TauLocal = max(0.0, TauXuv_GB(i,j,k,iBlock))
             PhotoRateSi = PhotoRateSi*exp(-TauLocal)

             if(UseXuvBodyShadow .and. UseBody .and. RBodySi > 0.0) then
                XSi = Xyz_DGB(1,i,j,k,iBlock)*No2Si_V(UnitX_)
                YSi = Xyz_DGB(2,i,j,k,iBlock)*No2Si_V(UnitX_)
                ZSi = Xyz_DGB(3,i,j,k,iBlock)*No2Si_V(UnitX_)
                IsShadowed = .false.
                RPerp2Si = YSi*YSi + ZSi*ZSi
                if(RPerp2Si < RBodySi*RBodySi) then
                   XOCCSi = sqrt(max(0.0, RBodySi*RBodySi - RPerp2Si))
                   IsShadowed = XSi <= XOCCSi
                end if
                if(IsShadowed) PhotoRateSi = 0.0
             end if
          end if

          ! Number-density conversion rates [1/m^3/s].
          if(UseChemSubcycling .and. nChemSubstep > 1) then
             nSubLocal = nChemSubstep
             DtSubSi = DtSi/real(nSubLocal)
             NHWorkSi = NHSi
             NHpWorkSi = NHpSi
             do iSub = 1, nSubLocal
                DnIonSi = PhotoRateSi*NHWorkSi
                DnRecombSi = RecombCoeffSi*NHpWorkSi*NHpWorkSi
                DnIonSi = min(DnIonSi, NHWorkSi/max(DtSubSi, Small))
                DnRecombSi = min(DnRecombSi, NHpWorkSi/max(DtSubSi, Small))
                NHWorkSi = max(0.0, NHWorkSi - (DnIonSi - DnRecombSi)*DtSubSi)
                NHpWorkSi = max(0.0, NHpWorkSi + (DnIonSi - DnRecombSi)*DtSubSi)
             end do
             DnNetSi = (NHpWorkSi - NHpSi)/DtSi
          else
             DnIonSi = PhotoRateSi*NHSi
             DnRecombSi = RecombCoeffSi*NHpSi*NESi
             DnIonSi = min(DnIonSi, NHSi /DtSi)
             DnRecombSi = min(DnRecombSi, NHpSi/DtSi)
             DnNetSi = DnIonSi - DnRecombSi
          end if

          ! Convert to normalized mass-density source for Hp/H species variables.
          RhoSourceNo = cProtonMass*DnNetSi*Si2No_V(UnitRho_)/Si2No_V(UnitT_)

          if(UseChemSourceLimiter) then
             ! Enforce positivity for species after explicit source update while
             ! preserving conservative chemistry transfer (+Hp, -H).
             SourceBaseHpNo = Source_VC(RhoHp_,i,j,k)
             SourceBaseHNo  = Source_VC(RhoH_, i,j,k)
             SourceMinNo = -SourceBaseHpNo - RhoHpNo/DtNo
             SourceMaxNo =  SourceBaseHNo  + RhoHNo /DtNo
             if(SourceMinNo > SourceMaxNo) SourceMinNo = SourceMaxNo
             RhoSourceNo = min(max(RhoSourceNo, SourceMinNo), SourceMaxNo)
          end if

          Source_VC(RhoHp_,i,j,k) = Source_VC(RhoHp_,i,j,k) + RhoSourceNo
          Source_VC(RhoH_, i,j,k) = Source_VC(RhoH_, i,j,k) - RhoSourceNo

          if(UseChemEnergyBookkeeping) then
             ! Optional bookkeeping: net ionization consumes internal energy,
             ! net recombination releases internal energy.
             DnNetSi = RhoSourceNo*Si2No_V(UnitT_) / &
                  (max(cProtonMass, Small)*Si2No_V(UnitRho_))
             ChemPowerSi = -ChemEnergyPerIonSi*DnNetSi
             ChemPowerNo = ChemPowerSi*Si2No_V(UnitEnergyDens_) / Si2No_V(UnitT_)
             if(UseNonConservative .and. &
                  .not.(nConservCrit > 0 .and. IsConserv_CB(i,j,k,iBlock))) then
                Source_VC(p_,i,j,k) = Source_VC(p_,i,j,k) + GammaMinus1*ChemPowerNo
             else
                Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) + ChemPowerNo
             end if
          end if
       end do; end do; end do

       call timing_stop('xuvneu_chem_src')
    end if

    call test_stop(NameSub, DoTest)

  end subroutine user_calc_sources_expl

  subroutine update_xuv_neutrals

    use BATL_lib, ONLY: nI, nJ, nK, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
         MaxBlock, nBlock, Unused_B, Xyz_DGB, CellSize_DB, CoordMax_D, &
         message_pass_cell
    use ModMain, ONLY: nStep, Dt, UseBody
    use ModPhysics, ONLY: No2Si_V, UnitX_, UnitRho_, UnitT_, rBody
    use ModAdvance, ONLY: State_VGB, nSpecies
    use ModBatsrusUtility, ONLY: stop_mpi
    use ModVarIndexes, ONLY: Rho_, RhoHp_, RhoH_
    use ModPhysics, ONLY: rBody
    use ModConst, ONLY: cProtonMass

    integer :: iBlock, i, j, k, iIter
    real :: XSi, YSi, ZSi, DxSi, TauUp, TauNew, TauOld
    real :: FluxInSi, AbsFrac, SigmaNdx, TauDiffMax, NeutralDiffMax
    real :: nTotSi, nHSi, nHpSi, nESi, fNeutralOld
    real :: RhoHpSi, RhoHSi
    real :: RBodySi, RPerp2Si, XOCCSi, XFaceRightSi
    logical :: IsShadowed
    logical :: DoTest
    character(len=*), parameter :: NameSub = 'update_xuv_neutrals'
    real, parameter :: Small = 1.0e-30
    real, parameter :: TauTol = 1.0e-12

    call test_start(NameSub, DoTest)

    call init_neutrals_arrays

    if(IsXuvNeutralsReady .and. nStepXuvNeutrals == nStep) then
       call test_stop(NameSub, DoTest)
       RETURN
    end if

    if(nSpecies < 2) then
       call stop_mpi('ModUserEarthXuvNeutrals requires Hp+H species state (nSpecies>=2)')
    end if

    if(.not.UseXuv) then
       TauXuv_GB = 0.0
       XuvHeat_GB = 0.0
    end if

    call timing_start('xuvneu_update')

    do iIter = 1, max(1, XuvMaxIter)
       if(UseXuv) call message_pass_cell(TauXuv_GB, nWidthIn=1, nProlongOrderIn=1)

       TauDiffMax = 0.0
       NeutralDiffMax = 0.0

       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE

          DxSi = CellSize_DB(1,iBlock)*No2Si_V(UnitX_)
          RBodySi = rBody*No2Si_V(UnitX_)

          do k = 1, nK
             do j = 1, nJ
                do i = nI, 1, -1

                   XSi = Xyz_DGB(1,i,j,k,iBlock)*No2Si_V(UnitX_)
                   YSi = Xyz_DGB(2,i,j,k,iBlock)*No2Si_V(UnitX_)
                   ZSi = Xyz_DGB(3,i,j,k,iBlock)*No2Si_V(UnitX_)

                   if(UseXuv) then
                      XFaceRightSi = XSi + 0.5*DxSi
                      if(XFaceRightSi >= CoordMax_D(1)*No2Si_V(UnitX_) - 10.0*Small) then
                         TauUp = 0.0
                      else
                         TauUp = max(0.0, TauXuv_GB(i+1,j,k,iBlock))
                      end if
                   else
                      TauUp = 0.0
                   end if

                   RhoHpSi = max(0.0, State_VGB(RhoHp_,i,j,k,iBlock)*No2Si_V(UnitRho_))
                   RhoHSi  = max(0.0, State_VGB(RhoH_, i,j,k,iBlock)*No2Si_V(UnitRho_))
                   nHpSi = RhoHpSi/max(cProtonMass, Small)
                   nHSi  = RhoHSi /max(cProtonMass, Small)
                   nTotSi = nHpSi + nHSi
                   if(nTotSi <= Small) then
                      fNeutralOld = min(max(NeutralFracWind, NeutralFracFloor), NeutralFracMax)
                   else
                      fNeutralOld = min(max(nHSi/nTotSi, NeutralFracFloor), NeutralFracMax)
                   end if

                   IsShadowed = .false.
                   if(UseXuvBodyShadow .and. UseBody .and. RBodySi > 0.0) then
                      RPerp2Si = YSi*YSi + ZSi*ZSi
                      if(RPerp2Si < RBodySi*RBodySi) then
                         XOCCSi = sqrt(max(0.0, RBodySi*RBodySi - RPerp2Si))
                         IsShadowed = XSi <= XOCCSi
                      end if
                   end if

                   if(UseXuv .and. .not.IsShadowed) then
                      FluxInSi = XuvFluxSi*exp(-TauUp)
                   else
                      FluxInSi = 0.0
                   end if

                   NeutralFrac_GB(i,j,k,iBlock) = fNeutralOld
                   NeutralDiffMax = max(NeutralDiffMax, 0.0)

                   nESi  = nHpSi
                   Nh_GB(i,j,k,iBlock) = nHSi
                   Nhp_GB(i,j,k,iBlock) = nHpSi
                   Ne_GB(i,j,k,iBlock) = nESi

                   if(UseXuv) then
                      TauOld = TauXuv_GB(i,j,k,iBlock)
                      SigmaNdx = XuvSigmaSi*nHSi*DxSi
                      TauNew = TauUp + SigmaNdx
                      TauXuv_GB(i,j,k,iBlock) = TauNew
                      TauDiffMax = max(TauDiffMax, abs(TauNew - TauOld))

                      if(IsShadowed) then
                         XuvHeat_GB(i,j,k,iBlock) = 0.0
                      else
                         AbsFrac = 1.0 - exp(-SigmaNdx)
                         XuvHeat_GB(i,j,k,iBlock) = FluxInSi*AbsFrac/max(DxSi, Small)
                      end if
                   else
                      TauXuv_GB(i,j,k,iBlock) = 0.0
                      XuvHeat_GB(i,j,k,iBlock) = 0.0
                   end if

                end do
             end do
          end do
       end do

       if(TauDiffMax < TauTol .and. NeutralDiffMax < TauTol) EXIT
    end do

    call timing_stop('xuvneu_update')

    IsXuvNeutralsReady = .true.
    nStepXuvNeutrals = nStep

    call test_stop(NameSub, DoTest)

  end subroutine update_xuv_neutrals


  subroutine user_set_ics(iBlock)

    use BATL_lib, ONLY: nI, nJ, nK, Used_GB, Xyz_DGB, test_start, test_stop
    use ModAdvance, ONLY: State_VGB, nSpecies
    use ModBatsrusUtility, ONLY: stop_mpi
    use ModVarIndexes, ONLY: Rho_, RhoHp_, RhoH_
    use ModPhysics, ONLY: rBody

    integer, intent(in) :: iBlock
    integer :: i, j, k
    real :: FracH, RhoTot, XLoc, YLoc, ZLoc, RLoc, T01
    real :: FracWind, FracPlanet, ScalePlanet, FracFloor, SmoothStep
    logical :: DoTest
    character(len=*), parameter :: NameSub = 'user_set_ics'

    call test_start(NameSub, DoTest, iBlock)

    if(.not.UseInitSpeciesFromRho) then
       call test_stop(NameSub, DoTest, iBlock)
       RETURN
    end if

    if(nSpecies < 2) call stop_mpi('UseInitSpeciesFromRho requires nSpecies>=2')

    FracFloor = min(max(NeutralFracFloor, 0.0), 1.0)
    FracWind = min(max(NeutralFracWind, FracFloor), NeutralFracMax)
    FracPlanet = min(max(NeutralFracPlanet, FracFloor), NeutralFracMax)
    ScalePlanet = max(NeutralPlanetScale, 1.0e-6)

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(.not.Used_GB(i,j,k,iBlock)) CYCLE

       ! Planet-centered radial profile: near-planet -> FracPlanet, far field -> FracWind.
       XLoc = Xyz_DGB(1,i,j,k,iBlock)
       YLoc = Xyz_DGB(2,i,j,k,iBlock)
       ZLoc = Xyz_DGB(3,i,j,k,iBlock)
       RLoc = sqrt(XLoc*XLoc + YLoc*YLoc + ZLoc*ZLoc)

       T01 = min(1.0, max(0.0, (RLoc - rBody)/ScalePlanet))
       SmoothStep = T01*T01*(3.0 - 2.0*T01)
       FracH = FracPlanet + (FracWind - FracPlanet)*SmoothStep
       FracH = min(max(FracH, FracFloor), NeutralFracMax)

       RhoTot = max(0.0, State_VGB(Rho_,i,j,k,iBlock))
       State_VGB(RhoH_, i,j,k,iBlock) = FracH*RhoTot
       State_VGB(RhoHp_,i,j,k,iBlock) = (1.0 - FracH)*RhoTot
    end do; end do; end do

    call test_stop(NameSub, DoTest, iBlock)

  end subroutine user_set_ics


  subroutine init_neutrals_arrays

    use BATL_lib, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK, MaxBlock

    if(.not.allocated(TauXuv_GB)) then
       allocate(TauXuv_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
       TauXuv_GB = 0.0
    end if
    if(.not.allocated(XuvHeat_GB)) then
       allocate(XuvHeat_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
       XuvHeat_GB = 0.0
    end if
    if(.not.allocated(NeutralFrac_GB)) then
       allocate(NeutralFrac_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
       NeutralFrac_GB = NeutralFracInit
    end if
    if(.not.allocated(Nh_GB)) then
       allocate(Nh_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
       Nh_GB = 0.0
    end if
    if(.not.allocated(Nhp_GB)) then
       allocate(Nhp_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
       Nhp_GB = 0.0
    end if
    if(.not.allocated(Ne_GB)) then
       allocate(Ne_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
       Ne_GB = 0.0
    end if

  end subroutine init_neutrals_arrays

end module ModUser
