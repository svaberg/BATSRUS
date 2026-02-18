module ModUser

  use BATL_lib, ONLY: test_start, test_stop

  use ModUserEmpty, ONLY: &
       user_set_boundary_cells, user_set_face_boundary, &
       user_set_cell_boundary, user_initial_perturbation, &
       user_set_ics, user_init_session, &
       user_specify_region, user_amr_criteria, &
       user_calc_sources_impl, &
       user_init_point_implicit, user_get_b0, user_update_states, &
       user_calc_timestep, user_normalization, user_io_units, &
       user_set_resistivity, user_material_properties, i_type_block_user

  include 'user_module.h'

  character(len=*), parameter :: NameUserFile = 'srcUser/ModUserEarlyEarthSkeleton.f90'
  character(len=*), parameter :: NameUserModule = 'EARLY EARTH SKELETON (EMPTY HOOKS)'

  real :: XuvSigmaSi = 1.0e-22
  integer :: XuvMaxIter = 12
  logical :: UseXuv = .false.
  logical :: UseXuvBodyShadow = .true.
  logical :: UseHeatingSource = .true.
  real :: XuvFluxSi = 4.64e-3

  logical :: IsTauXuvReady = .false.
  integer :: nStepTauXuv = -1
  real, allocatable :: TauXuv_GB(:,:,:,:)
  real, allocatable :: XuvHeat_GB(:,:,:,:)

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
       case('#USERINPUTEND')
          if(iProc == 0 .and. lVerbose > 0) then
             call write_prefix
             write(iUnitOut,*) 'User read_input EARLYEARTH ends'
          end if
          EXIT
       case default
          cycle
       end select
    end do

  end subroutine user_read_inputs

  subroutine user_action(NameAction)

    character(len=*), intent(in) :: NameAction

    select case(trim(NameAction))
    case('initial condition done')
       call timing_start('xuv_user_action')
       call update_tauxuv
       call timing_stop('xuv_user_action')
    case('timestep done')
       call timing_start('xuv_user_action')
       call update_tauxuv
       call timing_stop('xuv_user_action')
    case('load balance done')
       IsTauXuvReady = .false.
       nStepTauXuv = -1
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
    case default
       VarValue = -7777.0
    end select

  end subroutine user_get_log_var

  subroutine user_calc_sources_expl(iBlock)

    use BATL_lib, ONLY: nI, nJ, nK, Used_GB
    use ModAdvance, ONLY: Source_VC
    use ModConservative, ONLY: UseNonConservative, nConservCrit, IsConserv_CB
    use ModPhysics, ONLY: Si2No_V, UnitEnergyDens_, UnitT_, GammaMinus1
    use ModVarIndexes, ONLY: p_, Energy_

    integer, intent(in) :: iBlock
    integer :: i, j, k
    real :: HeatRateNo
    logical :: DoTest
    character(len=*), parameter :: NameSub = 'user_calc_sources_expl'

    call test_start(NameSub, DoTest)
    call timing_start('xuv_src_hook')

    if(UseHeatingSource .and. UseXuv .and. allocated(XuvHeat_GB)) then
       call timing_start('xuv_heat_src')
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(.not.Used_GB(i,j,k,iBlock)) CYCLE

          ! Volumetric XUV heating Q [W/m^3] from Beer-Lambert attenuation.
          HeatRateNo = XuvHeat_GB(i,j,k,iBlock) * &
               Si2No_V(UnitEnergyDens_) / Si2No_V(UnitT_)

          ! In non-conservative cells BATSRUS advances pressure directly:
          ! dp/dt = (gamma-1)Q. Conservative cells use dE/dt = Q.
          if(UseNonConservative .and. &
               .not.(nConservCrit > 0 .and. IsConserv_CB(i,j,k,iBlock))) then
             Source_VC(p_,i,j,k) = Source_VC(p_,i,j,k) + GammaMinus1*HeatRateNo
          else
             Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) + HeatRateNo
          end if
       end do; end do; end do
       call timing_stop('xuv_heat_src')
    end if

    call timing_stop('xuv_src_hook')
    call test_stop(NameSub, DoTest)

  end subroutine user_calc_sources_expl

  subroutine update_tauxuv

    use BATL_lib, ONLY: nI, nJ, nK, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
         MaxBlock, nBlock, Unused_B, Xyz_DGB, CellSize_DB, CoordMax_D, &
         message_pass_cell
    use ModMain, ONLY: nStep, UseBody
    use ModPhysics, ONLY: No2Si_V, UnitX_, UnitRho_, rBody
    use ModAdvance, ONLY: State_VGB
    use ModVarIndexes, ONLY: Rho_
    use ModConst, ONLY: cProtonMass

    integer :: iBlock, i, j, k, iIter
    real :: XSi, YSi, ZSi, DxSi, TauUp, TauNew, TauOld, TauDiffMax, NeutralProxySi
    real :: XFaceRightSi, SigmaNdx, FluxInSi, AbsFrac
    real :: RBodySi, RPerp2Si, XOCCSi
    logical :: IsShadowed
    logical :: DoTest
    character(len=*), parameter :: NameSub = 'update_tauxuv'
    real, parameter :: Small = 1.0e-30
    real, parameter :: TauMinOut = 1.0e-80
    real, parameter :: TauTol = 1.0e-12

    call test_start(NameSub, DoTest)

    if(IsTauXuvReady .and. nStepTauXuv == nStep) then
       call test_stop(NameSub, DoTest)
       RETURN
    end if

    if(.not.allocated(TauXuv_GB)) then
      allocate(TauXuv_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    end if
    if(.not.allocated(XuvHeat_GB)) then
      allocate(XuvHeat_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    end if

    if(.not.UseXuv) then
       TauXuv_GB = 0.0
       XuvHeat_GB = 0.0
       IsTauXuvReady = .true.
       nStepTauXuv = nStep
       call test_stop(NameSub, DoTest)
       RETURN
    end if

    call timing_start('xuv_tau_update')
    TauXuv_GB = 0.0
    XuvHeat_GB = 0.0

    do iIter = 1, max(1, XuvMaxIter)
       call message_pass_cell(TauXuv_GB, nWidthIn=1, nProlongOrderIn=1)
       TauDiffMax = 0.0

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
                   XFaceRightSi = XSi + 0.5*DxSi
                   if(XFaceRightSi >= CoordMax_D(1)*No2Si_V(UnitX_) - 10.0*Small) then
                      TauUp = 0.0
                   else
                      TauUp = max(0.0, TauXuv_GB(i+1,j,k,iBlock))
                   end if

                   TauOld = TauXuv_GB(i,j,k,iBlock)
                   NeutralProxySi = max(0.0, &
                        State_VGB(Rho_,i,j,k,iBlock)*No2Si_V(UnitRho_)/cProtonMass)
                   ! Optical depth increment for one cell:
                   ! dTau = sigma * n * ds, with n approximated by rho/mp.
                   SigmaNdx = XuvSigmaSi*NeutralProxySi*DxSi
                   TauNew = TauUp + SigmaNdx
                   if(abs(TauNew) < TauMinOut) TauNew = 0.0
                   TauXuv_GB(i,j,k,iBlock) = TauNew
                   TauDiffMax = max(TauDiffMax, abs(TauNew - TauOld))

                   IsShadowed = .false.
                   if(UseXuvBodyShadow .and. UseBody .and. RBodySi > 0.0) then
                      RPerp2Si = YSi*YSi + ZSi*ZSi
                      if(RPerp2Si < RBodySi*RBodySi) then
                         XOCCSi = sqrt(max(0.0, RBodySi*RBodySi - RPerp2Si))
                         IsShadowed = XSi <= XOCCSi
                      end if
                   end if

                   if(IsShadowed) then
                      XuvHeat_GB(i,j,k,iBlock) = 0.0
                   else
                      ! Beer-Lambert attenuation and absorbed power:
                      ! F_in = F0*exp(-Tau_up), AbsFrac = 1-exp(-dTau), Q = F_in*AbsFrac/ds.
                      FluxInSi = XuvFluxSi*exp(-TauUp)
                      AbsFrac = 1.0 - exp(-SigmaNdx)
                      XuvHeat_GB(i,j,k,iBlock) = FluxInSi*AbsFrac/max(DxSi, Small)
                      if(abs(XuvHeat_GB(i,j,k,iBlock)) < TauMinOut) &
                           XuvHeat_GB(i,j,k,iBlock) = 0.0
                   end if
                end do
             end do
          end do
       end do

       if(TauDiffMax < TauTol) EXIT
    end do

    call timing_stop('xuv_tau_update')
    IsTauXuvReady = .true.
    nStepTauXuv = nStep

    call test_stop(NameSub, DoTest)

  end subroutine update_tauxuv

end module ModUser
