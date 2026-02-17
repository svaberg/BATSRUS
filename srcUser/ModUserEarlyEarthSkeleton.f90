module ModUser

  use ModUserEmpty, ONLY: &
       user_set_boundary_cells, user_set_face_boundary, &
       user_set_cell_boundary, user_initial_perturbation, &
       user_set_ics, user_init_session, &
       user_specify_region, user_amr_criteria, user_get_log_var, &
       user_calc_sources_expl, user_calc_sources_impl, &
       user_init_point_implicit, user_get_b0, user_update_states, &
       user_calc_timestep, user_normalization, user_io_units, &
       user_set_resistivity, user_material_properties, i_type_block_user

  include 'user_module.h'

  character(len=*), parameter :: NameUserFile = 'srcUser/ModUserEarlyEarthSkeleton.f90'
  character(len=*), parameter :: NameUserModule = 'EARLY EARTH SKELETON (EMPTY HOOKS)'

  real :: XuvSigmaSi = 1.0e-22
  real :: XuvNeutralMinSi = 0.0
  integer :: XuvOctreeMaxIter = 12
  logical :: UseXuvOctree = .true.
  logical :: UseXuvBodyShadow = .true.
  real :: XuvFluxSi = 4.0
  logical :: UseTauDebugLog = .true.
  integer :: nTauUpdateCall = 0
  integer :: nTauUpdateCompute = 0

  logical :: IsTauOctreeReady = .false.
  integer :: nStepTauOctree = -1
  real, allocatable :: TauOctree_GB(:,:,:,:)
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
       case('#EARLYEARTHTAUXUV')
          call read_var('XuvSigmaSi', XuvSigmaSi)
          call read_var('XuvNeutralMinSi', XuvNeutralMinSi)
          call read_var('XuvFluxSi', XuvFluxSi)
       case('#EARLYEARTHXUVSWITCH')
          call read_var('UseXuvOctree', UseXuvOctree)
          call read_var('UseXuvBodyShadow', UseXuvBodyShadow)
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
       call update_tauxuv_octree
    case('timestep done')
       call update_tauxuv_octree
    case('load balance done')
       IsTauOctreeReady = .false.
       nStepTauOctree = -1
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
    case('tauxuv')
       if(allocated(TauOctree_GB)) then
          PlotVar_G = TauOctree_GB(:,:,:,iBlock)
       else
          PlotVar_G = 0.0
       end if
       PlotVarBody = 0.0
       UsePlotVarBody = .false.
       NameTecVar = 'TAUXUV'
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

  subroutine update_tauxuv_octree

    use BATL_lib, ONLY: nI, nJ, nK, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
         MaxBlock, nBlock, Unused_B, Xyz_DGB, CellSize_DB, CoordMax_D, &
         message_pass_cell, iProc
    use ModMain, ONLY: nStep, UseBody
    use ModPhysics, ONLY: No2Si_V, UnitX_, UnitRho_, rBody
    use ModAdvance, ONLY: State_VGB
    use ModVarIndexes, ONLY: Rho_
    use ModConst, ONLY: cProtonMass
    use ModIO, ONLY: write_prefix, iUnitOut

    integer :: iBlock, i, j, k, iIter
    real :: XSi, YSi, ZSi, DxSi, TauUp, TauNew, TauOld, TauDiffMax, NeutralProxySi
    real :: XFaceRightSi, SigmaNdx, FluxInSi, AbsFrac
    real :: RBodySi, RPerp2Si, XOCCSi
    logical :: IsShadowed
    real, parameter :: Small = 1.0e-30
    real, parameter :: TauMinOut = 1.0e-80
    real, parameter :: TauTol = 1.0e-12

    nTauUpdateCall = nTauUpdateCall + 1

    if(IsTauOctreeReady .and. nStepTauOctree == nStep) then
       RETURN
    end if

    nTauUpdateCompute = nTauUpdateCompute + 1
    if(UseTauDebugLog .and. iProc == 0) then
       call write_prefix
       write(iUnitOut,'(a,i8,a,i8,a,i8)') &
            'EARLYEARTH DEBUG: entering update_tauxuv_octree at nStep=', &
            nStep, ', call=', nTauUpdateCall, ', compute=', nTauUpdateCompute
    end if

    if(.not.allocated(TauOctree_GB)) then
      allocate(TauOctree_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    end if
    if(.not.allocated(XuvHeat_GB)) then
      allocate(XuvHeat_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    end if

    if(.not.UseXuvOctree) then
       TauOctree_GB = 0.0
       XuvHeat_GB = 0.0
       IsTauOctreeReady = .true.
       nStepTauOctree = nStep
       if(UseTauDebugLog .and. iProc == 0) then
          call write_prefix
          write(iUnitOut,'(a,i8)') &
               'EARLYEARTH DEBUG: update_tauxuv_octree disabled, nStep=', nStep
       end if
       RETURN
    end if

    TauOctree_GB = 0.0
    XuvHeat_GB = 0.0

    do iIter = 1, max(1, XuvOctreeMaxIter)
       call message_pass_cell(TauOctree_GB, nWidthIn=1, nProlongOrderIn=1)
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
                      TauUp = max(0.0, TauOctree_GB(i+1,j,k,iBlock))
                   end if

                   TauOld = TauOctree_GB(i,j,k,iBlock)
                   NeutralProxySi = max(0.0, &
                        State_VGB(Rho_,i,j,k,iBlock)*No2Si_V(UnitRho_)/cProtonMass)
                   NeutralProxySi = max(NeutralProxySi, XuvNeutralMinSi)
                   SigmaNdx = XuvSigmaSi*NeutralProxySi*DxSi
                   TauNew = TauUp + SigmaNdx
                   if(abs(TauNew) < TauMinOut) TauNew = 0.0
                   TauOctree_GB(i,j,k,iBlock) = TauNew
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

    IsTauOctreeReady = .true.
    nStepTauOctree = nStep

    if(UseTauDebugLog .and. iProc == 0) then
       call write_prefix
       write(iUnitOut,'(a,i8,a,1pe12.4,a,i8,a,i8)') &
            'EARLYEARTH DEBUG: exiting update_tauxuv_octree at nStep=', &
            nStep, ', TauDiffMax=', TauDiffMax, &
            ', call=', nTauUpdateCall, ', compute=', nTauUpdateCompute
    end if

  end subroutine update_tauxuv_octree

end module ModUser
