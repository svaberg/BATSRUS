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

  real :: ConstantPlotValue = 1.0
  real :: XuvSigmaSi = 1.0e-22
  real :: XuvNeutralN0Si = 1.0e15
  real :: XuvNeutralScaleHeightSi = 5.0e5
  real :: XuvNeutralR0Si = 6.371e6
  real :: XuvRayStepSi = 1.0e5
  real :: XuvRayMaxSi = -1.0
  real :: XuvNeutralMinSi = 0.0
  integer :: XuvOctreeMaxIter = 12
  logical :: UseTauDebugLog = .true.
  integer :: nTauUpdateCall = 0
  integer :: nTauUpdateCompute = 0

  logical :: IsTauOctreeReady = .false.
  integer :: nStepTauOctree = -1
  real, allocatable :: TauOctree_GB(:,:,:,:)

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
       case('#EARLYEARTHPLOT')
          call read_var('ConstantPlotValue', ConstantPlotValue)
       case('#EARLYEARTHTAUXUV')
          call read_var('XuvSigmaSi', XuvSigmaSi)
          call read_var('XuvNeutralN0Si', XuvNeutralN0Si)
          call read_var('XuvNeutralScaleHeightSi', XuvNeutralScaleHeightSi)
          call read_var('XuvNeutralR0Si', XuvNeutralR0Si)
          call read_var('XuvRayStepSi', XuvRayStepSi)
          call read_var('XuvRayMaxSi', XuvRayMaxSi)
          call read_var('XuvNeutralMinSi', XuvNeutralMinSi)
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

    use BATL_lib, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK, Xyz_DGB, CoordMax_D
    use ModPhysics, ONLY: No2Si_V, UnitX_

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
    integer :: i, j, k, iStep, nStep
    real :: XStartSi, XEndSi, XSi, YSi, ZSi, RSquareSi, RSi
    real :: DS, S, Tau, NeutralSi, Weight
    real, parameter :: Small = 1.0e-30
    real, parameter :: TauMinOut = 1.0e-80

    select case(NameVar)
    case('uconst')
       PlotVar_G = ConstantPlotValue
       PlotVarBody = ConstantPlotValue
       UsePlotVarBody = .true.
       NameTecVar = 'UCONST'
       NameTecUnit = '[none]'
       NameIdlUnit = 'none'
       IsFound = .true.
    case('tauxuv')
       PlotVar_G = 0.0
       do k = MinK, MaxK
          do j = MinJ, MaxJ
             do i = MinI, MaxI
                XStartSi = Xyz_DGB(1,i,j,k,iBlock)*No2Si_V(UnitX_)
                YSi = Xyz_DGB(2,i,j,k,iBlock)*No2Si_V(UnitX_)
                ZSi = Xyz_DGB(3,i,j,k,iBlock)*No2Si_V(UnitX_)

                if(XuvRayMaxSi > 0.0) then
                   XEndSi = XStartSi + XuvRayMaxSi
                else
                   XEndSi = CoordMax_D(1)*No2Si_V(UnitX_)
                end if

                if(XEndSi <= XStartSi + Small .or. XuvRayStepSi <= Small) then
                   PlotVar_G(i,j,k) = 0.0
                   CYCLE
                end if

                nStep = max(1, ceiling((XEndSi - XStartSi)/XuvRayStepSi))
                DS = (XEndSi - XStartSi)/real(nStep)
                Tau = 0.0

                do iStep = 0, nStep
                   S = real(iStep)*DS
                   XSi = XStartSi + S
                   RSquareSi = XSi*XSi + YSi*YSi + ZSi*ZSi
                   RSi = sqrt(max(Small, RSquareSi))

                   if(RSi <= XuvNeutralR0Si) then
                      NeutralSi = XuvNeutralN0Si
                   else
                      NeutralSi = XuvNeutralN0Si*exp( &
                           -(RSi - XuvNeutralR0Si)/max(Small, XuvNeutralScaleHeightSi))
                   end if
                   NeutralSi = max(NeutralSi, XuvNeutralMinSi)

                   Weight = 1.0
                   if(iStep == 0 .or. iStep == nStep) Weight = 0.5
                   Tau = Tau + Weight*NeutralSi*DS
                end do

                PlotVar_G(i,j,k) = XuvSigmaSi*Tau
                if(abs(PlotVar_G(i,j,k)) < TauMinOut) PlotVar_G(i,j,k) = 0.0
             end do
          end do
       end do

       PlotVarBody = 0.0
       UsePlotVarBody = .false.
       NameTecVar = 'TAUXUV'
       NameTecUnit = '[none]'
       NameIdlUnit = 'none'
       IsFound = .true.
    case('tauxuvo')
       if(allocated(TauOctree_GB)) then
          PlotVar_G = TauOctree_GB(:,:,:,iBlock)
       else
          PlotVar_G = 0.0
       end if
       PlotVarBody = 0.0
       UsePlotVarBody = .false.
       NameTecVar = 'TAUXUVO'
       NameTecUnit = '[none]'
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
    use ModMain, ONLY: nStep
    use ModPhysics, ONLY: No2Si_V, UnitX_
    use ModIO, ONLY: write_prefix, iUnitOut

    integer :: iBlock, i, j, k, iIter
    real :: XSi, YSi, ZSi, DxSi, TauUp, TauNew, TauOld, TauDiffMax
    real :: XFaceRightSi
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
    TauOctree_GB = 0.0

    do iIter = 1, max(1, XuvOctreeMaxIter)
       call message_pass_cell(TauOctree_GB, nWidthIn=1, nProlongOrderIn=1)
       TauDiffMax = 0.0

       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE

          DxSi = CellSize_DB(1,iBlock)*No2Si_V(UnitX_)
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
                   TauNew = TauUp + XuvSigmaSi*neutral_density_si(XSi,YSi,ZSi)*DxSi
                   if(abs(TauNew) < TauMinOut) TauNew = 0.0
                   TauOctree_GB(i,j,k,iBlock) = TauNew
                   TauDiffMax = max(TauDiffMax, abs(TauNew - TauOld))
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

  real function neutral_density_si(XSi, YSi, ZSi)

    real, intent(in) :: XSi, YSi, ZSi
    real :: RSi
    real, parameter :: Small = 1.0e-30

    RSi = sqrt(max(Small, XSi*XSi + YSi*YSi + ZSi*ZSi))
    if(RSi <= XuvNeutralR0Si) then
       neutral_density_si = XuvNeutralN0Si
    else
       neutral_density_si = XuvNeutralN0Si*exp( &
            -(RSi - XuvNeutralR0Si)/max(Small, XuvNeutralScaleHeightSi))
    end if
    neutral_density_si = max(neutral_density_si, XuvNeutralMinSi)

  end function neutral_density_si

end module ModUser
