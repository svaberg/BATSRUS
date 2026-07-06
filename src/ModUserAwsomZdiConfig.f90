module ModUserAwsomZdiConfig

  implicit none

  private

  ! Direct ZDI coefficient support (experimental, user-module path)
  logical, public :: UseZdiMagnetogram = .false.
  character(len=200), public :: NameZdiCoeffFile = ''
  real, public :: ZdiFieldScaleIo = 1.0
  real, public :: ZdiFieldScaleNo = 1.0
  real, public :: ZdiLonShiftDeg = 0.0
  real, public :: ZdiLonShift = 0.0

  logical, public :: UseZdiBoundary = .false.
  logical, public :: UseZdiBoundaryRadial = .false.
  character(len=20), public :: TypeZdiBoundary = 'off'   ! off/clamp/nudge/brabsb/absb
  character(len=20), public :: TypeZdiRamp = 'cosine'    ! none/linear/cosine
  real, public :: ZdiBcStrength = 1.0
  real, public :: ZdiBcScale = 1.0
  integer, public :: ZdiRampIterStart = 0
  integer, public :: ZdiRampIterStop  = 0
  real, public :: ZdiRampStart = -1.0
  real, public :: ZdiRampStop  = -1.0
  logical, public :: UseZdiBoundaryCheck = .false.
  real, public :: ZdiBoundaryCheckTol = 1.0e-8
  integer :: iZdiRampProgressLogged = -1

  public :: read_zdi_magnetogram_param
  public :: read_zdi_boundary_param
  public :: read_zdi_boundary_check_param
  public :: get_zdi_boundary_ramp_mix
  public :: get_zdi_boundary_mix
  public :: update_zdi_scale_cache

contains

  subroutine read_zdi_magnetogram_param()
    use ModReadParam, ONLY: read_var
    !--------------------------------------------------------------------------
    call read_var('UseZdiMagnetogram', UseZdiMagnetogram)
    if(.not.UseZdiMagnetogram) RETURN

    call read_var('NameZdiCoeffFile', NameZdiCoeffFile)
    call read_var('ZdiFieldScaleIo', ZdiFieldScaleIo)
    call read_var('ZdiLonShiftDeg', ZdiLonShiftDeg)
  end subroutine read_zdi_magnetogram_param

  subroutine read_zdi_boundary_param()
    use ModReadParam, ONLY: read_var
    !--------------------------------------------------------------------------
    call read_var('UseZdiBoundary', UseZdiBoundary)
    iZdiRampProgressLogged = -1
    if(.not.UseZdiBoundary) RETURN

    call read_var('UseZdiBoundaryRadial', UseZdiBoundaryRadial)
    call read_var('TypeZdiBoundary', TypeZdiBoundary)
    call read_var('TypeZdiRamp',     TypeZdiRamp)
    call read_var('ZdiBcStrength',   ZdiBcStrength)
    call read_var('ZdiBcScale',      ZdiBcScale)
    call read_var('ZdiRampIterStart', ZdiRampIterStart)
    call read_var('ZdiRampIterStop',  ZdiRampIterStop)
    call read_var('ZdiRampStart',    ZdiRampStart)
    call read_var('ZdiRampStop',     ZdiRampStop)
  end subroutine read_zdi_boundary_param

  subroutine read_zdi_boundary_check_param()
    use ModReadParam, ONLY: read_var
    !--------------------------------------------------------------------------
    call read_var('UseZdiBoundaryCheck', UseZdiBoundaryCheck)
    if(.not.UseZdiBoundaryCheck) RETURN

    call read_var('ZdiBoundaryCheckTol', ZdiBoundaryCheckTol)
  end subroutine read_zdi_boundary_check_param

  subroutine get_zdi_boundary_ramp_mix(nIteration, tSimulation, FrampZdi, MixZdi)
    use ModNumConst, ONLY: cPi

    integer, intent(in) :: nIteration
    real, intent(in) :: tSimulation
    real, intent(out) :: FrampZdi, MixZdi
    !--------------------------------------------------------------------------
    FrampZdi = 1.0

    if(ZdiRampIterStop > ZdiRampIterStart)then
       FrampZdi = real(nIteration - ZdiRampIterStart) &
            / real(ZdiRampIterStop - ZdiRampIterStart)
    else if(ZdiRampStop > ZdiRampStart)then
       FrampZdi = (tSimulation - ZdiRampStart)/(ZdiRampStop - ZdiRampStart)
    else if(ZdiRampIterStart > 0)then
       FrampZdi = merge(1.0, 0.0, nIteration >= ZdiRampIterStart)
    else if(ZdiRampStart >= 0.0)then
       FrampZdi = merge(1.0, 0.0, tSimulation >= ZdiRampStart)
    end if

    FrampZdi = min(1.0, max(0.0, FrampZdi))
    select case(trim(TypeZdiRamp))
    case('none')
       FrampZdi = merge(1.0, 0.0, FrampZdi > 0.0)
    case('linear')
       ! keep linear
    case default
       ! cosine ramp (default)
       FrampZdi = 0.5*(1.0 - cos(cPi*FrampZdi))
    end select

    select case(trim(TypeZdiBoundary))
    case('clamp', 'brabsb', 'br_absb', 'absb', 'abs_b')
       MixZdi = FrampZdi
    case('nudge')
       MixZdi = min(1.0, max(0.0, ZdiBcStrength*FrampZdi))
    case default
       MixZdi = 0.0
    end select
  end subroutine get_zdi_boundary_ramp_mix

  subroutine get_zdi_boundary_mix(nIteration, tSimulation, DoUseZdiBoundary, MixZdi)
    use ModZdiMagnetogram, ONLY: zdi_is_loaded

    integer, intent(in) :: nIteration
    real, intent(in) :: tSimulation
    logical, intent(out) :: DoUseZdiBoundary
    real, intent(out) :: MixZdi

    real :: FrampZdi
    !--------------------------------------------------------------------------
    DoUseZdiBoundary = UseZdiBoundary .and. zdi_is_loaded()
    if(.not.DoUseZdiBoundary)then
       MixZdi = 0.0
       RETURN
    end if

    call get_zdi_boundary_ramp_mix(nIteration, tSimulation, FrampZdi, MixZdi)
    call log_zdi_ramp_progress(nIteration, tSimulation, FrampZdi, MixZdi)
  end subroutine get_zdi_boundary_mix

  subroutine log_zdi_ramp_progress(nIteration, tSimulation, FrampZdi, MixZdi)
    use BATL_lib, ONLY: iProc
    use ModIO, ONLY: write_prefix, iUnitOut

    integer, intent(in) :: nIteration
    real, intent(in) :: tSimulation, FrampZdi, MixZdi

    integer :: iProgress
    real :: MixMax, Progress
    logical :: DoLog
    !--------------------------------------------------------------------------
    if(iProc /= 0) RETURN
    if(trim(TypeZdiBoundary) == 'off') RETURN
    if(.not.zdi_ramp_has_started(nIteration, tSimulation)) RETURN

    select case(trim(TypeZdiBoundary))
    case('clamp', 'brabsb', 'br_absb', 'absb', 'abs_b')
       MixMax = 1.0
    case('nudge')
       MixMax = min(1.0, max(0.0, ZdiBcStrength))
    case default
       MixMax = 0.0
    end select
    if(MixMax <= 0.0) RETURN

    Progress = min(1.0, max(0.0, MixZdi/MixMax))
    iProgress = min(10, int(10.0*Progress + 1.0e-5))
    DoLog = iProgress > iZdiRampProgressLogged
    if(.not.DoLog) RETURN

    call write_prefix
    write(iUnitOut,'(a,i3,a,2(a,es12.4),a,i10,a,es12.4,a,es12.4)') &
         'ZDI ramp progress=', 10*iProgress, '%', &
         ' MixZdi=', MixZdi, ' FrampZdi=', FrampZdi, &
         ' nIteration=', nIteration, ' tSimulation=', tSimulation, &
         ' MixMax=', MixMax

    iZdiRampProgressLogged = iProgress
  end subroutine log_zdi_ramp_progress

  logical function zdi_ramp_has_started(nIteration, tSimulation)
    integer, intent(in) :: nIteration
    real, intent(in) :: tSimulation
    !--------------------------------------------------------------------------
    if(ZdiRampIterStop > ZdiRampIterStart)then
       zdi_ramp_has_started = nIteration >= ZdiRampIterStart
    else if(ZdiRampStop > ZdiRampStart)then
       zdi_ramp_has_started = tSimulation >= ZdiRampStart
    else if(ZdiRampIterStart > 0)then
       zdi_ramp_has_started = nIteration >= ZdiRampIterStart
    else if(ZdiRampStart >= 0.0)then
       zdi_ramp_has_started = tSimulation >= ZdiRampStart
    else
       zdi_ramp_has_started = .true.
    end if
  end function zdi_ramp_has_started

  subroutine update_zdi_scale_cache(UnitBNoPerIo, DegToRad)
    real, intent(in) :: UnitBNoPerIo, DegToRad
    !--------------------------------------------------------------------------
    ZdiFieldScaleNo = ZdiFieldScaleIo * UnitBNoPerIo
    ZdiLonShift     = ZdiLonShiftDeg * DegToRad
  end subroutine update_zdi_scale_cache

end module ModUserAwsomZdiConfig
