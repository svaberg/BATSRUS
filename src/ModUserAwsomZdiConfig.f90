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
  character(len=20), public :: TypeZdiBoundary = 'off'   ! off/clamp/nudge
  character(len=20), public :: TypeZdiRamp = 'cosine'    ! none/linear/cosine
  real, public :: ZdiBcStrength = 1.0
  real, public :: ZdiBcScale = 1.0
  integer, public :: ZdiRampIterStart = 0
  integer, public :: ZdiRampIterStop  = 0
  real, public :: ZdiRampStart = -1.0
  real, public :: ZdiRampStop  = -1.0

  public :: read_zdi_magnetogram_param
  public :: read_zdi_boundary_param
  public :: get_zdi_boundary_ramp_mix
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
    case('clamp')
       MixZdi = FrampZdi
    case('nudge')
       MixZdi = min(1.0, max(0.0, ZdiBcStrength*FrampZdi))
    case default
       MixZdi = 0.0
    end select
  end subroutine get_zdi_boundary_ramp_mix

  subroutine update_zdi_scale_cache(UnitBNoPerIo, DegToRad)
    real, intent(in) :: UnitBNoPerIo, DegToRad
    !--------------------------------------------------------------------------
    ZdiFieldScaleNo = ZdiFieldScaleIo * UnitBNoPerIo
    ZdiLonShift     = ZdiLonShiftDeg * DegToRad
  end subroutine update_zdi_scale_cache

end module ModUserAwsomZdiConfig
