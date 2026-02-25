module ModUserAwsomZdiSelfTest

  use BATL_lib, ONLY: test_start, test_stop, iProc
  use ModIoUnit, ONLY: io_unit_new
  use ModUtilities, ONLY: open_file, close_file, CON_stop

  implicit none

  private

  integer, parameter, public :: MaxZdiSelfTestPoint = 512

  logical, public :: UseZdiSelfTest = .false.
  logical :: DoZdiSelfTestSingleCoeff = .true.
  character(len=200) :: NameZdiSelfTestPointFile = 'Param/ZDI/ZDI.selftest.points.dat'
  character(len=200) :: NameZdiSelfTestDumpFile = 'zdi_startup_points.out'
  character(len=200) :: NameZdiSelfTestImageFile = 'zdi_startup_field_2d.out'
  character(len=200) :: NameZdiSelfTestCoeffDumpFile = 'zdi_startup_singlecoeff_l2.out'
  character(len=200) :: NameZdiSelfTestCoeffImagePrefix = 'zdi_startup_singlecoeff'
  integer :: ZdiSelfTestCoeffLmax = 2
  real    :: ZdiSelfTestLonMinDeg = 0.0
  real    :: ZdiSelfTestLonMaxDeg = 360.0
  real    :: ZdiSelfTestDLonDeg   = 10.0
  real    :: ZdiSelfTestLatMinDeg = -90.0
  real    :: ZdiSelfTestLatMaxDeg = 90.0
  real    :: ZdiSelfTestDLatDeg   = 10.0
  integer :: nZdiSelfTestPoint = 0
  real    :: ZdiSelfTestLonDeg_I(MaxZdiSelfTestPoint) = 0.0
  real    :: ZdiSelfTestLatDeg_I(MaxZdiSelfTestPoint) = 0.0

  public :: read_zdi_selftest_param
  public :: run_zdi_selftest_startup_dump

contains

  subroutine read_zdi_selftest_param()
    use ModReadParam, ONLY: read_var
    !--------------------------------------------------------------------------
    call read_var('UseZdiSelfTest', UseZdiSelfTest)
    if(.not.UseZdiSelfTest) RETURN

    call read_var('NameZdiSelfTestPointFile', NameZdiSelfTestPointFile)
    call read_var('NameZdiSelfTestDumpFile', NameZdiSelfTestDumpFile)
    call read_var('NameZdiSelfTestImageFile', NameZdiSelfTestImageFile)
    call read_var('DoZdiSelfTestSingleCoeff', DoZdiSelfTestSingleCoeff)
    call read_var('NameZdiSelfTestCoeffDumpFile', NameZdiSelfTestCoeffDumpFile)
    call read_var('NameZdiSelfTestCoeffImagePrefix', NameZdiSelfTestCoeffImagePrefix)
    call read_var('ZdiSelfTestCoeffLmax', ZdiSelfTestCoeffLmax)
    call read_var('ZdiSelfTestLonMinDeg', ZdiSelfTestLonMinDeg)
    call read_var('ZdiSelfTestLonMaxDeg', ZdiSelfTestLonMaxDeg)
    call read_var('ZdiSelfTestDLonDeg',   ZdiSelfTestDLonDeg)
    call read_var('ZdiSelfTestLatMinDeg', ZdiSelfTestLatMinDeg)
    call read_var('ZdiSelfTestLatMaxDeg', ZdiSelfTestLatMaxDeg)
    call read_var('ZdiSelfTestDLatDeg',   ZdiSelfTestDLatDeg)
  end subroutine read_zdi_selftest_param

  subroutine run_zdi_selftest_startup_dump(NameZdiCoeffFile, ZdiFieldScaleIo, &
       ZdiLonShiftDeg, ZdiLonShift)

    use ModIO, ONLY: write_prefix, iUnitOut
    use ModZdiMagnetogram, ONLY: zdi_is_loaded

    character(len=*), intent(in) :: NameZdiCoeffFile
    real, intent(in) :: ZdiFieldScaleIo, ZdiLonShiftDeg, ZdiLonShift

    logical :: DoTest
    character(len=*), parameter :: NameSub = 'run_zdi_selftest_startup_dump'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    if(.not.zdi_is_loaded()) call CON_stop(NameSub//': ZDI coefficients not loaded')

    call read_zdi_selftest_point_file()

    if(iProc == 0)then
       call write_prefix; write(iUnitOut,*) 'ZDI self-test: nPoint =', nZdiSelfTestPoint
       call write_prefix; write(iUnitOut,*) 'ZDI self-test point file =', &
            trim(NameZdiSelfTestPointFile)
       call write_prefix; write(iUnitOut,*) 'ZDI self-test point dump =', &
            trim(NameZdiSelfTestDumpFile)
       call write_prefix; write(iUnitOut,*) 'ZDI self-test image dump =', &
            trim(NameZdiSelfTestImageFile)
       call write_prefix; write(iUnitOut,*) 'ZDI self-test grid [deg] lon=', &
            ZdiSelfTestLonMinDeg, ZdiSelfTestLonMaxDeg, ZdiSelfTestDLonDeg, &
            ' lat=', ZdiSelfTestLatMinDeg, ZdiSelfTestLatMaxDeg, ZdiSelfTestDLatDeg
    end if

    call dump_zdi_selftest_points(trim(NameZdiSelfTestDumpFile), 'rewind', &
         'Loaded coefficients (all sets active)', &
         NameZdiCoeffFile, ZdiFieldScaleIo, ZdiLonShiftDeg, ZdiLonShift)
    call dump_zdi_selftest_image(trim(NameZdiSelfTestImageFile), &
         'ZDI startup self-test field (all coeffs)', ZdiFieldScaleIo, ZdiLonShift)

    if(DoZdiSelfTestSingleCoeff)then
       call run_zdi_selftest_single_coeff_sweep(NameZdiCoeffFile, ZdiFieldScaleIo, &
            ZdiLonShiftDeg, ZdiLonShift)
    end if

    call test_stop(NameSub, DoTest)
  end subroutine run_zdi_selftest_startup_dump

  subroutine read_zdi_selftest_point_file()
    use ModIO, ONLY: write_prefix, iUnitOut

    integer :: iUnit, iError
    character(len=500) :: StringLine
    logical :: IsExist
    character(len=400) :: NamePointFileResolved

    character(len=*), parameter :: NameSub = 'read_zdi_selftest_point_file'
    !--------------------------------------------------------------------------
    nZdiSelfTestPoint = 0

    if(len_trim(NameZdiSelfTestPointFile) == 0)then
      call set_default_zdi_selftest_points()
      RETURN
    end if

    NamePointFileResolved = trim(NameZdiSelfTestPointFile)
    inquire(file=trim(NamePointFileResolved), exist=IsExist)
    if(.not.IsExist .and. len_trim(NamePointFileResolved) > 0)then
       if(NamePointFileResolved(1:1) /= '/')then
          inquire(file='../'//trim(NamePointFileResolved), exist=IsExist)
          if(IsExist) NamePointFileResolved = '../'//trim(NamePointFileResolved)
       end if
    end if
    if(.not.IsExist)then
       if(iProc == 0)then
          call write_prefix
          write(iUnitOut,*) 'ZDI self-test point file not found, using defaults: ', &
               trim(NameZdiSelfTestPointFile)
       end if
       call set_default_zdi_selftest_points()
       RETURN
    end if

    iUnit = io_unit_new()
    call open_file(iUnit, FILE=trim(NamePointFileResolved), STATUS='old', &
         NameCaller=NameSub)

    do
       read(iUnit, '(a)', iostat=iError) StringLine
       if(iError < 0) EXIT
       if(iError /= 0)then
          call close_file(iUnit, NameCaller=NameSub)
          call CON_stop(NameSub//': error reading point file')
       end if
       StringLine = adjustl(StringLine)
       if(len_trim(StringLine) == 0) CYCLE
       if(StringLine(1:1) == '#' .or. StringLine(1:1) == '!') CYCLE
       if(nZdiSelfTestPoint >= MaxZdiSelfTestPoint)then
          call close_file(iUnit, NameCaller=NameSub)
          call CON_stop(NameSub//': too many points in point file')
       end if
       nZdiSelfTestPoint = nZdiSelfTestPoint + 1
       read(StringLine, *, iostat=iError) ZdiSelfTestLonDeg_I(nZdiSelfTestPoint), &
            ZdiSelfTestLatDeg_I(nZdiSelfTestPoint)
       if(iError /= 0)then
          call close_file(iUnit, NameCaller=NameSub)
          call CON_stop(NameSub//': bad lon/lat line in point file')
       end if
    end do
    call close_file(iUnit, NameCaller=NameSub)

    if(nZdiSelfTestPoint == 0) call set_default_zdi_selftest_points()
  end subroutine read_zdi_selftest_point_file

  subroutine set_default_zdi_selftest_points()
    ! Default points chosen for easy ZDIpy cross-checking.
    !--------------------------------------------------------------------------
    nZdiSelfTestPoint = 10
    ZdiSelfTestLonDeg_I(1:nZdiSelfTestPoint) = &
         [0.0, 10.0, 30.0, 60.0, 90.0, 120.0, 180.0, 240.0, 300.0, 350.0]
    ZdiSelfTestLatDeg_I(1:nZdiSelfTestPoint) = &
         [60.0, 10.0, 30.0, 60.0, 0.0, -30.0, 10.0, -10.0, 45.0, -60.0]
  end subroutine set_default_zdi_selftest_points

  subroutine dump_zdi_selftest_points(NameFile, TypePosition, StringCase, &
       NameZdiCoeffFile, ZdiFieldScaleIo, ZdiLonShiftDeg, ZdiLonShift)

    use ModZdiMagnetogram, ONLY: eval_zdi_surface_field
    use ModNumConst, ONLY: cDegToRad, cTwoPi

    character(len=*), intent(in) :: NameFile, TypePosition, StringCase
    character(len=*), intent(in) :: NameZdiCoeffFile
    real, intent(in) :: ZdiFieldScaleIo, ZdiLonShiftDeg, ZdiLonShift

    integer :: iUnit, iPoint
    real :: LonDeg, LatDeg, LonEval, LatEval
    real :: BrRaw, BphiRaw, BthetaRaw
    real :: BrIo, BphiIo, BthetaIo

    character(len=*), parameter :: NameSub = 'dump_zdi_selftest_points'
    !--------------------------------------------------------------------------
    if(iProc == 0)then
       iUnit = io_unit_new()
       call open_file(iUnit, FILE=trim(NameFile), STATUS='unknown', POSITION=TypePosition, &
            NameCaller=NameSub)

       if(TypePosition /= 'append')then
          write(iUnit,'(a)') '# ZDI startup self-test point dump'
          write(iUnit,'(a,1x,a)') '# coeff_file', trim(NameZdiCoeffFile)
          write(iUnit,'(a,f12.6)') '# zdi_field_scale_io', ZdiFieldScaleIo
          write(iUnit,'(a,f12.6)') '# zdi_lon_shift_deg', ZdiLonShiftDeg
          write(iUnit,'(a)') '# columns: idx lon_deg lat_deg br_raw bphi_raw btheta_raw bmer_raw br_io bphi_io btheta_io bmer_io'
       end if

       write(iUnit,'(a)') ''
       write(iUnit,'(a,1x,a)') '# case', trim(StringCase)
    end if

    do iPoint = 1, nZdiSelfTestPoint
       LonDeg = ZdiSelfTestLonDeg_I(iPoint)
       LatDeg = ZdiSelfTestLatDeg_I(iPoint)
       LonEval = modulo(LonDeg*cDegToRad - ZdiLonShift, cTwoPi)
       LatEval = LatDeg*cDegToRad
       call eval_zdi_surface_field(LonEval, LatEval, BrRaw, BphiRaw, BthetaRaw)
       BrIo    = ZdiFieldScaleIo*BrRaw
       BphiIo  = ZdiFieldScaleIo*BphiRaw
       BthetaIo= ZdiFieldScaleIo*BthetaRaw
       if(iProc == 0)then
          write(iUnit,'(i4,1x,f9.3,1x,f8.3,1x,8(es14.6,1x))') &
               iPoint, LonDeg, LatDeg, BrRaw, BphiRaw, BthetaRaw, -BthetaRaw, &
               BrIo, BphiIo, BthetaIo, -BthetaIo
       end if
    end do

    if(iProc == 0) call close_file(iUnit, NameCaller=NameSub)
  end subroutine dump_zdi_selftest_points

  subroutine dump_zdi_selftest_image(NameFile, StringHeader, ZdiFieldScaleIo, ZdiLonShift)

    use ModZdiMagnetogram, ONLY: eval_zdi_surface_field
    use ModNumConst, ONLY: cDegToRad, cTwoPi
    use ModPlotFile, ONLY: save_plot_file

    character(len=*), intent(in) :: NameFile, StringHeader
    real, intent(in) :: ZdiFieldScaleIo, ZdiLonShift

    integer :: nLon, nLat, iLon, iLat
    real :: LonDeg, LatDeg, LonEval, LatEval
    real :: BrRaw, BphiRaw, BthetaRaw
    real, allocatable :: Lon_I(:), Lat_I(:), Field_VII(:,:,:)

    !--------------------------------------------------------------------------
    call get_zdi_selftest_grid_size(nLon, nLat)

    allocate(Lon_I(nLon), Lat_I(nLat), Field_VII(4,nLon,nLat))

    do iLon = 1, nLon
       Lon_I(iLon) = ZdiSelfTestLonMinDeg + (iLon - 1)*ZdiSelfTestDLonDeg
    end do
    do iLat = 1, nLat
       Lat_I(iLat) = ZdiSelfTestLatMinDeg + (iLat - 1)*ZdiSelfTestDLatDeg
    end do

    Field_VII = 0.0
    do iLat = 1, nLat
       LatDeg = Lat_I(iLat)
       LatEval = LatDeg*cDegToRad
       do iLon = 1, nLon
          LonDeg = Lon_I(iLon)
          LonEval = modulo(LonDeg*cDegToRad - ZdiLonShift, cTwoPi)
          call eval_zdi_surface_field(LonEval, LatEval, BrRaw, BphiRaw, BthetaRaw)
          Field_VII(1,iLon,iLat) = ZdiFieldScaleIo*BrRaw
          Field_VII(2,iLon,iLat) = ZdiFieldScaleIo*BphiRaw
          Field_VII(3,iLon,iLat) = ZdiFieldScaleIo*BthetaRaw
          Field_VII(4,iLon,iLat) = -ZdiFieldScaleIo*BthetaRaw
       end do
    end do

    if(iProc == 0)then
       call save_plot_file(NameFile=trim(NameFile), nDimIn=2, TypeFileIn='ascii', &
            TimeIn=0.0, Coord1In_I=Lon_I, Coord2In_I=Lat_I, &
            VarIn_VII=Field_VII, StringHeaderIn=trim(StringHeader), &
            NameVarIn='Longitude Latitude ZdiBr ZdiBphi ZdiBtheta ZdiBmer')
    end if

    deallocate(Lon_I, Lat_I, Field_VII)
  end subroutine dump_zdi_selftest_image

  subroutine get_zdi_selftest_grid_size(nLon, nLat)

    integer, intent(out) :: nLon, nLat

    real :: nLonFloat, nLatFloat
    character(len=*), parameter :: NameSub = 'get_zdi_selftest_grid_size'
    !--------------------------------------------------------------------------
    if(ZdiSelfTestDLonDeg <= 0.0 .or. ZdiSelfTestDLatDeg <= 0.0)then
       call CON_stop(NameSub//': dLon/dLat must be positive')
    end if
    if(ZdiSelfTestLonMaxDeg < ZdiSelfTestLonMinDeg .or. &
         ZdiSelfTestLatMaxDeg < ZdiSelfTestLatMinDeg)then
       call CON_stop(NameSub//': invalid lon/lat bounds')
    end if

    nLonFloat = (ZdiSelfTestLonMaxDeg - ZdiSelfTestLonMinDeg)/ZdiSelfTestDLonDeg
    nLatFloat = (ZdiSelfTestLatMaxDeg - ZdiSelfTestLatMinDeg)/ZdiSelfTestDLatDeg
    nLon = nint(nLonFloat) + 1
    nLat = nint(nLatFloat) + 1
    if(nLon < 2 .or. nLat < 2) call CON_stop(NameSub//': grid too small')
  end subroutine get_zdi_selftest_grid_size

  subroutine run_zdi_selftest_single_coeff_sweep(NameZdiCoeffFile, ZdiFieldScaleIo, &
       ZdiLonShiftDeg, ZdiLonShift)

    use ModIO, ONLY: write_prefix, iUnitOut
    use ModZdiMagnetogram, ONLY: ZdiCoefRe_III, ZdiCoefIm_III, nZdiSet, nZdiOrder

    character(len=*), intent(in) :: NameZdiCoeffFile
    real, intent(in) :: ZdiFieldScaleIo, ZdiLonShiftDeg, ZdiLonShift

    real, allocatable :: CoefReSave_III(:,:,:), CoefImSave_III(:,:,:)
    integer :: lMaxDump, l, m, iSet, iPart, iCase
    character(len=8)  :: NamePart
    character(len=12) :: NameSet
    character(len=64) :: StringCase, StringSuffix
    character(len=240) :: NameImage

    character(len=*), parameter :: NameSet_I(3) = ['alpha','beta ','gamma']
    !--------------------------------------------------------------------------
    if(.not.allocated(ZdiCoefRe_III)) RETURN

    lMaxDump = min(max(1, ZdiSelfTestCoeffLmax), nZdiOrder)

    allocate(CoefReSave_III(size(ZdiCoefRe_III,1), size(ZdiCoefRe_III,2), size(ZdiCoefRe_III,3)))
    allocate(CoefImSave_III(size(ZdiCoefIm_III,1), size(ZdiCoefIm_III,2), size(ZdiCoefIm_III,3)))
    CoefReSave_III = ZdiCoefRe_III
    CoefImSave_III = ZdiCoefIm_III

    call dump_zdi_selftest_points(trim(NameZdiSelfTestCoeffDumpFile), 'rewind', &
         'Single-coefficient sweep header (subsequent sections append)', &
         NameZdiCoeffFile, ZdiFieldScaleIo, ZdiLonShiftDeg, ZdiLonShift)

    iCase = 0
    do iSet = 1, min(nZdiSet,3)
       NameSet = adjustl(NameSet_I(iSet))
       do l = 1, lMaxDump
          do m = 0, l
             do iPart = 1, 2
                iCase = iCase + 1
                ZdiCoefRe_III = 0.0
                ZdiCoefIm_III = 0.0
                if(iPart == 1)then
                   ZdiCoefRe_III(iSet,l,m) = 1.0
                   NamePart = 're'
                else
                   ZdiCoefIm_III(iSet,l,m) = 1.0
                   NamePart = 'im'
                end if

                write(StringCase,'(a,i3.3,2x,a,2x,a,2x,''l='',i0,2x,''m='',i0)') &
                     'single_coeff', iCase, trim(NameSet), trim(NamePart), l, m
                call dump_zdi_selftest_points(trim(NameZdiSelfTestCoeffDumpFile), 'append', &
                     trim(StringCase), NameZdiCoeffFile, ZdiFieldScaleIo, ZdiLonShiftDeg, ZdiLonShift)

                write(StringSuffix,'(''c'',i3.3,''_'',a,''_l'',i0,''_m'',i0,''_'',a)') &
                     iCase, trim(NameSet), l, m, trim(NamePart)
                NameImage = trim(zdi_selftest_strip_extension(NameZdiSelfTestCoeffImagePrefix))// &
                     '_'//trim(StringSuffix)//'.out'
                call dump_zdi_selftest_image(trim(NameImage), trim(StringCase), &
                     ZdiFieldScaleIo, ZdiLonShift)
             end do
          end do
       end do
    end do

    ZdiCoefRe_III = CoefReSave_III
    ZdiCoefIm_III = CoefImSave_III
    deallocate(CoefReSave_III, CoefImSave_III)

    if(iProc == 0)then
       call write_prefix; write(iUnitOut,*) 'ZDI self-test single-coeff cases dumped =', iCase
       call write_prefix; write(iUnitOut,*) 'ZDI self-test coeff point dump =', &
            trim(NameZdiSelfTestCoeffDumpFile)
       call write_prefix; write(iUnitOut,*) 'ZDI self-test coeff image prefix =', &
            trim(NameZdiSelfTestCoeffImagePrefix)
    end if
  end subroutine run_zdi_selftest_single_coeff_sweep

  function zdi_selftest_strip_extension(NameIn) result(NameOut)

    character(len=*), intent(in) :: NameIn
    character(len=len(NameIn)) :: NameOut

    integer :: i, iDot
    !--------------------------------------------------------------------------
    NameOut = trim(NameIn)
    iDot = 0
    do i = len_trim(NameOut), 1, -1
       if(NameOut(i:i) == '.')then
          iDot = i
          EXIT
       end if
       if(NameOut(i:i) == '/' .or. NameOut(i:i) == '\\') EXIT
    end do
    if(iDot > 0) NameOut(iDot:) = ''
  end function zdi_selftest_strip_extension

end module ModUserAwsomZdiSelfTest
