module ModUserAwsomZdiRuntime

  implicit none

  private

  public :: run_zdi_user_init_runtime

contains

  subroutine run_zdi_user_init_runtime()

    use BATL_lib, ONLY: iProc
    use ModIO, ONLY: write_prefix, iUnitOut
    use ModBatsrusUtility, ONLY: stop_mpi
    use ModNumConst, ONLY: cDegToRad
    use ModPhysics, ONLY: Io2No_V, UnitB_
    use ModZdiMagnetogram, ONLY: read_zdi_coeff_file, zdi_is_loaded, &
         zdi_uses_donati_conjugation, nZdiOrder, nZdiCoeffPerSet, &
         StringZdiHeader, ZdiHeader_I
    use ModUserAwsomZdiConfig, ONLY: UseZdiMagnetogram, NameZdiCoeffFile, &
         ZdiFieldScaleIo, ZdiFieldScaleNo, ZdiLonShiftDeg, ZdiLonShift, &
         UseZdiBoundary, UseZdiBoundaryRadial, TypeZdiBoundary, TypeZdiRamp, &
         ZdiBcStrength, ZdiBcScale, ZdiRampIterStart, ZdiRampIterStop, &
         ZdiRampStart, ZdiRampStop, update_zdi_scale_cache
    use ModUserAwsomZdiB0, ONLY: UseZdiB0, rSourceSurfaceZdiB0, nOrderZdiB0, &
         nOrderB0, init_zdi_b0
    use ModUserAwsomZdiSelfTest, ONLY: UseZdiSelfTest, run_zdi_selftest_startup_dump

    !--------------------------------------------------------------------------
    call update_zdi_scale_cache(Io2No_V(UnitB_), cDegToRad)

    if(iProc == 0)then
       if(UseZdiMagnetogram)then
          call write_prefix; write(iUnitOut,*) &
               'Reading ZDI coefficient file: ', trim(NameZdiCoeffFile)
       end if
    end if

    if(UseZdiMagnetogram)then
       call read_zdi_coeff_file(trim(NameZdiCoeffFile))
       if(iProc == 0)then
          call write_prefix; write(iUnitOut,*) 'ZDI header: ', trim(StringZdiHeader)
          call write_prefix; write(iUnitOut,*) 'ZDI ints: ', ZdiHeader_I
          call write_prefix; write(iUnitOut,*) 'ZDI order / coeff per set: ', &
               nZdiOrder, nZdiCoeffPerSet
          call write_prefix; write(iUnitOut,*) &
               'ZDI conventions: normalized Y/X/Z basis, Btheta=co-latitude, '//&
               'Blat=-Btheta'
          call write_prefix; write(iUnitOut,*) &
               'ZDI Donati -3 conjugation applied = ', zdi_uses_donati_conjugation()
          call write_prefix; write(iUnitOut,*) 'ZDI scales (Io,No)=', &
               ZdiFieldScaleIo, ZdiFieldScaleNo
          call write_prefix; write(iUnitOut,*) 'ZDI lon shift [deg]=', ZdiLonShiftDeg
       end if
       if(UseZdiSelfTest) call run_zdi_selftest_startup_dump( &
            trim(NameZdiCoeffFile), ZdiFieldScaleIo, ZdiLonShiftDeg, ZdiLonShift)
    end if

    if(UseZdiB0)then
       if(.not.zdi_is_loaded()) &
            call stop_mpi('UseZdiB0 requires UseZdiMagnetogram and a valid file')
       call init_zdi_b0()
       if(iProc == 0)then
          call write_prefix; write(iUnitOut,*) &
               'ZDI PFSS/B0 enabled: rSourceSurface=', rSourceSurfaceZdiB0, &
               ' nOrderZdiB0=', nOrderZdiB0, ' nOrderUsed=', nOrderB0
       end if
    end if

    if(UseZdiBoundary .and. .not.zdi_is_loaded()) then
       call stop_mpi('UseZdiBoundary requires UseZdiMagnetogram and a valid file')
    end if

    if(UseZdiBoundary .and. iProc == 0)then
       call write_prefix; write(iUnitOut,*) 'ZDI boundary mode=', trim(TypeZdiBoundary), &
            ', ramp=', trim(TypeZdiRamp), ', strength=', ZdiBcStrength, &
            ', scale=', ZdiBcScale
       call write_prefix; write(iUnitOut,*) &
            'ZDI Br is imposed when UseZdiBoundary=T; legacy radial flag = ', &
            UseZdiBoundaryRadial
       call write_prefix; write(iUnitOut,*) 'ZDI ramp iter/time=', &
            ZdiRampIterStart, ZdiRampIterStop, ZdiRampStart, ZdiRampStop
    end if

  end subroutine run_zdi_user_init_runtime

end module ModUserAwsomZdiRuntime
