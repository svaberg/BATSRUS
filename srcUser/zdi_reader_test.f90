!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program zdi_reader_test

  use ModZdiMagnetogram
  use ModUtilities, ONLY: CON_stop

  implicit none

  character(len=500) :: NameFile
  integer :: nArg

  !--------------------------------------------------------------------------
  nArg = command_argument_count()
  if(nArg >= 1)then
     call get_command_argument(1, NameFile)
  else
     NameFile = 'zdi_coeff_test.dat'
  end if

  call read_zdi_coeff_file(trim(NameFile))

  write(*,'(a)')    'ZDI coefficient reader summary:'
  write(*,'(a,a)')  '  file   = ', trim(NameFile)
  write(*,'(a,a)')  '  header = ', trim(StringZdiHeader)
  write(*,'(a,3(i0,1x))') '  ints   = ', ZdiHeader_I
  write(*,'(a,i0)') '  lmax   = ', nZdiOrder
  write(*,'(a,i0)') '  nCoeff = ', nZdiCoeffPerSet

  if(index(trim(NameFile), 'zdi_coeff_test.dat') > 0) call check_synthetic_input

  call deallocate_zdi_coeff_arrays

contains
  !============================================================================
  subroutine check_synthetic_input
    !--------------------------------------------------------------------------
    if(trim(StringZdiHeader) /= 'General poloidal plus toroidal field') &
         call CON_stop('zdi_reader_test: wrong header string')

    if(any(ZdiHeader_I /= (/5,3,-3/))) &
         call CON_stop('zdi_reader_test: wrong integer header')

    if(nZdiOrder /= 2) call CON_stop('zdi_reader_test: wrong lmax')
    if(nZdiCoeffPerSet /= 5) call CON_stop('zdi_reader_test: wrong row count')

    call check_coef(1,1,0, 1.0, 0.0)
    call check_coef(1,2,2, 2.2,-2.2)
    call check_coef(2,1,1,11.1,-1.1)
    call check_coef(2,2,0,12.0, 1.2)
    call check_coef(3,2,1,22.1,-2.1)
    call check_coef(3,2,2,22.2,-2.2)

    write(*,'(a)') '  synthetic checks passed'

  end subroutine check_synthetic_input
  !============================================================================
  subroutine check_coef(iSet, l, m, ValueRe, ValueIm)

    integer, intent(in) :: iSet, l, m
    real,    intent(in) :: ValueRe, ValueIm

    real, parameter :: Tol = 1.0e-6
    !--------------------------------------------------------------------------
    if(abs(ZdiCoefRe_III(iSet,l,m) - ValueRe) > Tol) &
         call CON_stop('zdi_reader_test: wrong real coefficient')
    if(abs(ZdiCoefIm_III(iSet,l,m) - ValueIm) > Tol) &
         call CON_stop('zdi_reader_test: wrong imaginary coefficient')

  end subroutine check_coef
  !============================================================================

end program zdi_reader_test
