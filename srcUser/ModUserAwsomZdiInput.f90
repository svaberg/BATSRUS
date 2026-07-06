module ModUserAwsomZdiInput

  implicit none

  private

  public :: read_zdi_user_command

contains

  subroutine read_zdi_user_command(NameCommand, IsHandled)

    use ModUserAwsomZdiConfig, ONLY: read_zdi_magnetogram_param, &
         read_zdi_boundary_param, read_zdi_boundary_check_param
    use ModUserAwsomZdiB0, ONLY: read_zdi_b0_param
    use ModUserAwsomZdiSelfTest, ONLY: read_zdi_selftest_param

    character(len=*), intent(in) :: NameCommand
    logical, intent(out) :: IsHandled
    !--------------------------------------------------------------------------
    IsHandled = .true.

    select case(NameCommand)
    case('#ZDIMAGNETOGRAM')
       call read_zdi_magnetogram_param()
    case('#ZDIB0')
       call read_zdi_b0_param()
    case('#ZDIBOUNDARY')
       call read_zdi_boundary_param()
    case('#ZDIBOUNDARYCHECK')
       call read_zdi_boundary_check_param()
    case('#ZDISELFTEST')
       call read_zdi_selftest_param()
    case default
       IsHandled = .false.
    end select

  end subroutine read_zdi_user_command

end module ModUserAwsomZdiInput
