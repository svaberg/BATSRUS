module ModUserAwsomZdiInput

  implicit none

  private

  public :: read_zdi_user_command

contains

  subroutine read_zdi_user_command(NameCommand, IsHandled)

    use ModUserAwsomZdiConfig, ONLY: read_zdi_magnetogram_param, &
         read_zdi_boundary_param
    use ModUserAwsomZdiSelfTest, ONLY: read_zdi_selftest_param

    character(len=*), intent(in) :: NameCommand
    logical, intent(out) :: IsHandled
    !--------------------------------------------------------------------------
    IsHandled = .true.

    select case(NameCommand)
    case('#ZDIMAGNETOGRAM')
       call read_zdi_magnetogram_param()
    case('#ZDIBOUNDARY')
       call read_zdi_boundary_param()
    case('#ZDISELFTEST')
       call read_zdi_selftest_param()
    case default
       IsHandled = .false.
    end select

  end subroutine read_zdi_user_command

end module ModUserAwsomZdiInput
