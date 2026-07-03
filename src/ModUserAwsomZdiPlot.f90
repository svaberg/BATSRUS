module ModUserAwsomZdiPlot

  implicit none

  private

  public :: set_awsom_zdi_plot_var

contains

  subroutine set_awsom_zdi_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G, NameTecUnit, NameIdlUnit, IsFound)

    use ModAdvance, ONLY: State_VGB
    use ModB0, ONLY: B0_DGB, UseB0
    use ModConst, ONLY: cTwoPi
    use ModCoordTransform, ONLY: rot_xyz_rlonlat, xyz_to_rlonlat
    use ModGeometry, ONLY: r_GB
    use ModPhysics, ONLY: No2Io_V, UnitB_, NameIdlUnit_V
    use ModVarIndexes, ONLY: Bx_, Bz_
    use ModZdiMagnetogram, ONLY: zdi_is_loaded, eval_zdi_surface_field
    use ModUserAwsomZdiConfig, ONLY: UseZdiMagnetogram, ZdiFieldScaleNo, ZdiLonShift
    use BATL_lib, ONLY: Xyz_DGB, MinI, MaxI, MinJ, MaxJ, MinK, MaxK

    integer,          intent(in)    :: iBlock
    character(len=*), intent(in)    :: NameVar
    logical,          intent(in)    :: IsDimensional
    real,             intent(inout) :: PlotVar_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    character(len=*), intent(inout) :: NameTecUnit
    character(len=*), intent(inout) :: NameIdlUnit
    logical,          intent(out)   :: IsFound

    integer :: i, j, k
    real :: UnitB
    real :: FullB_D(3), Brlonlat_D(3), XyzRlonlat_DD(3,3)
    real :: rZdi, LonZdi, LatZdi
    real :: ZdiBr, ZdiBphi, ZdiBtheta
    real :: ZdiBthetaPol, ZdiBphiPol, ZdiBthetaTor, ZdiBphiTor
    !--------------------------------------------------------------------------
    IsFound = .true.

    select case(NameVar)
    case('blon', 'blat', 'bphi', 'btheta')
       if(IsDimensional)then
          UnitB = No2Io_V(UnitB_)
          NameIdlUnit = NameIdlUnit_V(UnitB_)
          NameTecUnit = '['//trim(NameIdlUnit)//']'
       else
          UnitB = 1.0
          NameIdlUnit = '-'
          NameTecUnit = '-'
       end if

       do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
          if(r_GB(i,j,k,iBlock) <= 0.0) CYCLE
          if(UseB0)then
             FullB_D = State_VGB(Bx_:Bz_,i,j,k,iBlock) + B0_DGB(:,i,j,k,iBlock)
          else
             FullB_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)
          end if
          XyzRlonlat_DD = rot_xyz_rlonlat(Xyz_DGB(:,i,j,k,iBlock))
          Brlonlat_D    = matmul(FullB_D, XyzRlonlat_DD)
          select case(NameVar)
          case('blon', 'bphi')
             PlotVar_G(i,j,k) = UnitB*Brlonlat_D(2)
          case('blat')
             PlotVar_G(i,j,k) = UnitB*Brlonlat_D(3)
          case('btheta')
             PlotVar_G(i,j,k) = -UnitB*Brlonlat_D(3)
          end select
       end do; end do; end do

    case('zdibr','zdibphi','zdibtheta','zdiblon','zdiblat', &
         'zdibphip','zdibthetap','zdibphit','zdibthetat')
       if(IsDimensional)then
          UnitB = ZdiFieldScaleNo*No2Io_V(UnitB_)
          NameIdlUnit = NameIdlUnit_V(UnitB_)
          NameTecUnit = '['//trim(NameIdlUnit)//']'
       else
          UnitB = ZdiFieldScaleNo
          NameIdlUnit = '-'
          NameTecUnit = '-'
       end if

       PlotVar_G = 0.0
       if(.not.UseZdiMagnetogram .or. .not.zdi_is_loaded()) RETURN

       do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
          if(r_GB(i,j,k,iBlock) <= 0.0) CYCLE
          call xyz_to_rlonlat(Xyz_DGB(:,i,j,k,iBlock), rZdi, LonZdi, LatZdi)
          LonZdi = modulo(LonZdi - ZdiLonShift, cTwoPi)
          call eval_zdi_surface_field(LonZdi, LatZdi, ZdiBr, ZdiBphi, ZdiBtheta, &
               ZdiBthetaPol, ZdiBphiPol, ZdiBthetaTor, ZdiBphiTor)
          select case(NameVar)
          case('zdibr')
             PlotVar_G(i,j,k) = UnitB*ZdiBr
          case('zdiblon','zdibphi')
             PlotVar_G(i,j,k) = UnitB*ZdiBphi
          case('zdiblat')
             PlotVar_G(i,j,k) = -UnitB*ZdiBtheta
          case('zdibtheta')
             PlotVar_G(i,j,k) = UnitB*ZdiBtheta
          case('zdibphip')
             PlotVar_G(i,j,k) = UnitB*ZdiBphiPol
          case('zdibthetap')
             PlotVar_G(i,j,k) = UnitB*ZdiBthetaPol
          case('zdibphit')
             PlotVar_G(i,j,k) = UnitB*ZdiBphiTor
          case('zdibthetat')
             PlotVar_G(i,j,k) = UnitB*ZdiBthetaTor
          end select
       end do; end do; end do

    case default
       IsFound = .false.
    end select

  end subroutine set_awsom_zdi_plot_var

end module ModUserAwsomZdiPlot
