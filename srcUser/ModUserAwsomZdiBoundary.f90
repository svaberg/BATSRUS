module ModUserAwsomZdiBoundary

  implicit none

  private

  public :: apply_zdi_boundary_target_cpu
  public :: compose_zdi_inner_b1_bc

contains

  subroutine apply_zdi_boundary_target_cpu(XyzFace_D, B0Face_D, Runit_D, MixZdi, &
       Br1_D, Bt1_D)

    use ModCoordTransform, ONLY: rot_xyz_rlonlat, xyz_to_rlonlat
    use ModNumConst, ONLY: cPi, cTwoPi
    use ModZdiMagnetogram, ONLY: eval_zdi_surface_field
    use ModUserAwsomZdiConfig, ONLY: ZdiBcScale, ZdiFieldScaleNo, ZdiLonShift

    real, intent(in)    :: XyzFace_D(3), B0Face_D(3), Runit_D(3), MixZdi
    real, intent(inout) :: Br1_D(3), Bt1_D(3)

    real :: rZdi, LonZdi, LatZdi, LatZdiEval
    real :: ZdiBr, ZdiBphi, ZdiBtheta
    real :: ZdiRlonLat_D(3), XyzRlonLat_DD(3,3)
    real :: B1tTarget_D(3), B1Target_D(3)
    !--------------------------------------------------------------------------
    call xyz_to_rlonlat(XyzFace_D, rZdi, LonZdi, LatZdi)
    LonZdi = modulo(LonZdi - ZdiLonShift, cTwoPi)

    ! Avoid exact pole singularities in tangential basis evaluation.
    LatZdiEval = min(0.5*cPi - 1.0e-8, max(-0.5*cPi + 1.0e-8, LatZdi))
    call eval_zdi_surface_field(LonZdi, LatZdiEval, ZdiBr, ZdiBphi, ZdiBtheta)

    ! ZDI evaluator returns (Br, Bphi, Btheta) with theta = co-latitude.
    ! rot_xyz_rlonlat uses (Br, BLon, BLat), so BLat = -Btheta.
    XyzRlonLat_DD = rot_xyz_rlonlat(XyzFace_D)
    if(MixZdi > 0.0)then
       ZdiRlonLat_D = [ZdiBcScale*ZdiFieldScaleNo*ZdiBr, &
            ZdiBcScale*ZdiFieldScaleNo*ZdiBphi, &
            -ZdiBcScale*ZdiFieldScaleNo*ZdiBtheta]
       B1Target_D = matmul(ZdiRlonLat_D, transpose(XyzRlonLat_DD))
       B1Target_D = B1Target_D - B0Face_D
       Br1_D = sum(B1Target_D*Runit_D)*Runit_D
       B1tTarget_D = B1Target_D - sum(B1Target_D*Runit_D)*Runit_D
       Bt1_D = (1.0 - MixZdi)*Bt1_D + MixZdi*B1tTarget_D
    else
       ! Free tangential mode: impose only ZDI Br and leave Bt unchanged.
       ZdiRlonLat_D = [ZdiBcScale*ZdiFieldScaleNo*ZdiBr, 0.0, 0.0]
       B1Target_D = matmul(ZdiRlonLat_D, transpose(XyzRlonLat_DD))
       B1Target_D = B1Target_D - B0Face_D
       Br1_D = sum(B1Target_D*Runit_D)*Runit_D
    end if

  end subroutine apply_zdi_boundary_target_cpu

  subroutine compose_zdi_inner_b1_bc(DoUseZdiBoundary, Br1_D, Bt1_D, B1Bc_D)
    logical, intent(in) :: DoUseZdiBoundary
    real, intent(in) :: Br1_D(3), Bt1_D(3)
    real, intent(out) :: B1Bc_D(3)
    !--------------------------------------------------------------------------
    if(DoUseZdiBoundary)then
       B1Bc_D = Br1_D + Bt1_D
    else
       ! Default AWSoM inner BC keeps B1r=0 (tangential B1 only).
       B1Bc_D = Bt1_D
    end if
  end subroutine compose_zdi_inner_b1_bc

end module ModUserAwsomZdiBoundary
