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
    use ModUserAwsomZdiConfig, ONLY: TypeZdiBoundary, ZdiBcScale, &
         ZdiFieldScaleNo, ZdiLonShift

    real, intent(in)    :: XyzFace_D(3), B0Face_D(3), Runit_D(3), MixZdi
    real, intent(inout) :: Br1_D(3), Bt1_D(3)

    real, parameter :: EpsB = 1.0e-30
    real :: rZdi, LonZdi, LatZdi, LatZdiEval
    real :: ZdiBr, ZdiBphi, ZdiBtheta
    real :: ZdiRlonLat_D(3), XyzRlonLat_DD(3,3)
    real :: ZdiTotal_D(3), ZdiBr_D(3), ZdiBt_D(3)
    real :: B1tTarget_D(3), B1Target_D(3), Btotal_D(3)
    real :: BtCurrent_D(3), BtDir_D(3), Bdir_D(3)
    real :: ZdiBrScalar, BtCurrentMag, BtTargetMag, BtNewMag
    real :: ZdiBtMag, ZdiMag, BcurrentMag, BnewMag
    !--------------------------------------------------------------------------
    call xyz_to_rlonlat(XyzFace_D, rZdi, LonZdi, LatZdi)
    LonZdi = modulo(LonZdi - ZdiLonShift, cTwoPi)

    ! Avoid exact pole singularities in tangential basis evaluation.
    LatZdiEval = min(0.5*cPi - 1.0e-8, max(-0.5*cPi + 1.0e-8, LatZdi))
    call eval_zdi_surface_field(LonZdi, LatZdiEval, ZdiBr, ZdiBphi, ZdiBtheta)

    ! ZDI evaluator returns (Br, Bphi, Btheta) with theta = co-latitude.
    ! rot_xyz_rlonlat uses (Br, BLon, BLat), so BLat = -Btheta.
    XyzRlonLat_DD = rot_xyz_rlonlat(XyzFace_D)

    ZdiRlonLat_D = [ZdiBcScale*ZdiFieldScaleNo*ZdiBr, &
         ZdiBcScale*ZdiFieldScaleNo*ZdiBphi, &
         -ZdiBcScale*ZdiFieldScaleNo*ZdiBtheta]
    ZdiTotal_D = matmul(ZdiRlonLat_D, transpose(XyzRlonLat_DD))
    ZdiBrScalar = sum(ZdiTotal_D*Runit_D)
    ZdiBr_D = ZdiBrScalar*Runit_D
    ZdiBt_D = ZdiTotal_D - ZdiBr_D

    Btotal_D = B0Face_D + Br1_D + Bt1_D
    BtCurrent_D = Btotal_D - sum(Btotal_D*Runit_D)*Runit_D

    select case(trim(TypeZdiBoundary))
    case('clamp', 'nudge')
       ! Full-vector ZDI target: Br is fixed and Bt is mixed toward ZDI Bt.
       B1Target_D = ZdiTotal_D - B0Face_D
       Br1_D = sum(B1Target_D*Runit_D)*Runit_D
       B1tTarget_D = B1Target_D - sum(B1Target_D*Runit_D)*Runit_D
       Bt1_D = (1.0 - MixZdi)*Bt1_D + MixZdi*B1tTarget_D

    case('brabsb', 'br_absb')
       ! Fix ZDI Br and ramp |Bt| so that the total |B| approaches ZDI |B|.
       ! The tangential direction is inherited from the current boundary state.
       ZdiMag = sqrt(sum(ZdiTotal_D**2))
       BtTargetMag = sqrt(max(ZdiMag**2 - ZdiBrScalar**2, 0.0))
       BtCurrentMag = sqrt(sum(BtCurrent_D**2))
       if(BtCurrentMag > EpsB)then
          BtDir_D = BtCurrent_D/BtCurrentMag
       else
          ZdiBtMag = sqrt(sum(ZdiBt_D**2))
          if(ZdiBtMag > EpsB)then
             BtDir_D = ZdiBt_D/ZdiBtMag
          else
             BtDir_D = 0.0
          end if
       end if

       BtNewMag = (1.0 - MixZdi)*BtCurrentMag + MixZdi*BtTargetMag
       B1Target_D = ZdiBr_D + BtNewMag*BtDir_D - B0Face_D
       Br1_D = sum(B1Target_D*Runit_D)*Runit_D
       Bt1_D = B1Target_D - Br1_D

    case('absb', 'abs_b')
       ! Ramp only the total field magnitude; preserve the current direction.
       ZdiMag = sqrt(sum(ZdiTotal_D**2))
       BcurrentMag = sqrt(sum(Btotal_D**2))
       if(BcurrentMag > EpsB)then
          Bdir_D = Btotal_D/BcurrentMag
       else if(ZdiMag > EpsB)then
          Bdir_D = ZdiTotal_D/ZdiMag
       else
          Bdir_D = 0.0
       end if

       BnewMag = (1.0 - MixZdi)*BcurrentMag + MixZdi*ZdiMag
       B1Target_D = BnewMag*Bdir_D - B0Face_D
       Br1_D = sum(B1Target_D*Runit_D)*Runit_D
       Bt1_D = B1Target_D - Br1_D

    case default
       ! Free tangential mode: impose only ZDI Br and leave Bt unchanged.
       B1Target_D = ZdiBr_D - B0Face_D
       Br1_D = sum(B1Target_D*Runit_D)*Runit_D
    end select

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
