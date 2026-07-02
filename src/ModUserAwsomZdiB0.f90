module ModUserAwsomZdiB0

  implicit none

  private

  logical, public :: UseZdiB0 = .false.
  real, public :: rSourceSurfaceZdiB0 = 25.0
  integer, public :: nOrderZdiB0 = -1

  logical :: IsZdiB0Initialized = .false.
  integer, public :: nOrderB0 = -1

  real, allocatable :: g_II(:,:), h_II(:,:)
  real, allocatable :: p_II(:,:), Dp_II(:,:)
  real, allocatable :: SinPhi_I(:), CosPhi_I(:)
  real, allocatable :: Sqrt_I(:), SqrtRatio_I(:)
  real, allocatable :: rRsPower_I(:), RmRsPower_I(:), RmRPower_I(:)

  public :: read_zdi_b0_param
  public :: init_zdi_b0
  public :: clean_zdi_b0
  public :: get_zdi_b0

contains

  subroutine read_zdi_b0_param()
    use ModReadParam, ONLY: read_var
    !--------------------------------------------------------------------------
    call read_var('UseZdiB0', UseZdiB0)
    if(.not.UseZdiB0) RETURN

    call read_var('rSourceSurfaceZdiB0', rSourceSurfaceZdiB0)
    call read_var('nOrderZdiB0', nOrderZdiB0)

  end subroutine read_zdi_b0_param

  subroutine init_zdi_b0()

    use ModUtilities, ONLY: CON_stop
    use ModZdiMagnetogram, ONLY: zdi_is_loaded, nZdiOrder, ZdiCoefRe_III, &
         ZdiCoefIm_III

    integer :: l, m, MaxInt
    real :: ConversionFactor, Denom, SourceRatio

    character(len=*), parameter :: NameSub = 'init_zdi_b0'
    !--------------------------------------------------------------------------
    if(.not.UseZdiB0) RETURN

    if(.not.zdi_is_loaded()) call CON_stop( &
         NameSub//': UseZdiB0 requires loaded ZDI coefficients')

    if(rSourceSurfaceZdiB0 <= 1.0) call CON_stop( &
         NameSub//': rSourceSurfaceZdiB0 must be larger than 1')

    call clean_zdi_b0()

    nOrderB0 = nZdiOrder
    if(nOrderZdiB0 > 0) nOrderB0 = min(nOrderB0, nOrderZdiB0)
    if(nOrderB0 < 1) call CON_stop(NameSub//': no usable ZDI harmonics')

    allocate(g_II(0:nOrderB0,0:nOrderB0), h_II(0:nOrderB0,0:nOrderB0))
    allocate(p_II(0:nOrderB0,0:nOrderB0), Dp_II(0:nOrderB0,0:nOrderB0))
    allocate(SinPhi_I(0:nOrderB0), CosPhi_I(0:nOrderB0))
    allocate(rRsPower_I(-1:nOrderB0+2), RmRsPower_I(0:nOrderB0+2), &
         RmRPower_I(0:nOrderB0+2))

    MaxInt = max(nOrderB0**2, 5*nOrderB0, 10)
    allocate(Sqrt_I(MaxInt), SqrtRatio_I(nOrderB0+1))

    do m = 1, MaxInt
       Sqrt_I(m) = sqrt(real(m))
    end do

    SqrtRatio_I(1) = 1.0
    do m = 1, nOrderB0
       SqrtRatio_I(m+1) = SqrtRatio_I(m)*Sqrt_I(2*m-1)/Sqrt_I(2*m)
    end do

    g_II = 0.0
    h_II = 0.0
    do l = 1, nOrderB0
       SourceRatio = (1.0/rSourceSurfaceZdiB0)**(2*l + 1)
       Denom = real(l + 1) + real(l)*SourceRatio
       do m = 0, l
          ConversionFactor = zdi_to_pfss_conversion(l, m)
          g_II(l,m) =  ConversionFactor*ZdiCoefRe_III(1,l,m)/Denom
          h_II(l,m) = -ConversionFactor*ZdiCoefIm_III(1,l,m)/Denom
       end do
    end do

    IsZdiB0Initialized = .true.

  end subroutine init_zdi_b0

  subroutine clean_zdi_b0()
    !--------------------------------------------------------------------------
    if(allocated(g_II)) deallocate(g_II, h_II, p_II, Dp_II, SinPhi_I, &
         CosPhi_I, Sqrt_I, SqrtRatio_I, rRsPower_I, RmRsPower_I, RmRPower_I)
    IsZdiB0Initialized = .false.
    nOrderB0 = -1

  end subroutine clean_zdi_b0

  subroutine get_zdi_b0(Xyz_D, B0_D)

    use ModCoordTransform, ONLY: rot_xyz_sph, xyz_to_rlonlat
    use ModNumConst, ONLY: cHalfPi, cTwoPi
    use ModUserAwsomZdiConfig, ONLY: ZdiFieldScaleNo, ZdiLonShift
    use ModUtilities, ONLY: CON_stop

    real, intent(in) :: Xyz_D(3)
    real, intent(out):: B0_D(3)

    real :: RLonLat_D(3), r, rEval, Theta, Phi, PhiEval
    real :: Bsph_D(3), XyzSph_DD(3,3)

    character(len=*), parameter :: NameSub = 'get_zdi_b0'
    !--------------------------------------------------------------------------
    B0_D = 0.0
    if(.not.UseZdiB0) RETURN
    if(.not.IsZdiB0Initialized) call CON_stop( &
         NameSub//': ZDI B0 was requested but not initialized')

    call xyz_to_rlonlat(Xyz_D, RLonLat_D)
    r = max(RLonLat_D(1), 1.0e-10)
    Phi = RLonLat_D(2)
    Theta = cHalfPi - RLonLat_D(3)
    PhiEval = modulo(Phi - ZdiLonShift, cTwoPi)

    rEval = min(r, rSourceSurfaceZdiB0)
    call get_zdi_pfss_sph(rEval, Theta, PhiEval, Bsph_D)
    if(r > rSourceSurfaceZdiB0) &
         Bsph_D = (rSourceSurfaceZdiB0/r)**2 * Bsph_D

    XyzSph_DD = rot_xyz_sph(Theta, Phi)
    B0_D = matmul(XyzSph_DD, Bsph_D)*ZdiFieldScaleNo

  end subroutine get_zdi_b0

  subroutine get_zdi_pfss_sph(r, Theta, Phi, Bsph_D)

    real, intent(in) :: r, Theta, Phi
    real, intent(out):: Bsph_D(3)

    integer :: l, m
    real :: Br, Btheta, Bphi
    real :: Coef1, Coef2, Coef3, Coef4
    real :: CosTheta, SinTheta
    !--------------------------------------------------------------------------
    call calc_radial_functions(r)

    CosTheta = cos(Theta)
    SinTheta = max(sin(Theta), 1.0e-10)
    call calc_legendre_polynomial(SinTheta, CosTheta)
    call calc_azimuthal_functions(Phi)

    Br = 0.0
    Btheta = 0.0
    Bphi = 0.0

    do m = 0, nOrderB0
       do l = max(1,m), nOrderB0
          Coef1 = real(l + 1)*RmRPower_I(l+2) &
               + RmRsPower_I(l+2)*real(l)*rRsPower_I(l-1)
          Coef3 = RmRPower_I(l+2) - RmRsPower_I(l+2)*rRsPower_I(l-1)
          Coef2 = g_II(l,m)*CosPhi_I(m) + h_II(l,m)*SinPhi_I(m)
          Coef4 = g_II(l,m)*SinPhi_I(m) - h_II(l,m)*CosPhi_I(m)

          Br     = Br     + p_II(l,m)*Coef1*Coef2
          Btheta = Btheta - Dp_II(l,m)*Coef2*Coef3
          Bphi   = Bphi   + p_II(l,m)*real(m)/SinTheta*Coef3*Coef4
       end do
    end do

    Bsph_D = [Br, Btheta, Bphi]

  end subroutine get_zdi_pfss_sph

  subroutine calc_radial_functions(r)

    real, intent(in) :: r

    integer :: l
    real :: RmRs, RmR, rRs
    !--------------------------------------------------------------------------
    RmRs = 1.0/rSourceSurfaceZdiB0
    RmR  = 1.0/r
    rRs  = r/rSourceSurfaceZdiB0

    rRsPower_I(-1) = 1.0/rRs
    rRsPower_I(0)  = 1.0
    RmRsPower_I(0) = 1.0
    RmRPower_I(0)  = 1.0

    do l = 1, nOrderB0 + 2
       RmRsPower_I(l) = RmRsPower_I(l-1)*RmRs
       RmRPower_I(l)  = RmRPower_I(l-1)*RmR
       rRsPower_I(l)  = rRsPower_I(l-1)*rRs
    end do

  end subroutine calc_radial_functions

  subroutine calc_legendre_polynomial(SinTheta, CosTheta)

    real, intent(in) :: SinTheta, CosTheta

    integer :: l, m
    real :: SinThetaM, SinThetaM1
    real :: Coef1, Coef2, Coef3
    !--------------------------------------------------------------------------
    SinThetaM  = 1.0
    SinThetaM1 = 1.0
    p_II  = 0.0
    Dp_II = 0.0

    do m = 0, nOrderB0
       if(m == 0)then
          Coef1 = Sqrt_I(2*m + 1)
       else
          Coef1 = Sqrt_I(2*(2*m + 1))
       end if

       p_II(m,m) = SqrtRatio_I(m+1)*Coef1*SinThetaM
       if(m < nOrderB0) &
            p_II(m+1,m) = p_II(m,m)*Sqrt_I(2*m + 3)*CosTheta

       Dp_II(m,m) = SqrtRatio_I(m+1)*Coef1*m*CosTheta*SinThetaM1
       if(m < nOrderB0) &
            Dp_II(m+1,m) = Sqrt_I(2*m + 3) &
            *(CosTheta*Dp_II(m,m) - SinTheta*p_II(m,m))

       SinThetaM1 = SinThetaM
       SinThetaM  = SinThetaM*SinTheta
    end do

    do m = 0, nOrderB0 - 2
       do l = m + 2, nOrderB0
          Coef1 = Sqrt_I(2*l + 1)/Sqrt_I(l**2 - m**2)
          Coef2 = Sqrt_I(2*l - 1)
          Coef3 = Sqrt_I((l - 1)**2 - m**2)/Sqrt_I(2*l - 3)

          p_II(l,m) = Coef1*(Coef2*CosTheta*p_II(l-1,m) &
               - Coef3*p_II(l-2,m))

          Dp_II(l,m) = Coef1*(Coef2*(CosTheta*Dp_II(l-1,m) &
               - SinTheta*p_II(l-1,m)) - Coef3*Dp_II(l-2,m))
       end do
    end do

    do m = 0, nOrderB0
       do l = m, nOrderB0
          Coef1 = 1.0/Sqrt_I(2*l + 1)
          p_II(l,m)  = p_II(l,m)*Coef1
          Dp_II(l,m) = Dp_II(l,m)*Coef1
       end do
    end do

  end subroutine calc_legendre_polynomial

  subroutine calc_azimuthal_functions(Phi)

    real, intent(in) :: Phi

    integer :: m
    complex :: z, zM
    !--------------------------------------------------------------------------
    CosPhi_I(0) = 1.0
    SinPhi_I(0) = 0.0
    if(nOrderB0 < 1) RETURN

    z = exp(cmplx(0.0, Phi))
    zM = z
    CosPhi_I(1) = real(zM)
    SinPhi_I(1) = aimag(zM)

    do m = 2, nOrderB0
       zM = zM*z
       CosPhi_I(m) = real(zM)
       SinPhi_I(m) = aimag(zM)
    end do

  end subroutine calc_azimuthal_functions

  real function zdi_to_pfss_conversion(l, m)
    use ModNumConst, ONLY: cPi

    integer, intent(in) :: l, m

    real :: ComplexToReal, CondonShortley
    !--------------------------------------------------------------------------
    ComplexToReal = merge(1.0, sqrt(2.0), m == 0)
    CondonShortley = merge(1.0, -1.0, mod(m,2) == 0)

    zdi_to_pfss_conversion = sqrt(real(2*l + 1)) &
         /(CondonShortley*ComplexToReal*sqrt(4.0*cPi))

  end function zdi_to_pfss_conversion

end module ModUserAwsomZdiB0
