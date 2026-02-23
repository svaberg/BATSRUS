!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModZdiMagnetogram

  use ModIoUnit,    ONLY: io_unit_new
  use ModUtilities, ONLY: CON_stop, open_file, close_file
  use ModNumConst,  ONLY: cHalfPi, cPi

  implicit none

  private

  integer, parameter, public :: nZdiSet = 3
  integer, public :: nZdiCoeffPerSet = 0
  integer, public :: nZdiOrder = -1
  character(len=200), public :: StringZdiHeader = ''
  integer, public :: ZdiHeader_I(3) = 0

  ! Real and imaginary parts of the three ZDI coefficient sets:
  ! 1=radial(alpha), 2=poloidal(beta), 3=toroidal(gamma)
  !
  ! Conventions (matched to ZDIpy/core/magneticGeom.py):
  ! - Input is the Donati ZDI coefficient text format with 3 blocks.
  ! - If header flag nPotential == -3 (third integer on line 2), we conjugate
  !   alpha/beta/gamma after reading (same behavior as ZDIpy).
  ! - eval_zdi_surface_field returns (Br, Bphi, Btheta) where theta is
  !   co-latitude (southward). Therefore Blat = -Btheta.
  ! - The spherical-harmonic normalization is
  !     sqrt((2l+1)/(4*pi) * (l-m)!/(l+m)!)
  !   and tangential terms include the standard /(l+1) factor.
  real, allocatable, public :: ZdiCoefRe_III(:,:,:)
  real, allocatable, public :: ZdiCoefIm_III(:,:,:)

  public :: read_zdi_coeff_file
  public :: deallocate_zdi_coeff_arrays
  public :: zdi_is_loaded
  public :: eval_zdi_surface_field
  public :: zdi_uses_donati_conjugation

contains
  !============================================================================
  subroutine read_zdi_coeff_file(NameFile)

    character(len=*), intent(in) :: NameFile

    integer :: iUnit, iError, iSet, iCoeff, l, m
    integer :: lExpect, mExpect
    real    :: CoefRe, CoefIm
    character(len=500) :: StringLine
    logical :: IsEof

    character(len=*), parameter :: NameSub = 'read_zdi_coeff_file'
    !--------------------------------------------------------------------------
    call deallocate_zdi_coeff_arrays

    iUnit = io_unit_new()
    call open_file(iUnit, FILE=NameFile, STATUS='old', NameCaller=NameSub)

    ! First header line is free text
    read(iUnit, '(a)', iostat=iError) StringZdiHeader
    if(iError /= 0)then
       call close_file(iUnit, NameCaller=NameSub)
       call CON_stop(NameSub//': could not read header from '//trim(NameFile))
    end if

    ! Second header line contains three integers; keep them as metadata for now
    call read_next_nonempty_line(iUnit, StringLine, IsEof)
    if(IsEof)then
       call close_file(iUnit, NameCaller=NameSub)
       call CON_stop(NameSub//': missing integer header line in '//trim(NameFile))
    end if
    read(StringLine, *, iostat=iError) ZdiHeader_I(1), ZdiHeader_I(2), ZdiHeader_I(3)
    if(iError /= 0)then
       call close_file(iUnit, NameCaller=NameSub)
       call CON_stop(NameSub//': bad integer header line in '//trim(NameFile))
    end if

    nZdiCoeffPerSet = ZdiHeader_I(1)
    nZdiOrder = zdi_order_from_ncoeff(nZdiCoeffPerSet)
    if(nZdiOrder < 1)then
       call close_file(iUnit, NameCaller=NameSub)
       call CON_stop(NameSub//': invalid coefficient count in '//trim(NameFile))
    end if

    allocate(ZdiCoefRe_III(nZdiSet,0:nZdiOrder,0:nZdiOrder))
    allocate(ZdiCoefIm_III(nZdiSet,0:nZdiOrder,0:nZdiOrder))
    ZdiCoefRe_III = 0.0
    ZdiCoefIm_III = 0.0

    do iSet = 1, nZdiSet
       lExpect = 1
       mExpect = 0
       do iCoeff = 1, nZdiCoeffPerSet
          call read_next_nonempty_line(iUnit, StringLine, IsEof)
          if(IsEof)then
             call close_file(iUnit, NameCaller=NameSub)
             call CON_stop(NameSub//': premature end of file in '//trim(NameFile))
          end if

          read(StringLine, *, iostat=iError) l, m, CoefRe, CoefIm
          if(iError /= 0)then
             call close_file(iUnit, NameCaller=NameSub)
             call CON_stop(NameSub//': bad coefficient row in '//trim(NameFile))
          end if

          if(l /= lExpect .or. m /= mExpect)then
             call close_file(iUnit, NameCaller=NameSub)
             call CON_stop(NameSub//': unexpected (l,m) ordering in '//trim(NameFile))
          end if

          ZdiCoefRe_III(iSet,l,m) = CoefRe
          ZdiCoefIm_III(iSet,l,m) = CoefIm

          mExpect = mExpect + 1
          if(mExpect > lExpect)then
             lExpect = lExpect + 1
             mExpect = 0
          end if
       end do
    end do

    ! Donati "-3" files store the imaginary parts with the opposite sign
    ! relative to the complex-product convention used in the field formulas.
    if(zdi_uses_donati_conjugation()) ZdiCoefIm_III = -ZdiCoefIm_III

    ! Only blank lines are allowed after the coefficients.
    do
       read(iUnit, '(a)', iostat=iError) StringLine
       if(iError < 0) EXIT
       if(iError > 0)then
          call close_file(iUnit, NameCaller=NameSub)
          call CON_stop(NameSub//': read error after coefficients in '//trim(NameFile))
       end if
       if(len_trim(StringLine) > 0)then
          call close_file(iUnit, NameCaller=NameSub)
          call CON_stop(NameSub//': unexpected extra content in '//trim(NameFile))
       end if
    end do

    call close_file(iUnit, NameCaller=NameSub)

  end subroutine read_zdi_coeff_file
  !============================================================================
  subroutine deallocate_zdi_coeff_arrays
    !--------------------------------------------------------------------------
    if(allocated(ZdiCoefRe_III)) deallocate(ZdiCoefRe_III)
    if(allocated(ZdiCoefIm_III)) deallocate(ZdiCoefIm_III)

    nZdiCoeffPerSet = 0
    nZdiOrder = -1
    StringZdiHeader = ''
    ZdiHeader_I = 0

  end subroutine deallocate_zdi_coeff_arrays
  !============================================================================
  logical function zdi_is_loaded()
    !--------------------------------------------------------------------------
    zdi_is_loaded = allocated(ZdiCoefRe_III) .and. allocated(ZdiCoefIm_III) &
         .and. nZdiOrder > 0

  end function zdi_is_loaded
  !============================================================================
  logical function zdi_uses_donati_conjugation()
    !--------------------------------------------------------------------------
    zdi_uses_donati_conjugation = ZdiHeader_I(3) == -3

  end function zdi_uses_donati_conjugation
  !============================================================================
  subroutine eval_zdi_surface_field(Lon, Lat, Br, Bphi, Btheta, &
       BthetaPol, BphiPol, BthetaTor, BphiTor)

    real, intent(in)  :: Lon, Lat
    real, intent(out) :: Br, Bphi, Btheta
    real, optional, intent(out) :: BthetaPol, BphiPol, BthetaTor, BphiTor

    integer :: l, m
    real    :: Theta, SinTheta, CosTheta, InvSinTheta
    real    :: P_II(0:max(0,nZdiOrder),0:max(0,nZdiOrder))
    real    :: dPdTheta_II(0:max(0,nZdiOrder),0:max(0,nZdiOrder))
    real    :: CosMLon, SinMLon
    real    :: AlphaRe, BetaRe, BetaIm, GammaRe, GammaIm
    real    :: aAlpha, aBeta, aGamma, iBeta, iGamma
    real    :: P, dP, PoverSin, NormLm, NormTang, MPoverSin
    real    :: BthetaPolLoc, BphiPolLoc, BthetaTorLoc, BphiTorLoc

    character(len=*), parameter :: NameSub = 'eval_zdi_surface_field'
    !--------------------------------------------------------------------------
    if(.not.zdi_is_loaded()) call CON_stop(NameSub//': ZDI coefficients not loaded')

    ! Input uses latitude; internal formulas use co-latitude theta.
    Theta = cHalfPi - Lat
    SinTheta = sin(Theta)
    CosTheta = cos(Theta)
    InvSinTheta = 0.0
    if(abs(SinTheta) > 1.0e-6) InvSinTheta = 1.0/SinTheta

    call calc_assoc_legendre_theta(Theta, P_II, dPdTheta_II)

    Br = 0.0
    BthetaPolLoc = 0.0
    BphiPolLoc   = 0.0
    BthetaTorLoc = 0.0
    BphiTorLoc   = 0.0

    do l = 1, nZdiOrder
       do m = 0, l
          CosMLon = cos(real(m)*Lon)
          SinMLon = sin(real(m)*Lon)

          AlphaRe = ZdiCoefRe_III(1,l,m); aAlpha = &
               AlphaRe*CosMLon - ZdiCoefIm_III(1,l,m)*SinMLon
          BetaRe  = ZdiCoefRe_III(2,l,m); BetaIm = ZdiCoefIm_III(2,l,m)
          GammaRe = ZdiCoefRe_III(3,l,m); GammaIm = ZdiCoefIm_III(3,l,m)

          aBeta  = BetaRe*CosMLon  - BetaIm*SinMLon
          aGamma = GammaRe*CosMLon - GammaIm*SinMLon
          ! Real(i*c*exp(i m phi)) = -Imag(c*exp(i m phi))
          iBeta  = -(BetaRe*SinMLon  + BetaIm*CosMLon)
          iGamma = -(GammaRe*SinMLon + GammaIm*CosMLon)

          P  = P_II(l,m)
          dP = dPdTheta_II(l,m)
          PoverSin = P*InvSinTheta
          NormLm   = zdi_norm_lm(l,m)
          NormTang = NormLm/real(l+1)
          MPoverSin = real(m)*PoverSin

          ! Matched to ZDIpy/core/magneticGeom.py:
          !   Br    = Re(alpha*Y)
          !   Bclat = -Re(beta*Z + gamma*X)
          !   Blon  = -Re(beta*X - gamma*Z)
          Br = Br + NormLm*aAlpha*P

          BthetaPolLoc = BthetaPolLoc - NormTang*(aBeta*dP)
          BphiPolLoc   = BphiPolLoc   - NormTang*(MPoverSin*iBeta)

          BthetaTorLoc = BthetaTorLoc - NormTang*(MPoverSin*iGamma)
          BphiTorLoc   = BphiTorLoc   + NormTang*(aGamma*dP)
       end do
    end do

    Btheta = BthetaPolLoc + BthetaTorLoc
    Bphi   = BphiPolLoc   + BphiTorLoc

    if(present(BthetaPol)) BthetaPol = BthetaPolLoc
    if(present(BphiPol))   BphiPol   = BphiPolLoc
    if(present(BthetaTor)) BthetaTor = BthetaTorLoc
    if(present(BphiTor))   BphiTor   = BphiTorLoc

  end subroutine eval_zdi_surface_field
  !============================================================================
  subroutine read_next_nonempty_line(iUnit, StringLine, IsEof)

    integer,          intent(in)  :: iUnit
    character(len=*), intent(out) :: StringLine
    logical,          intent(out) :: IsEof

    integer :: iError
    !--------------------------------------------------------------------------
    IsEof = .false.
    StringLine = ''

    do
       read(iUnit, '(a)', iostat=iError) StringLine
       if(iError < 0)then
          IsEof = .true.
          RETURN
       end if
       if(iError > 0) call CON_stop('read_next_nonempty_line: read error')
       if(len_trim(StringLine) == 0) CYCLE
       RETURN
    end do

  end subroutine read_next_nonempty_line
  !============================================================================
  integer function zdi_order_from_ncoeff(nCoeff)

    integer, intent(in) :: nCoeff

    integer :: nCount, lMax
    !--------------------------------------------------------------------------
    zdi_order_from_ncoeff = -1
    if(nCoeff <= 0) RETURN

    nCount = 0
    do lMax = 1, 10000
       nCount = nCount + lMax + 1
       if(nCount == nCoeff) then
          zdi_order_from_ncoeff = lMax
          RETURN
       end if
       if(nCount > nCoeff) then
          RETURN
       end if
    end do

  end function zdi_order_from_ncoeff
  !============================================================================
  real function zdi_norm_lm(l, m)

    integer, intent(in) :: l, m

    integer :: i
    real    :: FactorialRatio
    !--------------------------------------------------------------------------
    if(l < 0 .or. m < 0 .or. m > l) call CON_stop('zdi_norm_lm: bad (l,m)')

    FactorialRatio = 1.0
    do i = l - m + 1, l + m
       FactorialRatio = FactorialRatio/real(i)
    end do

    zdi_norm_lm = sqrt((2.0*real(l) + 1.0)/(4.0*cPi) * FactorialRatio)

  end function zdi_norm_lm
  !============================================================================
  subroutine calc_assoc_legendre_theta(Theta, P_II, dPdTheta_II)

    real, intent(in) :: Theta
    real, intent(out):: P_II(0:,0:)
    real, intent(out):: dPdTheta_II(0:,0:)

    integer :: lMax, l, m
    real    :: SinTheta, CosTheta
    real    :: PrevP
    real    :: Tiny
    !--------------------------------------------------------------------------
    lMax = ubound(P_II,1)
    if(ubound(P_II,2) /= lMax) call CON_stop('calc_assoc_legendre_theta: bad shape')

    P_II = 0.0
    dPdTheta_II = 0.0

    P_II(0,0) = 1.0
    if(lMax < 1) RETURN

    SinTheta = sin(Theta)
    CosTheta = cos(Theta)

    ! Diagonal terms P_m^m (Condon-Shortley phase included)
    do m = 1, lMax
       P_II(m,m) = -(2*m - 1)*SinTheta*P_II(m-1,m-1)
    end do

    ! First off-diagonal P_{m+1}^m
    do m = 0, lMax - 1
       P_II(m+1,m) = (2*m + 1)*CosTheta*P_II(m,m)
    end do

    ! Upward recursion in l
    do m = 0, lMax
       do l = m + 2, lMax
          P_II(l,m) = ((2*l - 1)*CosTheta*P_II(l-1,m) - (l + m - 1)*P_II(l-2,m)) &
               / real(l - m)
       end do
    end do

    Tiny = 1.0e-6
    if(abs(SinTheta) <= Tiny) RETURN

    do l = 1, lMax
       do m = 0, l
          PrevP = 0.0
          if(l-1 >= m) PrevP = P_II(l-1,m)
          dPdTheta_II(l,m) = (real(l)*CosTheta*P_II(l,m) - real(l + m)*PrevP) &
               / SinTheta
       end do
    end do

  end subroutine calc_assoc_legendre_theta
  !============================================================================

end module ModZdiMagnetogram
