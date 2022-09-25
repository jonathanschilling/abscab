module abscab
use mod_cel
use mod_compsum

use, intrinsic :: ieee_arithmetic, only: IEEE_VALUE, IEEE_SIGNALING_NAN

implicit none

!> vacuum magnetic permeability in Vs/Am (CODATA-2018)
real(wp), parameter :: MU_0 = 1.25663706212e-6_wp

!> vacuum magnetic permeability, divided by pi
real(wp), parameter :: MU_0_BY_PI = MU_0 / PI

!> vacuum magnetic permeability, divided by 2 pi
real(wp), parameter :: MU_0_BY_2_PI = MU_0 / (2.0_wp * PI)

!> vacuum magnetic permeability, divided by 4 pi
real(wp), parameter :: MU_0_BY_4_PI = MU_0 / (4.0_wp * PI)

contains

!------------ A_z of straight wire segment

!> Compute the normalized axial component of magnetic vector potential of straight wire segment,
!> evaluated along axis of wire segment (rho = 0).
!> This is a special case for points away from the wire ("far-field") for zP < -1 or zP >= 2.
!>
!> @param zP normalized axial coordinate of evaluation location; must not be in [0, 1] (on wire segment)
!> @return normalized axial component of magnetic vector potential
function sws_A_z_ax_f(zP)
    real(wp) :: sws_A_z_ax_f
    real(wp) :: zP
    sws_A_z_ax_f = atanh(1.0_wp / (abs(zP) + abs(1.0_wp - zP)))
end function ! sws_A_z_ax_f

!> Compute the normalized axial component of magnetic vector potential of straight wire segment,
!> evaluated along axis of wire segment (rhoP = 0).
!> This is a special case for points close to the wire ("near-field") for -1 <= zP < 2.
!>
!> @param zP normalized axial coordinate of evaluation location; must not be in [0, 1] (on wire segment)
!> @return normalized axial component of magnetic vector potential
function sws_A_z_ax_n(zP)
    real(wp) :: sws_A_z_ax_n
    real(wp) :: zP
    ! Two negative signs must be able to cancel each other here!
    sws_A_z_ax_n = sign(1.0_wp, zP) * log(zP / (zP - 1.0_wp)) / 2.0_wp
end function ! sws_A_z_ax_n

!> Compute the normalized axial component of magnetic vector potential of straight wire segment,
!> evaluated along axis of wire segment (rho = 0).
!>
!> @param zP normalized axial coordinate of evaluation location; must not be in [0, 1] (on wire segment)
!> @return normalized axial component of magnetic vector potential
function sws_A_z_ax(zP)
    real(wp) :: sws_A_z_ax
    real(wp) :: zP
    if (zP .lt. -1.0_wp .or. zP .ge. 2.0_wp) then
        sws_A_z_ax = sws_A_z_ax_f(zP)
    else
        sws_A_z_ax = sws_A_z_ax_n(zP)
    end if
end function ! sws_A_z_ax

!> Compute the normalized axial component of the magnetic vector potential of a straight wire segment,
!> evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
!> This is a special case for points away from the wire ("far-field") for rhoP > 1.
!>
!> @param rhoP normalized radial coordinate of evaluation location; must not be zero (on wire segment)
!> @return normalized axial component of magnetic vector potential
function sws_A_z_rad_f(rhoP)
    real(wp) :: sws_A_z_rad_f
    real(wp) :: rhoP
    sws_A_z_rad_f = atanh(1.0_wp / (rhoP + hypot(rhoP, 1.0_wp)))
end function ! sws_A_z_rad_f

!> Compute the normalized axial component of the magnetic vector potential of a straight wire segment,
!> evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
!> This is a special case for points close to the wire ("near-field") for rhoP <= 1.
!>
!> @param rhoP normalized radial coordinate of evaluation location; must not be zero (on wire segment)
!> @return normalized axial component of magnetic vector potential
function sws_A_z_rad_n(rhoP)
    real(wp) :: sws_A_z_rad_n
    real(wp) :: rhoP
    real(wp) :: cat, sat, rc, num, den
    cat = 1.0_wp / hypot(rhoP, 1.0_wp) ! cos(atan(...)  )
    sat = sin(atan(rhoP) / 2.0_wp)     ! sin(atan(...)/2)
    rc = rhoP * cat
    num = rc + 1.0_wp + cat
    den = rc + 2.0_wp * sat * sat
    sws_A_z_rad_n = log(num / den) / 2.0_wp
end function ! sws_A_z_rad_n

!> Compute the normalized axial component of the magnetic vector potential of a straight wire segment,
!> evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
!>
!> @param rhoP normalized radial coordinate of evaluation location; must not be zero (on wire segment)
!> @return normalized axial component of magnetic vector potential
function sws_A_z_rad(rhoP)
    real(wp) :: sws_A_z_rad
    real(wp) :: rhoP
    if (rhoP .gt. 1.0_wp) then
        sws_A_z_rad = sws_A_z_rad_f(rhoP)
    else
        sws_A_z_rad = sws_A_z_rad_n(rhoP)
    end if
end function ! sws_A_z_rad

!> Compute the normalized axial component of the magnetic vector potential of a straight wire segment.
!> This formulation is useful for points away from the wire ("far-field")
!> at rhoP >= 1 or zP <= -1 or zP > 2.
!>
!> @param rhoP normalized radial coordinate of evaluation location
!> @param zP normalized axial coordinate of evaluation location
!> @return normalized axial component of magnetic vector potential
function sws_A_z_f(rhoP, zP)
    real(wp) :: sws_A_z_f
    real(wp) :: rhoP
    real(wp) :: zP
    real(wp) :: r_i, r_f
    r_i = hypot(rhoP, zP)
    r_f = hypot(rhoP, 1.0_wp - zP)
    sws_A_z_f = atanh(1.0_wp / (r_i + r_f))
end function ! sws_A_z_f

!> Compute the normalized axial component of the magnetic vector potential of a straight wire segment.
!> This formulation is useful for points close to the wire ("near-field")
!> at rhoP < 1 and -1 < zP <= 2.
!>
!> @param rhoP normalized radial coordinate of evaluation location
!> @param zP normalized axial coordinate of evaluation location
!> @return normalized axial component of magnetic vector potential
function sws_A_z_n(rhoP, zP)
    real(wp) :: sws_A_z_n
    real(wp) :: rhoP
    real(wp) :: zP
    real(wp) :: omz, r_i, r_f, alpha, sinAlphaHalf, &
                beta, sinBetaHalf, Ri_zP, Rf_p_zM1, n

    omz = 1.0_wp - zP

    r_i = hypot(rhoP, zP)
    r_f = hypot(rhoP, omz)

    alpha = atan2(rhoP, zP)
    sinAlphaHalf = sin(alpha / 2.0_wp)

    beta = atan2(rhoP, omz)
    sinBetaHalf = sin(beta / 2.0_wp)

    Ri_zP    = r_i * sinAlphaHalf * sinAlphaHalf ! 0.5 * (r_i - z')
    Rf_p_zM1 = r_f * sinBetaHalf  * sinBetaHalf  ! 0.5 * (r_f - (1 - z'))

    n = Ri_zP + Rf_p_zM1

    sws_A_z_n = (log(1.0_wp + n) - log(n)) / 2.0_wp
end function ! sws_A_z_n

!------------ B_phi of straight wire segment

!> Compute the normalized tangential component of the magnetic field of a straight wire segment,
!> evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
!>
!> @param rhoP normalized radial coordinate of evaluation location
!> @return normalized tangential component of magnetic field
function sws_B_phi_rad(rhoP)
    real(wp) :: sws_B_phi_rad
    real(wp) :: rhoP
    sws_B_phi_rad = 1.0_wp / (rhoP * hypot(rhoP, 1.0_wp))
end function ! sws_B_phi_rad

!> Compute the normalized tangential component of the magnetic field of a straight wire segment.
!> This formulation is useful for points away from the wire ("far-field")
!> at rhoP >= 1 or zP <= 0 or zP >= 1 or rhoP/(1-zP) >= 1 or rhoP/zP >= 1.
!>
!> @param rhoP normalized radial coordinate of evaluation location
!> @param zP normalized axial coordinate of evaluation location
!> @return normalized tangential component of magnetic field
function sws_B_phi_f(rhoP, zP)
    real(wp) :: sws_B_phi_f
    real(wp) :: rhoP
    real(wp) :: zP
    real(wp) :: omz, r_i, r_f, num, den

    omz = 1.0_wp - zP

    r_i = hypot(rhoP, zP)
    r_f = hypot(rhoP, omz)

    num = rhoP * (1.0_wp/r_i + 1.0_wp/r_f)
    den = rhoP * rhoP - zP * omz + r_i * r_f

    sws_B_phi_f = num / den
end function ! sws_B_phi_f

!> Compute the normalized tangential component of the magnetic field of a straight wire segment.
!> This formulation is useful for points close to the wire ("near-field")
!> at rhoP < 1 and 0 < zP < 1 and rhoP/(1-zP) < 1 and rhoP/zP < 1.
!>
!> @param rhoP normalized radial coordinate of evaluation location
!> @param zP normalized axial coordinate of evaluation location
!> @return normalized tangential component of magnetic field
function sws_B_phi_n(rhoP, zP)
    real(wp) :: sws_B_phi_n
    real(wp) :: rhoP
    real(wp) :: zP
    real(wp) :: omz, r_i, r_f, num, den, rfb_omza, &
                alpha, sinAlphaHalf, beta, sinBetaHalf

    omz = 1.0_wp - zP

    r_i = hypot(rhoP, zP)
    r_f = hypot(rhoP, omz)

    num = rhoP * (1.0_wp/r_i + 1.0_wp/r_f)

    alpha = atan2(rhoP, zP)
    sinAlphaHalf = sin(alpha / 2.0_wp)

    beta = atan2(rhoP, omz)
    sinBetaHalf = sin(beta / 2.0_wp)

    ! r_f * sin^2(beta/2) + (1 - z') * sin^2(alpha/2)
    rfb_omza =   r_f * sinBetaHalf  * sinBetaHalf  &
               + omz * sinAlphaHalf * sinAlphaHalf

    !     r_i * r_f - z' * (1 - z')
    ! =   r_i * r_f - r_i * (1 - z') + r_i * (1 - z') - z' * (1 - z')
    ! =   r_i * r_f - r_i * r_f * cos(beta)
    !   + r_i * (1 - z') + (1 - z') * r_i * cos(alpha)
    ! =   r_i *    r_f   * (1 - cos(beta))
    !   + r_i * (1 - z') * (1 - cos(alpha))
    ! = 2 * r_i * [ r_f * sin^2(beta/2) + (1 - z') * sin^2(alpha/2) ]
    den = rhoP * rhoP + 2.0_wp * r_i * rfb_omza

    sws_B_phi_n = num / den
end function ! sws_B_phi_n

!------------ A_phi of circular wire loop

!> Compute the normalized tangential component of the magnetic vector potential of a circular wire loop.
!> This formulation is useful for points away from the wire ("far-field")
!> at rhoP < 1/2 or rhoP > 2 or |zP| >= 1.
!>
!> @param rhoP normalized radial coordinate of evaluation location
!> @param zP normalized axial coordinate of evaluation location
!> @return normalized tangential component of magnetic vector potential
function cwl_A_phi_f(rhoP, zP)
    real(wp) :: cwl_A_phi_f
    real(wp) :: rhoP
    real(wp) :: zP
    real(wp) :: sqrt_kCSqNum, sqrt_kCSqDen, kC, kSq, &
                kCp1, arg1, arg2, C

    sqrt_kCSqNum = hypot(zP, 1.0_wp - rhoP)
    sqrt_kCSqDen = hypot(zP, 1.0_wp + rhoP)

    kC = sqrt_kCSqNum / sqrt_kCSqDen
    kSq = 4.0_wp * rhoP / (sqrt_kCSqDen * sqrt_kCSqDen)

    kCp1 = 1.0_wp + kC
    arg1 = 2.0_wp * sqrt(kC) / kCp1
    arg2 = 2.0_wp / (kCp1 * kCp1 * kCp1)
    C = cel(arg1, 1.0_wp, 0.0_wp, arg2)

    cwl_A_phi_f = kSq/sqrt_kCSqDen * C
end function ! cwl_A_phi_f

!> Compute the normalized tangential component of the magnetic vector potential of a circular wire loop.
!> This formulation is useful for points close to the wire ("near-field")
!> at 1/2 <= rhoP <= 2 and |zP| < 1.
!>
!> @param rhoP normalized radial coordinate of evaluation location
!> @param zP normalized axial coordinate of evaluation location
!> @return normalized tangential component of magnetic vector potential
function cwl_A_phi_n(rhoP, zP)
    real(wp) :: cwl_A_phi_n
    real(wp) :: rhoP
    real(wp) :: zP
    real(wp) :: rhoP_m_1, n, m, num, den, kCSq, prefac, celPart

    rhoP_m_1 = rhoP - 1.0_wp

    n = zP / rhoP_m_1
    m = 1.0_wp + 2.0_wp / rhoP_m_1

    num = n * n + 1.0_wp
    den = n * n + m * m

    kCSq = num / den

    prefac = 1.0_wp / (abs(rhoP - 1.0_wp) * sqrt(den))
    celPart = cel(sqrt(kCSq), 1.0_wp, -1.0_wp, 1.0_wp)

    cwl_A_phi_n = prefac * celPart;
end function ! cwl_A_phi_n

!> Compute the normalized tangential component of the magnetic vector potential of a circular wire loop.
!> This formulation is useful for points along rhoP=1 with |zP| < 1.
!>
!> @param zP normalized axial coordinate of evaluation location
!> @return normalized tangential component of magnetic vector potential
function cwl_A_phi_v(zP)
    real(wp) :: cwl_A_phi_v
    real(wp) :: zP
    real(wp) :: absZp, kCInv

    absZp = abs(zP)

    ! 1/k_c
    kCInv = sqrt(4.0_wp + zP * zP) / absZp

    cwl_A_phi_v = cel(kCInv, 1.0_wp, 1.0_wp, -1.0_wp) / absZp
end function ! cwl_A_phi_v

!------------ B_rho of circular wire loop

!> Compute the normalized radial component of the magnetic field of a circular wire loop.
!> This formulation is useful for points away from the wire ("far-field")
!> at rhoP < 1/2 or rhoP > 2 or |zP| >= 1.
!>
!> @param rhoP normalized radial coordinate of evaluation location
!> @param zP normalized axial coordinate of evaluation location
!> @return normalized radial component of magnetic field
function cwl_B_rho_f(rhoP, zP)
    real(wp) :: cwl_B_rho_f
    real(wp) :: rhoP
    real(wp) :: zP
    real(wp) :: sqrt_kCSqNum, sqrt_kCSqDen, kCSqNum, kCSqDen, &
                kCSq, kC, D, kCp1, arg1, arg2, C, prefac

    sqrt_kCSqNum = hypot(zP, 1.0_wp - rhoP)
    sqrt_kCSqDen = hypot(zP, 1.0_wp + rhoP)

    kCSqNum = sqrt_kCSqNum * sqrt_kCSqNum
    kCSqDen = sqrt_kCSqDen * sqrt_kCSqDen

    kCSq = kCSqNum / kCSqDen
    kC = sqrt(kCSq)

    D = cel(kC, 1.0_wp, 0.0_wp, 1.0_wp)

    kCp1 = 1.0_wp + kC
    arg1 = 2.0_wp * sqrt(kC) / kCp1
    arg2 = 2.0_wp / (kCp1 * kCp1 * kCp1)
    C = cel(arg1, 1.0_wp, 0.0_wp, arg2)

    prefac = 4.0_wp * rhoP / (kCSqDen * sqrt_kCSqDen * kCSqNum)

    cwl_B_rho_f = prefac * zP * (D - C)
end function ! cwl_B_rho_f

!> Compute the normalized radial component of the magnetic field of a circular wire loop.
!> This formulation is useful for points close to the wire ("near-field")
!> at 1/2 <= rhoP <= 2 and |zP| < 1.
!>
!> @param rhoP normalized radial coordinate of evaluation location
!> @param zP normalized axial coordinate of evaluation location
!> @return normalized radial component of magnetic field
function cwl_B_rho_n(rhoP, zP)
    real(wp) :: cwl_B_rho_n
    real(wp) :: rhoP
    real(wp) :: zP
    real(wp) :: rhoP_m_1, rd2, n, m, &
                sqrt_kCSqNum, sqrt_kCSqDen, kCSqNum, kCSqDen, kC, &
                D, kCp1, arg1, arg2, C, zP_rd5, prefac

    rhoP_m_1 = rhoP - 1.0_wp
    rd2 = rhoP_m_1 * rhoP_m_1

    n = zP / rhoP_m_1
    m = 1.0_wp + 2.0_wp / rhoP_m_1

    sqrt_kCSqNum = hypot(n, 1.0_wp)
    sqrt_kCSqDen = hypot(n, m)

    kCSqNum = sqrt_kCSqNum * sqrt_kCSqNum
    kCSqDen = sqrt_kCSqDen * sqrt_kCSqDen

    kC = sqrt_kCSqNum / sqrt_kCSqDen

    D = cel(kC, 1.0_wp, 0.0_wp, 1.0_wp)

    kCp1 = 1.0_wp + kC
    arg1 = 2.0_wp * sqrt(kC) / kCp1
    arg2 = 2.0_wp / (kCp1 * kCp1 * kCp1)
    C = arg2 * cel(arg1, 1.0_wp, 0.0_wp, 1.0_wp)

    ! z' / |rho' - 1|^5
    zP_rd5 = zP / (abs(rhoP_m_1) * rd2 * rd2)

    prefac = 4.0_wp * rhoP / (kCSqDen * sqrt_kCSqDen * kCSqNum)

    cwl_B_rho_n = prefac * zP_rd5 * (D - C)
end function ! cwl_B_rho_n

!> Compute the normalized radial component of the magnetic field of a circular wire loop.
!> This formulation is useful for points along rhoP=1 with |zP| < 1.
!>
!> @param zP normalized axial coordinate of evaluation location
!> @return normalized radial component of magnetic field
function cwl_B_rho_v(zP)
    real(wp) :: cwl_B_rho_v
    real(wp) :: zP
    real(wp) :: zPSq, kCSq, kC, K, E

    zPSq = zP * zP

    kCSq = 1.0_wp / (1.0_wp + 4.0_wp / zPSq)
    kC = sqrt(kCSq)

    K = cel(kC, 1.0_wp, 1.0_wp, 1.0_wp)
    E = cel(kC, 1.0_wp, 1.0_wp, kCSq)

    cwl_B_rho_v = sign(kC / 2.0_wp * ((2.0_wp / zPSq + 1.0_wp) * E - K), zP)
end function ! cwl_B_rho_v

!------------ B_z of circular wire loop

!> Compute the normalized vertical component of the magnetic field of a circular wire loop.
!> This formulation is useful for certain points away from the wire ("far-field")
!> at rhoP < 1/2 or (rhoP <= 2 and |zP| >= 1).
!>
!> @param rhoP normalized radial coordinate of evaluation location
!> @param zP normalized axial coordinate of evaluation location
!> @return normalized vertical component of magnetic field
function cwl_B_z_f1(rhoP, zP)
    real(wp) :: cwl_B_z_f1
    real(wp) :: rhoP
    real(wp) :: zP
    real(wp) :: sqrt_kCSqNum, sqrt_kCSqDen, kC, K, E, D, prefac, comb

    sqrt_kCSqNum = hypot(zP, 1.0_wp - rhoP)
    sqrt_kCSqDen = hypot(zP, 1.0_wp + rhoP)

    kC = sqrt_kCSqNum / sqrt_kCSqDen

    K = cel(kC, 1.0_wp, 1.0_wp, 1.0_wp)
    E = cel(kC, 1.0_wp, 1.0_wp, kC * kC)
    D = cel(kC, 1.0_wp, 0.0_wp, 1.0_wp)

    prefac = 1.0_wp / (sqrt_kCSqDen * sqrt_kCSqNum * sqrt_kCSqNum)
    comb = E - 2.0_wp * K + 2.0_wp * D

    cwl_B_z_f1 = prefac * (E + rhoP * comb)
end function ! cwl_B_z_f1

!> Compute the normalized vertical component of the magnetic field of a circular wire loop.
!> This formulation is useful for certain other points away from the wire ("far-field")
!> at rhoP > 2.
!>
!> @param rhoP normalized radial coordinate of evaluation location
!> @param zP normalized axial coordinate of evaluation location
!> @return normalized vertical component of magnetic field
function cwl_B_z_f2(rhoP, zP)
    real(wp) :: cwl_B_z_f2
    real(wp) :: rhoP
    real(wp) :: zP
    real(wp) :: sqrt_kCSqNum, sqrt_kCSqDen, kC, kCSq, &
                zPSqP1, rhoPSq, t1, t2, a, b, &
                prefac, cdScale, E, D, kCP1, arg1, arg2, C

    sqrt_kCSqNum = hypot(zP, 1.0_wp - rhoP)
    sqrt_kCSqDen = hypot(zP, 1.0_wp + rhoP)

    kC = sqrt_kCSqNum / sqrt_kCSqDen
    kCSq = kC * kC

    zPSqP1 = zP * zP + 1.0_wp
    rhoPSq = rhoP * rhoP
    t1 = zPSqP1 / rhoPSq + 1.0_wp
    t2 = 2.0_wp / rhoP

    ! a is sqrt_kCSqDen normalized to rho'^2
    ! b is sqrt_kCSqNum normalized to rho'^2
    ! a == (z'^2 + (1 + rho')^2) / rho'^2 = (z'^2 + 1)/rho'^2 + 1  +  2/rho'
    ! b == (z'^2 + (1 - rho')^2) / rho'^2 = (z'^2 + 1)/rho'^2 + 1  -  2/rho'
    a = t1 + t2
    b = t1 - t2

    ! 1/prefac = sqrt( z'^2 + (1 + rho')^2)           * (z'^2 + (1 - rho')^2)
    !          = sqrt((z'^2 + (1 + rho')^2) / rho'^2) * (z'^2 + (1 - rho')^2) / rho'^2 * rho'^3
    !          = sqrt(a)                              * b                              * rho'^3
    prefac = 1.0_wp / (sqrt(a) * b * rhoPSq * rhoP)

    cdScale = 1.0_wp + (2.0_wp + zPSqP1 / rhoP) / rhoP

    E = cel(kC, 1.0_wp, 1.0_wp, kCSq)
    D = cel(kC, 1.0_wp, 0.0_wp, 1.0_wp)

    kCP1 = 1.0_wp + kC
    arg1 = 2.0_wp * sqrt(kC) / kCP1
    arg2 = 2.0_wp / (kCP1 * kCP1 * kCP1)
    C = arg2 * cel(arg1, 1.0_wp, 0.0_wp, 1.0_wp)

    ! use C - D for (2 * D - E)/kSq
    cwl_B_z_f2 = prefac * (E + 4.0_wp * (C - D) / cdScale)
end function ! cwl_B_z_f2

!> Compute the normalized vertical component of the magnetic field of a circular wire loop.
!> This formulation is useful for points close to the wire ("near-field")
!> at 1/2 <= rhoP <= 2, but not rhoP=1, and |zP| <= 1.
!>
!> @param rhoP normalized radial coordinate of evaluation location
!> @param zP normalized axial coordinate of evaluation location
!> @return normalized vertical component of magnetic field
function cwl_B_z_n(rhoP, zP)
    real(wp) :: cwl_B_z_n
    real(wp) :: rhoP
    real(wp) :: zP
    real(wp) :: rp1, n, m, prefac, &
                sqrt_kCSqNum, sqrt_kCSqDen, kCSqDen, kC

    rp1 = rhoP - 1.0_wp

    n = zP / rp1
    m = 1.0_wp + 2.0_wp / rp1

    sqrt_kCSqNum = hypot(n, 1.0_wp)
    sqrt_kCSqDen = hypot(n, m)

    kCSqDen = sqrt_kCSqDen * sqrt_kCSqDen

    kC = sqrt_kCSqNum / sqrt_kCSqDen

    prefac = 1.0_wp / (abs(rp1) * rp1 * rp1 * kCSqDen * sqrt_kCSqDen)

    cwl_B_z_n = prefac * cel(kC, kC * kC, 1.0_wp + rhoP, 1.0_wp - rhoP)
end function ! cwl_B_z_n

!> Compute the normalized vertical component of the magnetic field of a circular wire loop.
!> This formulation is useful for points along rhoP=1 with |zP| <= 1.
!>
!> @param zP normalized axial coordinate of evaluation location
!> @return normalized vertical component of magnetic field
function cwl_B_z_v(zP)
    real(wp) :: cwl_B_z_v
    real(wp) :: zP
    real(wp) :: kCSq, kC, f, prefac

    kCSq = zP * zP / (4.0_wp + zP * zP)
    kC = sqrt(kCSq)

    f = zP * zP + 4.0_wp
    prefac = 1.0_wp / (f * sqrt(f))

    cwl_B_z_v = prefac * cel(kC, kCSq, 2.0_wp, 0.0_wp)
end function ! cwl_B_z_v

! --------------------------------------------------

!> Compute the normalized axial component of the magnetic vector potential of a straight wire segment.
!>
!> @param rhoP normalized radial coordinate of evaluation location
!> @param zP normalized axial coordinate of evaluation location
!> @return normalized axial component of magnetic vector potential
function straightWireSegment_A_z(rhoP, zP)
    real(wp) :: straightWireSegment_A_z
    real(wp) :: rhoP
    real(wp) :: zP
    if (rhoP .eq. 0.0_wp) then
        if (zP .lt. 0.0_wp .or. zP .gt. 1.0_wp) then
            straightWireSegment_A_z = sws_A_z_ax(zP)
        else
            write(*,*) "evaluation locations on the wire segment (rho'=", &
              rhoP, " z'=", zP, ") are not allowed"
            straightWireSegment_A_z = IEEE_VALUE(straightWireSegment_A_z, &
                                                 IEEE_SIGNALING_NAN)
        end if
    else if (zP .eq. 0.0_wp .or. zP .eq. 1.0_wp) then
        straightWireSegment_A_z = sws_A_z_rad(rhoP)
    else if (rhoP .ge. 1.0_wp .or. zP .le. -1.0_wp .or. zP .gt. 2.0_wp) then
        straightWireSegment_A_z = sws_A_z_f(rhoP, zP)
    else
        straightWireSegment_A_z = sws_A_z_n(rhoP, zP)
    end if
end function ! straightWireSegment_A_z

!> Compute the normalized tangential component of the magnetic field of a straight wire segment.
!>
!> @param rhoP normalized radial coordinate of evaluation location
!> @param zP normalized axial coordinate of evaluation location
!> @return normalized tangential component of magnetic field
function straightWireSegment_B_phi(rhoP, zP)
    real(wp) :: straightWireSegment_B_phi
    real(wp) :: rhoP
    real(wp) :: zP
    if (rhoP .eq. 0.0_wp) then
        if (zP .lt. 0.0_wp .or. zP .gt. 1.0_wp) then
            straightWireSegment_B_phi = 0.0_wp
        else
            write(*,*) "evaluation locations on the wire segment (rho'=", &
              rhoP, " z'=", zP, ") are not allowed"
            straightWireSegment_B_phi = IEEE_VALUE(straightWireSegment_B_phi, &
                                                   IEEE_SIGNALING_NAN)
        end if
    else if (zP .eq. 0.0_wp .or. zP .eq. 1.0_wp) then
        straightWireSegment_B_phi = sws_B_phi_rad(rhoP)
    else if (rhoP .ge. zP .or. rhoP .ge. 1.0_wp - zP .or. &
             zP .lt. 0.0_wp .or. zP .gt. 1.0_wp) then
        straightWireSegment_B_phi = sws_B_phi_f(rhoP, zP)
    else
        straightWireSegment_B_phi = sws_B_phi_n(rhoP, zP)
    end if
end function ! straightWireSegment_B_phi

!> Geometric part of magnetic vector potential computation for circular wire
!> loop at rho'=1, z'=0 (normalized coordinates). This routine selects special
!> case routines to get the most accurate formulation for given evaluation
!> coordinates.
!>
!> @param rhoP normalized radial evaluation position
!> @param zP   normalized vertical evaluation position
!> @return A_phi: toroidal component of magnetic vector potential: geometric
!>         part (no mu0*I/pi factor included)
function circularWireLoop_A_phi(rhoP, zP)
    real(wp) :: circularWireLoop_A_phi
    real(wp) :: rhoP
    real(wp) :: zP
    if (rhoP .eq. 0.0_wp) then
        circularWireLoop_A_phi = 0.0_wp
    else if (rhoP .lt. 0.5_wp .or. rhoP .gt. 2.0_wp .or. &
             abs(zP) .ge. 1.0_wp) then
        circularWireLoop_A_phi = cwl_A_phi_f(rhoP, zP)
    else if (rhoP .ne. 1.0_wp) then
        circularWireLoop_A_phi = cwl_A_phi_n(rhoP, zP)
    else
        if (zP .ne. 0.0_wp) then
            circularWireLoop_A_phi = cwl_A_phi_v(zP)
        else
            write(*,*) "evaluation at location of wire loop (rho' = 1, z' = 0) is not defined"
            circularWireLoop_A_phi = IEEE_VALUE(circularWireLoop_A_phi, &
                                                IEEE_SIGNALING_NAN)
        end if
    end if
end function ! circularWireLoop_A_phi

!> Geometric part of radial magnetic field computation for circular wire loop at
!> rho'=1, z'=0 (normalized coordinates). This routine selects special case
!> routines to get the most accurate formulation for given evaluation
!> coordinates.
!>
!> @param rhoP normalized radial evaluation position
!> @param zP   normalized vertical evaluation position
!> @return B_rho: radial component of magnetic field: geometric part (no
!>         mu0*I/(pi*a) factor included)
function circularWireLoop_B_rho(rhoP, zP)
    real(wp) :: circularWireLoop_B_rho
    real(wp) :: rhoP
    real(wp) :: zP
    if (rhoP .eq. 0.0_wp .or. zP .eq. 0.0_wp) then
        if (rhoP .ne. 1.0_wp) then
            circularWireLoop_B_rho = 0.0_wp
        else
            write(*,*) "evaluation at location of wire loop (rho' = 1, z' = 0) is not defined"
            circularWireLoop_B_rho = IEEE_VALUE(circularWireLoop_B_rho, &
                                                IEEE_SIGNALING_NAN)
        end if
    else if (rhoP .lt. 0.5_wp .or. rhoP .gt. 2.0_wp .or. &
             abs(zP) .ge. 1.0_wp) then
        circularWireLoop_B_rho = cwl_B_rho_f(rhoP, zP)
    else if (rhoP .ne. 1.0_wp) then
        circularWireLoop_B_rho = cwl_B_rho_n(rhoP, zP)
    else
        circularWireLoop_B_rho = cwl_B_rho_v(zP)
    end if
end function ! circularWireLoop_B_rho

!> Geometric part of vertical magnetic field computation for circular wire loop
!> at rho'=1, z'=0 (normalized coordinates). This routine selects special case
!> routines to get the most accurate formulation for given evaluation
!> coordinates.
!>
!> @param rhoP normalized radial evaluation position
!> @param zP   normalized vertical evaluation position
!> @return B_z: vertical component of magnetic field: geometric part (no
!>         mu0*I/(pi*a) factor included)
function circularWireLoop_B_z(rhoP, zP)
    real(wp) :: circularWireLoop_B_z
    real(wp) :: rhoP
    real(wp) :: zP
    if (rhoP .lt. 0.5_wp .or. &
        (rhoP .le. 2.0_wp .and. abs(zP) .gt. 1.0_wp)) then
        circularWireLoop_B_z = cwl_B_z_f1(rhoP, zP)
    else if (rhoP .gt. 2.0_wp) then
        circularWireLoop_B_z = cwl_B_z_f2(rhoP, zP)
    else if (rhoP .ne. 1.0_wp) then
        circularWireLoop_B_z = cwl_B_z_n(rhoP, zP)
    else
        if (zP .ne. 0.0_wp) then
            circularWireLoop_B_z = cwl_B_z_v(zP)
        else
            write(*,*) "evaluation at location of wire loop (rho' = 1, z' = 0) is not defined"
            circularWireLoop_B_z = IEEE_VALUE(circularWireLoop_B_z, &
                                              IEEE_SIGNALING_NAN)
        end if
    end if
end function ! circularWireLoop_B_z

! --------------------------------------------------

!> Compute the magnetic vector potential of a circular wire loop.
!>
!> @param center  [3: x, y, z] origin of loop (in meters)
!> @param normal  [3: x, y, z] normal vector of loop (in meters); will be
!>                normalized internally
!> @param radius  radius of the wire loop (in meters)
!> @param current loop current (in A)
!> @param evalPos [3: x, y, z][nEvalPos] evaluation locations (in meters)
!> @param vectorPotential [3: A_x, A_y, A_z][nEvalPos] Cartesian components of the magnetic
!>         vector potential evaluated at the given locations (in Tm); has to be allocated on entry
subroutine vectorPotentialCircularFilament(center, normal, radius, current, &
        nEvalPos, evalPos, vectorPotential)

    real(wp), intent(in),  dimension(3)           :: center
    real(wp), intent(in),  dimension(3)           :: normal
    real(wp), intent(in)                          :: radius
    real(wp), intent(in)                          :: current
    integer,  intent(in)                          :: nEvalPos
    real(wp), intent(in),  dimension(3, nEvalPos) :: evalPos
    real(wp), intent(out), dimension(3, nEvalPos) :: vectorPotential

    integer  :: idxEval
    real(wp) :: aPrefactor, nLen2, nLen, eX, eY, eZ, &
                r0x, r0y, r0z, alignedZ, zP, &
                rParallelX, rParallelY, rParallelZ, &
                rPerpX, rPerpY, rPerpZ, &
                alignedRSq, alignedR, eRX, eRY, eRZ, &
                rhoP, aPhi, ePhiX, ePhiY, ePhiZ

    if (.not. ieee_is_finite(radius) .or. radius .le. 0.0_wp) then
        print *, "radius must be finite and positive, but is ", radius
        return
    end if

    aPrefactor = MU_0_BY_PI * current

    ! squared length of normal vector
    nLen2 = normal(1) * normal(1) + normal(2) * normal(2) + normal(3) * normal(3)

    if (nLen2 .eq. 0.0_wp) then
        print *, "length of normal vector must not be zero"
        return
    end if

    ! length of normal vector
    nLen = sqrt(nLen2)

    ! unit normal vector of wire loop
    eX = normal(1) / nLen
    eY = normal(2) / nLen
    eZ = normal(3) / nLen

    do idxEval = 1, nEvalPos

        ! vector from center of wire loop to eval pos
        r0x = evalPos(1, idxEval) - center(1)
        r0y = evalPos(2, idxEval) - center(2)
        r0z = evalPos(3, idxEval) - center(3)

        ! z position along normal of wire loop
        alignedZ = eX * r0x + eY * r0y + eZ * r0z

        ! normalized z component of evaluation location in coordinate system of wire loop
        zP = alignedZ / radius

        ! r0 projected onto axis of wire loop
        rParallelX = alignedZ * eX
        rParallelY = alignedZ * eY
        rParallelZ = alignedZ * eZ

        ! vector perpendicular to axis of wire loop, pointing at evaluation pos
        rPerpX = r0x - rParallelX
        rPerpY = r0y - rParallelY
        rPerpZ = r0z - rParallelZ

        ! perpendicular distance squared between evalPos and axis of wire loop
        alignedRSq = rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ

        ! prevent division-by-zero when computing radial unit vector
        ! A_phi is zero anyway on-axis --> no contribution expected
        if (alignedRSq .gt. 0.0_wp) then

            ! perpendicular distance between evalPos and axis of wire loop
            alignedR = sqrt(alignedRSq)

            ! unit vector in radial direction
            eRX = rPerpX / alignedR
            eRY = rPerpY / alignedR
            eRZ = rPerpZ / alignedR

            ! normalized rho component of evaluation location in coordinate system of wire loop
            rhoP = alignedR / radius

            ! compute tangential component of magnetic vector potential, including current and mu_0
            aPhi = aPrefactor * circularWireLoop_A_phi(rhoP, zP)

            ! compute cross product between e_z and e_rho to get e_phi
            ePhiX = eRY * eZ - eRZ * eY
            ePhiY = eRZ * eX - eRX * eZ
            ePhiZ = eRX * eY - eRY * eX

            ! add contribution from wire loop to result
            vectorPotential(1, idxEval) = aPhi * ePhiX
            vectorPotential(2, idxEval) = aPhi * ePhiY
            vectorPotential(3, idxEval) = aPhi * ePhiZ
        end if ! alignedRSq .gt. 0.0
    end do ! idxEval = 1, nEvalPos
end subroutine ! vectorPotentialCircularFilament

!> Compute the magnetic field of a circular wire loop.
!>
!> @param center  [3: x, y, z] origin of loop (in meters)
!> @param normal  [3: x, y, z] normal vector of loop (in meters); will be
!>                normalized internally
!> @param radius  radius of the wire loop (in meters)
!> @param current loop current (in A)
!> @param evalPos [3: x, y, z][nEvalPos] evaluation locations (in meters)
!> @param magneticField [3: B_x, B_y, B_z][nEvalPos] Cartesian components of the magnetic
!>         field evaluated at the given locations (in T); has to be allocated on entry
subroutine magneticFieldCircularFilament(center, normal, radius, current, &
        nEvalPos, evalPos, magneticField)

    real(wp), intent(in),  dimension(3)           :: center
    real(wp), intent(in),  dimension(3)           :: normal
    real(wp), intent(in)                          :: radius
    real(wp), intent(in)                          :: current
    integer,  intent(in)                          :: nEvalPos
    real(wp), intent(in),  dimension(3, nEvalPos) :: evalPos
    real(wp), intent(out), dimension(3, nEvalPos) :: magneticField

    integer  :: idxEval
    real(wp) :: bPrefactor, nLen2, nLen, eX, eY, eZ, &
                r0x, r0y, r0z, alignedZ, zP, &
                rParallelX, rParallelY, rParallelZ, &
                rPerpX, rPerpY, rPerpZ, &
                alignedRSq, alignedR, eRX, eRY, eRZ, &
                rhoP, bRho, bZ, ePhiX, ePhiY, ePhiZ

    if (.not. ieee_is_finite(radius) .or. radius .le. 0.0_wp) then
        print *, "radius must be finite and positive, but is ", radius
        return
    end if

    bPrefactor = MU_0_BY_PI * current / radius

    ! squared length of normal vector
    nLen2 = normal(1) * normal(1) + normal(2) * normal(2) + normal(3) * normal(3)

    if (nLen2 .eq. 0.0_wp) then
        print *, "length of normal vector must not be zero"
        return
    end if

    ! length of normal vector
    nLen = sqrt(nLen2)

    ! unit normal vector of wire loop
    eX = normal(1) / nLen
    eY = normal(2) / nLen
    eZ = normal(3) / nLen

    do idxEval = 1, nEvalPos

        ! vector from center of wire loop to eval pos
        r0x = evalPos(1, idxEval) - center(1)
        r0y = evalPos(2, idxEval) - center(2)
        r0z = evalPos(3, idxEval) - center(3)

        ! z position along normal of wire loop
        alignedZ = eX * r0x + eY * r0y + eZ * r0z

        ! normalized z component of evaluation location in coordinate system of wire loop
        zP = alignedZ / radius

        ! r0 projected onto axis of wire loop
        rParallelX = alignedZ * eX
        rParallelY = alignedZ * eY
        rParallelZ = alignedZ * eZ

        ! vector perpendicular to axis of wire loop, pointing at evaluation pos
        rPerpX = r0x - rParallelX
        rPerpY = r0y - rParallelY
        rPerpZ = r0z - rParallelZ

        ! perpendicular distance squared between evalPos and axis of wire loop
        alignedRSq = rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ

        if (alignedRSq .gt. 0.0_wp) then
            ! radial unit vector is only defined if evaluation pos is off-axis

            ! perpendicular distance between evalPos and axis of wire loop
            alignedR = sqrt(alignedRSq)

            ! unit vector in radial direction
            eRX = rPerpX / alignedR
            eRY = rPerpY / alignedR
            eRZ = rPerpZ / alignedR

            ! normalized rho component of evaluation location in coordinate system of wire loop
            rhoP = alignedR / radius

            ! compute radial component of normalized magnetic field
            ! and scale by current and mu_0
            bRho = bPrefactor * circularWireLoop_B_rho(rhoP, zP)

            ! add contribution from B_rho of wire loop to result
            magneticField(1, idxEval) = bRho * eRX
            magneticField(2, idxEval) = bRho * eRY
            magneticField(3, idxEval) = bRho * eRZ
        else
            rhoP = 0.0_wp
        end if

        ! compute vertical component of normalized magnetic field
        ! and scale by current and mu_0
        bZ = bPrefactor * circularWireLoop_B_z(rhoP, zP)

        ! add contribution from B_z of wire loop to result
        magneticField(1, idxEval) = magneticField(1, idxEval) + bZ * eX
        magneticField(2, idxEval) = magneticField(2, idxEval) + bZ * eY
        magneticField(3, idxEval) = magneticField(3, idxEval) + bZ * eZ
    end do ! idxEval = 1, nEvalPos
end subroutine ! magneticFieldCircularFilament

!> Compute the magnetic vector potential of a polygon filament
!> at a number of evaluation locations.
!>
!> @param vertices [3: x, y, z][numVertices] points along polygon; in m
!> @param current current along polygon; in A
!> @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
!> @param vectorPotential [3: x, y, z][numEvalPos] target array for magnetic vector potential at evaluation locations; in Tm
!> @param idxSourceStart first index in {@code vertices} to take into account
!> @param idxSourceEnd last index in {@code vertices} to take into account
!> @param idxEvalStart first index in {@code evalPos} to take into account
!> @param idxEvalEnd last index in {@code evalPos} to take into account
!> @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
!>                                of the contributions from the polygon vertices; otherwise, use standard += summation
subroutine kernelVectorPotentialPolygonFilament ( &
        vertices, current, evalPos, vectorPotential, &
        idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd, &
        useCompensatedSummation)

    real(wp), intent(in),  dimension(3, *) :: vertices
    real(wp), intent(in)                   :: current
    real(wp), intent(in),  dimension(3, *) :: evalPos
    real(wp), intent(out), dimension(3, *) :: vectorPotential
    integer,  intent(in)                   :: idxSourceStart
    integer,  intent(in)                   :: idxSourceEnd
    integer,  intent(in)                   :: idxEvalStart
    integer,  intent(in)                   :: idxEvalEnd
    logical,  intent(in)                   :: useCompensatedSummation

    real(wp) :: aPrefactor, x_i, y_i, z_i, x_f, y_f, z_f, &
                dx, dy, dz, l2, l, eX, eY, eZ, r0x, r0y, r0z, &
                alignedZ, zP, rPerpX, rPerpY, rPerpZ, &
                alignedR, rhoP, aParallel

    integer :: istat, idxEval, numEvalPos, idxSource
    real(wp), dimension(:,:), allocatable :: aXSum, aYSum, aZSum

    aPrefactor = MU_0_BY_2_PI * current

    ! setup compensated summation
    if (useCompensatedSummation) then
        numEvalPos = idxEvalEnd - idxEvalStart + 1

        ! need three values (s, cs, ccs) per eval pos --> see mod_compsum
        allocate(aXSum(3, numEvalPos), &
                 aYSum(3, numEvalPos), &
                 aZSum(3, numEvalPos), stat=istat)
        if (istat .ne. 0) then
            print *, "failed to allocate compensated summation buffers: stat=", istat
            return
        end if

        ! initialize target array to zero
        aXSum(:,:) = 0.0_wp
        aYSum(:,:) = 0.0_wp
        aZSum(:,:) = 0.0_wp
    else
        ! initialize target array to zero
        vectorPotential(:, idxEvalStart:idxEvalEnd) = 0.0_wp
    end if ! useCompensatedSummation

    x_i = vertices(1, idxSourceStart)
    y_i = vertices(2, idxSourceStart)
    z_i = vertices(3, idxSourceStart)

    do idxSource = idxSourceStart, idxSourceEnd

        x_f = vertices(1, idxSource + 1)
        y_f = vertices(2, idxSource + 1)
        z_f = vertices(3, idxSource + 1)

        ! vector from start to end of i:th wire segment
        dx = x_f - x_i
        dy = y_f - y_i
        dz = z_f - z_i

        ! squared length of wire segment
        l2 = dx * dx + dy * dy + dz * dz
        if (l2 .eq. 0.0_wp) then
            ! skip zero-length segments: no contribution
            cycle
        end if

        ! length of wire segment
        l = sqrt(l2)

        ! unit vector parallel to wire segment
        eX = dx / l
        eY = dy / l
        eZ = dz / l

        do idxEval = idxEvalStart, idxEvalEnd

            ! vector from start of wire segment to eval pos
            r0x = evalPos(1, idxEval) - x_i
            r0y = evalPos(2, idxEval) - y_i
            r0z = evalPos(3, idxEval) - z_i

            ! z position along axis of wire segment
            alignedZ = eX * r0x + eY * r0y + eZ * r0z

            ! normalized z component of evaluation location in coordinate system of wire segment
            zP = alignedZ / l

            ! vector perpendicular to axis of wire segment, pointing at evaluation pos
            rPerpX = r0x - alignedZ * eX
            rPerpY = r0y - alignedZ * eY
            rPerpZ = r0z - alignedZ * eZ

            ! perpendicular distance between evalPos and axis of wire segment
            alignedR = sqrt(rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ)

            ! normalized rho component of evaluation location in coordinate system of wire segment
            rhoP = alignedR / l

            ! compute parallel component of magnetic vector potential, including current and mu_0
            aParallel = aPrefactor * straightWireSegment_A_z(rhoP, zP)

            ! add contribution from wire segment to result
            if (useCompensatedSummation) then
                call compAdd(aParallel * eX, aXSum(:, idxEval - idxEvalStart + 1))
                call compAdd(aParallel * eY, aYSum(:, idxEval - idxEvalStart + 1))
                call compAdd(aParallel * eZ, aZSum(:, idxEval - idxEvalStart + 1))
            else
                vectorPotential(1, idxEval) = vectorPotential(1, idxEval) + aParallel * eX
                vectorPotential(2, idxEval) = vectorPotential(2, idxEval) + aParallel * eY
                vectorPotential(3, idxEval) = vectorPotential(3, idxEval) + aParallel * eZ
            end if ! useCompensatedSummation
        end do ! idxEval = idxEvalStart, idxEvalEnd

        ! shift to next point
        x_i = x_f
        y_i = y_f
        z_i = z_f
    end do ! idxSource = idxSourceStart, idxSourceEnd

    if (useCompensatedSummation) then
        ! obtain compensated sums from summation objects
        do idxEval = idxEvalStart, idxEvalEnd
            vectorPotential(1, idxEval) = sum(aXSum(:, idxEval - idxEvalStart + 1))
            vectorPotential(2, idxEval) = sum(aYSum(:, idxEval - idxEvalStart + 1))
            vectorPotential(3, idxEval) = sum(aZSum(:, idxEval - idxEvalStart + 1))
        end do

        deallocate(aXSum, aYSum, aZSum, stat=istat)
        if (istat .ne. 0) then
            print *, "failed to deallocate compensated summation buffers: stat=", istat
            return
        end if
    end if
end subroutine ! kernelVectorPotentialPolygonFilament

!> Compute the magnetic vector potential of a polygon filament
!> at a number of evaluation locations.
!>
!> @param vertexSupplier callback to put i-th current carrier polygon vertex into pointData as [3: x, y, z]; in m
!> @param current current along polygon; in A
!> @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
!> @param vectorPotential [3: x, y, z][numEvalPos] target array for magnetic vector potential at evaluation locations; in Tm
!> @param idxSourceStart first index in {@code vertices} to take into account
!> @param idxSourceEnd last index in {@code vertices} to take into account
!> @param idxEvalStart first index in {@code evalPos} to take into account
!> @param idxEvalEnd last index in {@code evalPos} to take into account
!> @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
!>                                of the contributions from the polygon vertices; otherwise, use standard += summation
subroutine kernelVectorPotentialPolygonFilamentVertexSupplier ( &
        vertexSupplier, current, evalPos, vectorPotential, &
        idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd, &
        useCompensatedSummation)

    interface
        subroutine vertexSupplier(i, pointData)
            use mod_kinds, only: wp => dp
            implicit none
            integer,  intent(in)                :: i
            real(wp), intent(out), dimension(3) :: pointData
        end subroutine
    end interface

    real(wp), intent(in)                   :: current
    real(wp), intent(in),  dimension(3, *) :: evalPos
    real(wp), intent(out), dimension(3, *) :: vectorPotential
    integer,  intent(in)                   :: idxSourceStart
    integer,  intent(in)                   :: idxSourceEnd
    integer,  intent(in)                   :: idxEvalStart
    integer,  intent(in)                   :: idxEvalEnd
    logical,  intent(in)                   :: useCompensatedSummation

    real(wp) :: aPrefactor, x_i, y_i, z_i, x_f, y_f, z_f, &
                dx, dy, dz, l2, l, eX, eY, eZ, r0x, r0y, r0z, &
                alignedZ, zP, rPerpX, rPerpY, rPerpZ, &
                alignedR, rhoP, aParallel
    real(wp), dimension(3) :: pointData

    integer :: istat, idxEval, numEvalPos, idxSource
    real(wp), dimension(:,:), allocatable :: aXSum, aYSum, aZSum

    aPrefactor = MU_0_BY_2_PI * current

    ! setup compensated summation
    if (useCompensatedSummation) then
        numEvalPos = idxEvalEnd - idxEvalStart + 1

        ! need three values (s, cs, ccs) per eval pos --> see mod_compsum
        allocate(aXSum(3, numEvalPos), &
                 aYSum(3, numEvalPos), &
                 aZSum(3, numEvalPos), stat=istat)
        if (istat .ne. 0) then
            print *, "failed to allocate compensated summation buffers: stat=", istat
            return
        end if

        ! initialize target array to zero
        aXSum(:,:) = 0.0_wp
        aYSum(:,:) = 0.0_wp
        aZSum(:,:) = 0.0_wp
    else
        ! initialize target array to zero
        vectorPotential(:, idxEvalStart:idxEvalEnd) = 0.0_wp
    end if ! useCompensatedSummation

    call vertexSupplier(idxSourceStart, pointData)
    x_i = pointData(1)
    y_i = pointData(2)
    z_i = pointData(3)

    do idxSource = idxSourceStart, idxSourceEnd

        call vertexSupplier(idxSource + 1, pointData)
        x_f = pointData(1)
        y_f = pointData(2)
        z_f = pointData(3)

        ! vector from start to end of i:th wire segment
        dx = x_f - x_i
        dy = y_f - y_i
        dz = z_f - z_i

        ! squared length of wire segment
        l2 = dx * dx + dy * dy + dz * dz
        if (l2 .eq. 0.0_wp) then
            ! skip zero-length segments: no contribution
            cycle
        end if

        ! length of wire segment
        l = sqrt(l2)

        ! unit vector parallel to wire segment
        eX = dx / l
        eY = dy / l
        eZ = dz / l

        do idxEval = idxEvalStart, idxEvalEnd

            ! vector from start of wire segment to eval pos
            r0x = evalPos(1, idxEval) - x_i
            r0y = evalPos(2, idxEval) - y_i
            r0z = evalPos(3, idxEval) - z_i

            ! z position along axis of wire segment
            alignedZ = eX * r0x + eY * r0y + eZ * r0z

            ! normalized z component of evaluation location in coordinate system of wire segment
            zP = alignedZ / l

            ! vector perpendicular to axis of wire segment, pointing at evaluation pos
            rPerpX = r0x - alignedZ * eX
            rPerpY = r0y - alignedZ * eY
            rPerpZ = r0z - alignedZ * eZ

            ! perpendicular distance between evalPos and axis of wire segment
            alignedR = sqrt(rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ)

            ! normalized rho component of evaluation location in coordinate system of wire segment
            rhoP = alignedR / l

            ! compute parallel component of magnetic vector potential, including current and mu_0
            aParallel = aPrefactor * straightWireSegment_A_z(rhoP, zP)

            ! add contribution from wire segment to result
            if (useCompensatedSummation) then
                call compAdd(aParallel * eX, aXSum(:, idxEval - idxEvalStart + 1))
                call compAdd(aParallel * eY, aYSum(:, idxEval - idxEvalStart + 1))
                call compAdd(aParallel * eZ, aZSum(:, idxEval - idxEvalStart + 1))
            else
                vectorPotential(1, idxEval) = vectorPotential(1, idxEval) + aParallel * eX
                vectorPotential(2, idxEval) = vectorPotential(2, idxEval) + aParallel * eY
                vectorPotential(3, idxEval) = vectorPotential(3, idxEval) + aParallel * eZ
            end if ! useCompensatedSummation
        end do ! idxEval = idxEvalStart, idxEvalEnd

        ! shift to next point
        x_i = x_f
        y_i = y_f
        z_i = z_f
    end do ! idxSource = idxSourceStart, idxSourceEnd

    if (useCompensatedSummation) then
        ! obtain compensated sums from summation objects
        do idxEval = idxEvalStart, idxEvalEnd-1
            vectorPotential(1, idxEval) = sum(aXSum(:, idxEval - idxEvalStart + 1))
            vectorPotential(2, idxEval) = sum(aYSum(:, idxEval - idxEvalStart + 1))
            vectorPotential(3, idxEval) = sum(aZSum(:, idxEval - idxEvalStart + 1))
        end do

        deallocate(aXSum, aYSum, aZSum, stat=istat)
        if (istat .ne. 0) then
            print *, "failed to deallocate compensated summation buffers: stat=", istat
            return
        end if
    end if
end subroutine ! kernelVectorPotentialPolygonFilamentVertexSupplier

!> Compute the magnetic field of a polygon filament
!> at a number of evaluation locations.
!>
!> @param vertices [3: x, y, z][numVertices] points along polygon; in m
!> @param current current along polygon; in A
!> @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
!> @param magneticField [3: x, y, z][numEvalPos] target array for magnetic field at evaluation locations; in T
!> @param idxSourceStart first index in {@code vertices} to take into account
!> @param idxSourceEnd (last+1) index in {@code vertices} to take into account
!> @param idxEvalStart first index in {@code evalPos} to take into account
!> @param idxEvalEnd (last+1) index in {@code evalPos} to take into account
!> @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
!>                                of the contributions from the polygon vertices; otherwise, use standard += summation
subroutine kernelMagneticFieldPolygonFilament ( &
        vertices, current, evalPos, magneticField, &
        idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd, &
        useCompensatedSummation)

    real(wp), intent(in),  dimension(3, *) :: vertices
    real(wp), intent(in)                   :: current
    real(wp), intent(in),  dimension(3, *) :: evalPos
    real(wp), intent(out), dimension(3, *) :: magneticField
    integer,  intent(in)                   :: idxSourceStart
    integer,  intent(in)                   :: idxSourceEnd
    integer,  intent(in)                   :: idxEvalStart
    integer,  intent(in)                   :: idxEvalEnd
    logical,  intent(in)                   :: useCompensatedSummation

    real(wp) :: bPrefactorL, x_i, y_i, z_i, x_f, y_f, z_f, &
                dx, dy, dz, l2, l, eX, eY, eZ, r0x, r0y, r0z, &
                alignedZ, zP, rPerpX, rPerpY, rPerpZ, &
                alignedR, alignedRSq, rhoP, bPrefactor, bPhi, &
                eRX, eRY, eRZ, ePhiX, ePhiY, ePhiZ

    integer :: istat, idxEval, numEvalPos, idxSource
    real(wp), dimension(:,:), allocatable :: bXSum, bYSum, bZSum

    ! needs additional division by length of wire segment!
    bPrefactorL = MU_0_BY_4_PI * current

    ! setup compensated summation
    if (useCompensatedSummation) then
        numEvalPos = idxEvalEnd - idxEvalStart + 1

        ! need three values (s, cs, ccs) per eval pos --> see mod_compsum
        allocate(bXSum(3, numEvalPos), &
                 bYSum(3, numEvalPos), &
                 bZSum(3, numEvalPos), stat=istat)
        if (istat .ne. 0) then
            print *, "failed to allocate compensated summation buffers: stat=", istat
            return
        end if

        ! initialize target array to zero
        bXSum(:,:) = 0.0_wp
        bYSum(:,:) = 0.0_wp
        bZSum(:,:) = 0.0_wp
    else
        ! initialize target array to zero
        magneticField(:, idxEvalStart:idxEvalEnd) = 0.0_wp
    end if ! useCompensatedSummation

    x_i = vertices(1, idxSourceStart)
    y_i = vertices(2, idxSourceStart)
    z_i = vertices(3, idxSourceStart)

    do idxSource = idxSourceStart, idxSourceEnd

        x_f = vertices(1, idxSource + 1)
        y_f = vertices(2, idxSource + 1)
        z_f = vertices(3, idxSource + 1)

        ! vector from start to end of i:th wire segment
        dx = x_f - x_i
        dy = y_f - y_i
        dz = z_f - z_i

        ! squared length of wire segment
        l2 = dx * dx + dy * dy + dz * dz
        if (l2 .eq. 0.0_wp) then
            ! skip zero-length segments: no contribution
            cycle
        end if

        ! length of wire segment
        l = sqrt(l2)

        ! assemble full prefactor for B_phi
        bPrefactor = bPrefactorL / l;

        ! unit vector parallel to wire segment
        eX = dx / l
        eY = dy / l
        eZ = dz / l

        do idxEval = idxEvalStart, idxEvalEnd

            ! vector from start of wire segment to eval pos
            r0x = evalPos(1, idxEval) - x_i
            r0y = evalPos(2, idxEval) - y_i
            r0z = evalPos(3, idxEval) - z_i

            ! z position along axis of wire segment
            alignedZ = eX * r0x + eY * r0y + eZ * r0z

            ! vector perpendicular to axis of wire segment, pointing at evaluation pos
            rPerpX = r0x - alignedZ * eX
            rPerpY = r0y - alignedZ * eY
            rPerpZ = r0z - alignedZ * eZ

            ! perpendicular distance squared between evalPos and axis of wire segment
            alignedRSq = rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ

            if (alignedRSq .gt. 0.0_wp) then

                ! perpendicular distance between evalPos and axis of wire segment
                alignedR = sqrt(alignedRSq)

                ! normalized rho component of evaluation location in coordinate system of wire segment
                rhoP = alignedR / l

                ! normalized z component of evaluation location in coordinate system of wire segment
                zP = alignedZ / l

                ! compute parallel component of magnetic vector potential, including current and mu_0
                bPhi = bPrefactor * straightWireSegment_B_phi(rhoP, zP)

                ! unit vector in radial direction
                eRX = rPerpX / alignedR
                eRY = rPerpY / alignedR
                eRZ = rPerpZ / alignedR

                ! compute cross product between e_z and e_rho to get e_phi
                ePhiX = eY * eRZ - eZ * eRY
                ePhiY = eZ * eRX - eX * eRZ
                ePhiZ = eX * eRY - eY * eRX

                ! add contribution from wire segment to result
                if (useCompensatedSummation) then
                    call compAdd(bPhi * ePhiX, bXSum(:, idxEval - idxEvalStart + 1))
                    call compAdd(bPhi * ePhiY, bYSum(:, idxEval - idxEvalStart + 1))
                    call compAdd(bPhi * ePhiZ, bZSum(:, idxEval - idxEvalStart + 1))
                else
                    magneticField(1, idxEval) = magneticField(1, idxEval) + bPhi * ePhiX
                    magneticField(2, idxEval) = magneticField(2, idxEval) + bPhi * ePhiY
                    magneticField(3, idxEval) = magneticField(3, idxEval) + bPhi * ePhiZ
                end if ! useCompensatedSummation
            end if ! alignedRSq .gt. 0.0_wp
        end do ! idxEval = idxEvalStart, idxEvalEnd

        ! shift to next point
        x_i = x_f
        y_i = y_f
        z_i = z_f
    end do ! idxSource = idxSourceStart, idxSourceEnd

    if (useCompensatedSummation) then
        ! obtain compensated sums from summation objects
        do idxEval = idxEvalStart, idxEvalEnd
            magneticField(1, idxEval) = sum(bXSum(:, idxEval - idxEvalStart + 1))
            magneticField(2, idxEval) = sum(bYSum(:, idxEval - idxEvalStart + 1))
            magneticField(3, idxEval) = sum(bZSum(:, idxEval - idxEvalStart + 1))
        end do

        deallocate(bXSum, bYSum, bZSum, stat=istat)
        if (istat .ne. 0) then
            print *, "failed to deallocate compensated summation buffers: stat=", istat
            return
        end if
    end if
end subroutine ! kernelMagneticFieldPolygonFilament

!> Compute the magnetic field of a polygon filament
!> at a number of evaluation locations.
!>
!> @param vertexSupplier callback to put i-th current carrier polygon vertex into pointData as [3: x, y, z]; in m
!> @param current current along polygon; in A
!> @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
!> @param magneticField [3: x, y, z][numEvalPos] target array for magnetic field at evaluation locations; in T
!> @param idxSourceStart first index in {@code vertices} to take into account
!> @param idxSourceEnd (last+1) index in {@code vertices} to take into account
!> @param idxEvalStart first index in {@code evalPos} to take into account
!> @param idxEvalEnd (last+1) index in {@code evalPos} to take into account
!> @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
!>                                of the contributions from the polygon vertices; otherwise, use standard += summation
subroutine kernelMagneticFieldPolygonFilamentVertexSupplier ( &
        vertexSupplier, current, evalPos, magneticField, &
        idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd, &
        useCompensatedSummation)

    interface
        subroutine vertexSupplier(i, pointData)
            use mod_kinds, only: wp => dp
            implicit none
            integer,  intent(in)                :: i
            real(wp), intent(out), dimension(3) :: pointData
        end subroutine
    end interface

    real(wp), intent(in)                   :: current
    real(wp), intent(in),  dimension(3, *) :: evalPos
    real(wp), intent(out), dimension(3, *) :: magneticField
    integer,  intent(in)                   :: idxSourceStart
    integer,  intent(in)                   :: idxSourceEnd
    integer,  intent(in)                   :: idxEvalStart
    integer,  intent(in)                   :: idxEvalEnd
    logical,  intent(in)                   :: useCompensatedSummation

    real(wp) :: bPrefactorL, x_i, y_i, z_i, x_f, y_f, z_f, &
                dx, dy, dz, l2, l, eX, eY, eZ, r0x, r0y, r0z, &
                alignedZ, zP, rPerpX, rPerpY, rPerpZ, &
                alignedR, alignedRSq, rhoP, bPrefactor, bPhi, &
                eRX, eRY, eRZ, ePhiX, ePhiY, ePhiZ
    real(wp), dimension(3) :: pointData

    integer :: istat, idxEval, numEvalPos, idxSource
    real(wp), dimension(:,:), allocatable :: bXSum, bYSum, bZSum

    ! needs additional division by length of wire segment!
    bPrefactorL = MU_0_BY_4_PI * current

    ! setup compensated summation
    if (useCompensatedSummation) then
        numEvalPos = idxEvalEnd - idxEvalStart + 1

        ! need three values (s, cs, ccs) per eval pos --> see mod_compsum
        allocate(bXSum(3, numEvalPos), &
                 bYSum(3, numEvalPos), &
                 bZSum(3, numEvalPos), stat=istat)
        if (istat .ne. 0) then
            print *, "failed to allocate compensated summation buffers: stat=", istat
            return
        end if

        ! initialize target array to zero
        bXSum(:,:) = 0.0_wp
        bYSum(:,:) = 0.0_wp
        bZSum(:,:) = 0.0_wp
    else
        ! initialize target array to zero
        magneticField(:, idxEvalStart:idxEvalEnd) = 0.0_wp
    end if ! useCompensatedSummation

    call vertexSupplier(idxSourceStart, pointData)
    x_i = pointData(1)
    y_i = pointData(2)
    z_i = pointData(3)

    do idxSource = idxSourceStart, idxSourceEnd

        call vertexSupplier(idxSource + 1, pointData)
        x_f = pointData(1)
        y_f = pointData(2)
        z_f = pointData(3)

        ! vector from start to end of i:th wire segment
        dx = x_f - x_i
        dy = y_f - y_i
        dz = z_f - z_i

        ! squared length of wire segment
        l2 = dx * dx + dy * dy + dz * dz
        if (l2 .eq. 0.0_wp) then
            ! skip zero-length segments: no contribution
            cycle
        end if

        ! length of wire segment
        l = sqrt(l2)

        ! assemble full prefactor for B_phi
        bPrefactor = bPrefactorL / l;

        ! unit vector parallel to wire segment
        eX = dx / l
        eY = dy / l
        eZ = dz / l

        do idxEval = idxEvalStart, idxEvalEnd

            ! vector from start of wire segment to eval pos
            r0x = evalPos(1, idxEval) - x_i
            r0y = evalPos(2, idxEval) - y_i
            r0z = evalPos(3, idxEval) - z_i

            ! z position along axis of wire segment
            alignedZ = eX * r0x + eY * r0y + eZ * r0z

            ! vector perpendicular to axis of wire segment, pointing at evaluation pos
            rPerpX = r0x - alignedZ * eX
            rPerpY = r0y - alignedZ * eY
            rPerpZ = r0z - alignedZ * eZ

            ! perpendicular distance squared between evalPos and axis of wire segment
            alignedRSq = rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ

            if (alignedRSq .gt. 0.0_wp) then

                ! perpendicular distance between evalPos and axis of wire segment
                alignedR = sqrt(alignedRSq)

                ! normalized rho component of evaluation location in coordinate system of wire segment
                rhoP = alignedR / l

                ! normalized z component of evaluation location in coordinate system of wire segment
                zP = alignedZ / l

                ! compute parallel component of magnetic vector potential, including current and mu_0
                bPhi = bPrefactor * straightWireSegment_B_phi(rhoP, zP)

                ! unit vector in radial direction
                eRX = rPerpX / alignedR
                eRY = rPerpY / alignedR
                eRZ = rPerpZ / alignedR

                ! compute cross product between e_z and e_rho to get e_phi
                ePhiX = eY * eRZ - eZ * eRY
                ePhiY = eZ * eRX - eX * eRZ
                ePhiZ = eX * eRY - eY * eRX

                ! add contribution from wire segment to result
                if (useCompensatedSummation) then
                    call compAdd(bPhi * ePhiX, bXSum(:, idxEval - idxEvalStart + 1))
                    call compAdd(bPhi * ePhiY, bYSum(:, idxEval - idxEvalStart + 1))
                    call compAdd(bPhi * ePhiZ, bZSum(:, idxEval - idxEvalStart + 1))
                else
                    magneticField(1, idxEval) = magneticField(1, idxEval) + bPhi * ePhiX
                    magneticField(2, idxEval) = magneticField(2, idxEval) + bPhi * ePhiY
                    magneticField(3, idxEval) = magneticField(3, idxEval) + bPhi * ePhiZ
                end if ! useCompensatedSummation
            end if ! alignedRSq .gt. 0.0_wp
        end do ! idxEval = idxEvalStart, idxEvalEnd

        ! shift to next point
        x_i = x_f
        y_i = y_f
        z_i = z_f
    end do ! idxSource = idxSourceStart, idxSourceEnd

    if (useCompensatedSummation) then
        ! obtain compensated sums from summation objects
        do idxEval = idxEvalStart, idxEvalEnd
            magneticField(1, idxEval) = sum(bXSum(:, idxEval - idxEvalStart + 1))
            magneticField(2, idxEval) = sum(bYSum(:, idxEval - idxEvalStart + 1))
            magneticField(3, idxEval) = sum(bZSum(:, idxEval - idxEvalStart + 1))
        end do

        deallocate(bXSum, bYSum, bZSum, stat=istat)
        if (istat .ne. 0) then
            print *, "failed to deallocate compensated summation buffers: stat=", istat
            return
        end if
    end if
end subroutine ! kernelMagneticFieldPolygonFilamentVertexSupplier

! --------------------------------------------------

!> Compute the magnetic vector potential of a polygon filament
!> at a number of evaluation locations.
!>
!> @param vertices [3: x, y, z][numVertices] points along polygon; in m
!> @param current current along polygon; in A
!> @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
!> @param vectorPotential [3: x, y, z][numEvalPos] target array for magnetic vector potential at evaluation locations; in Tm
!> @param numProcessors number of processors to use for parallelization
!> @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
!>                                of the contributions from the polygon vertices; otherwise, use standard += summation
subroutine vectorPotentialPolygonFilament( &
        numVertices, vertices, current, &
        numEvalPos, evalPos, vectorPotential, &
        numProcessors, useCompensatedSummation)

    integer,  intent(in)                             :: numVertices
    real(wp), intent(in),  dimension(3, numVertices) :: vertices
    real(wp), intent(in)                             :: current
    integer,  intent(in)                             :: numEvalPos
    real(wp), intent(in),  dimension(3, numEvalPos)  :: evalPos
    real(wp), intent(out), dimension(3, numEvalPos)  :: vectorPotential
    integer,  intent(in), optional                   :: numProcessors
    logical,  intent(in), optional            :: useCompensatedSummation

    integer :: actNumProc
    logical :: actUseCompSum

    integer :: i, idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd, &
               nThreads, nSourcePerThread, nEvalPerThread, istat, idxThread
    real(wp), dimension(3) :: sumX, sumY, sumZ
    real(wp), dimension(:, :, :), allocatable :: vecPotContribs

    ! handle optional arguments
    if (present(numProcessors)) then
        actNumProc = numProcessors
    else
        ! default: single-threaded
        actNumProc = 1
    end if ! present(numProcessors)

    if (present(useCompensatedSummation)) then
        actUseCompSum = useCompensatedSummation
    else
        ! default: use compensated summation
        actUseCompSum = .true.
    end if

    if (numVertices .lt. 2) then
        print *, "need at least 2 vertices, but got only ", numVertices
        return
    end if

    if (actNumProc .lt. 1) then
        print *, "need at least 1 processor, but got only ", actNumProc
        return
    end if

    if (current .eq. 0.0_wp) then
        vectorPotential(:, :) = 0.0_wp
        return
    end if

    if (actNumProc .eq. 1) then
        ! single-threaded call
        idxSourceStart = 1
        idxSourceEnd   = numVertices-1
        idxEvalStart   = 1
        idxEvalEnd     = numEvalPos
        call kernelVectorPotentialPolygonFilament( &
            vertices, current, &
            evalPos, &
            vectorPotential, &
            idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd, &
            actUseCompSum)
    else
        ! use multithreading

        if (numVertices-1 .gt. numEvalPos) then
            ! parallelize over nSource-1

            ! Note that each thread needs its own copy of the vectorPotential array,
            ! so this approach might need quite some memory in case the number of
            ! threads and the number of evaluation points is large.

            if (numVertices-1 .lt. actNumProc) then
                nThreads       = numVertices-1
                nSourcePerThread = 1
            else
                nThreads = actNumProc

                ! It is better that many threads do more than one thread needs to do more.
                nSourcePerThread = int(ceiling(real(numVertices-1, kind=wp) / &
                                          nThreads))
            end if

            allocate(vecPotContribs(3, numEvalPos, nThreads), stat=istat)
            if (istat .ne. 0) then
                print *, "could not allocate vecPotContribs: stat=", istat
                return
            end if

            ! parallelized evaluation
!$ifdef _OPENMP
!$omp parallel do private(idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd)
!$endif // _OPENMP
            do idxThread = 0, nThreads-1
                idxSourceStart =      idxThread    * nSourcePerThread + 1
                idxSourceEnd   = min((idxThread+1) * nSourcePerThread, &
                                     numVertices-1)
                idxEvalStart   = 1
                idxEvalEnd     = numEvalPos

                call kernelVectorPotentialPolygonFilament( &
                        vertices, current, &
                        evalPos, &
                        vecPotContribs(:, :, idxThread+1), &
                        idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd, &
                        actUseCompSum)
            end do ! idxThread = 0, nThreads-1

            ! sum up contributions from source chunks
            if (actUseCompSum) then
                do i = 1, numEvalPos
                    ! TODO: bad memory access pattern here --> potential bottleneck !!!
                    sumX(:) = 0.0_wp
                    sumY(:) = 0.0_wp
                    sumZ(:) = 0.0_wp
                    do idxThread = 0, nThreads-1
                        call compAdd(vecPotContribs(1, i, idxThread+1), sumX)
                        call compAdd(vecPotContribs(2, i, idxThread+1), sumY)
                        call compAdd(vecPotContribs(3, i, idxThread+1), sumZ)
                    end do ! idxThread = 0, nThreads-1
                    vectorPotential(1, i) = sum(sumX)
                    vectorPotential(2, i) = sum(sumY)
                    vectorPotential(3, i) = sum(sumZ)
                end do ! i = 1, numEvalPos

                deallocate(vecPotContribs, stat=istat)
                if (istat .ne. 0) then
                    print *, "could not deallocate vecPotContribs: stat=", istat
                    return
                end if
            else
                vectorPotential(:, :) = 0.0_wp
                do idxThread = 0, nThreads-1
                    do i = 1, numEvalPos
                        vectorPotential(:, i) = vectorPotential(:, i) &
                                              + vecPotContribs(:, i, idxThread+1)
                    end do ! i = 1, numEvalPos
                end do ! idxThread = 0, nThreads-1
            end if ! actUseCompSum
        else ! numVertices-1 .gt. numEvalPos
            ! parallelize over nEval

            if (numEvalPos .lt. actNumProc) then
                nThreads       = numEvalPos
                nEvalPerThread = 1
            else
                nThreads = actNumProc
                nEvalPerThread = int(ceiling(real(numEvalPos, kind=wp) / &
                                          nThreads))
            end if

            ! parallelized evaluation
!$ifdef _OPENMP
!$omp parallel do private(idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd)
!$endif // _OPENMP
            do idxThread = 0, nThreads-1
                idxSourceStart = 1
                idxSourceEnd   = numVertices-1
                idxEvalStart   =      idxThread    * nEvalPerThread + 1
                idxEvalEnd     = min((idxThread+1) * nEvalPerThread, &
                                     numEvalPos)

                call kernelVectorPotentialPolygonFilament( &
                        vertices, current, &
                        evalPos, &
                        vectorPotential, &
                        idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd, &
                        actUseCompSum)
            end do ! idxThread = 0, nThreads-1
        end if ! numVertices-1 .gt. numEvalPos
    end if ! actNumProc .eq. 1
end subroutine ! vectorPotentialPolygonFilament

!> Compute the magnetic vector potential of a polygon filament
!> at a number of evaluation locations.
!>
!> @param void (*vertexSupplier)(int i, double *point): callback to put i-th current carrier polygon vertex into point as [3: x, y, z]; in m
!> @param current current along polygon; in A
!> @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
!> @param vectorPotential [3: x, y, z][numEvalPos] target array for magnetic vector potential at evaluation locations; in Tm
!> @param numProcessors number of processors to use for parallelization
!> @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
!>                                of the contributions from the polygon vertices; otherwise, use standard += summation
subroutine vectorPotentialPolygonFilamentVertexSupplier( &
        numVertices, vertexSupplier, current, &
        numEvalPos, evalPos, vectorPotential, &
        numProcessors, useCompensatedSummation)

    interface
        subroutine vertexSupplier(i, pointData)
            use mod_kinds, only: wp => dp
            implicit none
            integer,  intent(in)                :: i
            real(wp), intent(out), dimension(3) :: pointData
        end subroutine
    end interface

    integer,  intent(in)                             :: numVertices
    real(wp), intent(in)                             :: current
    integer,  intent(in)                             :: numEvalPos
    real(wp), intent(in),  dimension(3, numEvalPos)  :: evalPos
    real(wp), intent(out), dimension(3, numEvalPos)  :: vectorPotential
    integer,  intent(in), optional                   :: numProcessors
    logical,  intent(in), optional            :: useCompensatedSummation

    integer :: actNumProc
    logical :: actUseCompSum

    integer :: i, idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd, &
               nThreads, nSourcePerThread, nEvalPerThread, istat, idxThread
    real(wp), dimension(3) :: sumX, sumY, sumZ
    real(wp), dimension(:, :, :), allocatable :: vecPotContribs

    ! handle optional arguments
    if (present(numProcessors)) then
        actNumProc = numProcessors
    else
        ! default: single-threaded
        actNumProc = 1
    end if ! present(numProcessors)

    if (present(useCompensatedSummation)) then
        actUseCompSum = useCompensatedSummation
    else
        ! default: use compensated summation
        actUseCompSum = .true.
    end if

    if (numVertices .lt. 2) then
        print *, "need at least 2 vertices, but got only ", numVertices
        return
    end if

    if (actNumProc .lt. 1) then
        print *, "need at least 1 processor, but got only ", actNumProc
        return
    end if

    if (current .eq. 0.0_wp) then
        vectorPotential(:, :) = 0.0_wp
        return
    end if

    if (actNumProc .eq. 1) then
        ! single-threaded call
        idxSourceStart = 1
        idxSourceEnd   = numVertices-1
        idxEvalStart   = 1
        idxEvalEnd     = numEvalPos
        call kernelVectorPotentialPolygonFilamentVertexSupplier( &
            vertexSupplier, current, &
            evalPos, &
            vectorPotential, &
            idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd, &
            actUseCompSum)
    else
        ! use multithreading

        if (numVertices-1 .gt. numEvalPos) then
            ! parallelize over nSource-1

            ! Note that each thread needs its own copy of the vectorPotential array,
            ! so this approach might need quite some memory in case the number of
            ! threads and the number of evaluation points is large.

            if (numVertices-1 .lt. actNumProc) then
                nThreads       = numVertices-1
                nSourcePerThread = 1
            else
                nThreads = actNumProc

                ! It is better that many threads do more than one thread needs to do more.
                nSourcePerThread = int(ceiling(real(numVertices-1, kind=wp) / &
                                          nThreads))
            end if

            allocate(vecPotContribs(3, numEvalPos, nThreads), stat=istat)
            if (istat .ne. 0) then
                print *, "could not allocate vecPotContribs: stat=", istat
                return
            end if

            ! parallelized evaluation
!$ifdef _OPENMP
!$omp parallel do private(idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd)
!$endif // _OPENMP
            do idxThread = 0, nThreads-1
                idxSourceStart =      idxThread    * nSourcePerThread + 1
                idxSourceEnd   = min((idxThread+1) * nSourcePerThread, &
                                     numVertices-1)
                idxEvalStart   = 1
                idxEvalEnd     = numEvalPos

                call kernelVectorPotentialPolygonFilamentVertexSupplier( &
                        vertexSupplier, current, &
                        evalPos, &
                        vecPotContribs(:, :, idxThread+1), &
                        idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd, &
                        actUseCompSum)
            end do ! idxThread = 0, nThreads-1

            ! sum up contributions from source chunks
            if (actUseCompSum) then
                do i = 1, numEvalPos
                    ! TODO: bad memory access pattern here --> potential bottleneck !!!
                    sumX(:) = 0.0_wp
                    sumY(:) = 0.0_wp
                    sumZ(:) = 0.0_wp
                    do idxThread = 0, nThreads-1
                        call compAdd(vecPotContribs(1, i, idxThread+1), sumX)
                        call compAdd(vecPotContribs(2, i, idxThread+1), sumY)
                        call compAdd(vecPotContribs(3, i, idxThread+1), sumZ)
                    end do ! idxThread = 0, nThreads-1
                    vectorPotential(1, i) = sum(sumX)
                    vectorPotential(2, i) = sum(sumY)
                    vectorPotential(3, i) = sum(sumZ)
                end do ! i = 1, numEvalPos

                deallocate(vecPotContribs, stat=istat)
                if (istat .ne. 0) then
                    print *, "could not deallocate vecPotContribs: stat=", istat
                    return
                end if
            else
                vectorPotential(:, :) = 0.0_wp
                do idxThread = 0, nThreads-1
                    do i = 1, numEvalPos
                        vectorPotential(:, i) = vectorPotential(:, i) &
                                              + vecPotContribs(:, i, idxThread+1)
                    end do ! i = 1, numEvalPos
                end do ! idxThread = 0, nThreads-1
            end if ! actUseCompSum
        else ! numVertices-1 .gt. numEvalPos
            ! parallelize over nEval

            if (numEvalPos .lt. actNumProc) then
                nThreads       = numEvalPos
                nEvalPerThread = 1
            else
                nThreads = actNumProc
                nEvalPerThread = int(ceiling(real(numEvalPos, kind=wp) / &
                                          nThreads))
            end if

            ! parallelized evaluation
!$ifdef _OPENMP
!$omp parallel do private(idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd)
!$endif // _OPENMP
            do idxThread = 0, nThreads-1
                idxSourceStart = 1
                idxSourceEnd   = numVertices-1
                idxEvalStart   =      idxThread    * nEvalPerThread + 1
                idxEvalEnd     = min((idxThread+1) * nEvalPerThread, &
                                     numEvalPos)

                call kernelVectorPotentialPolygonFilamentVertexSupplier( &
                        vertexSupplier, current, &
                        evalPos, &
                        vectorPotential, &
                        idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd, &
                        actUseCompSum)
            end do ! idxThread = 0, nThreads-1
        end if ! numVertices-1 .gt. numEvalPos
    end if ! actNumProc .eq. 1
end subroutine ! vectorPotentialPolygonFilamentVertexSupplier

!> Compute the magnetic field of a polygon filament
!> at a number of evaluation locations.
!>
!> @param vertices [3: x, y, z][numVertices] points along polygon; in m
!> @param current current along polygon; in A
!> @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
!> @param magneticField [3: x, y, z][numEvalPos] target array for magnetic field at evaluation locations; in T
!> @param numProcessors number of processors to use for parallelization
!> @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
!>                                of the contributions from the polygon vertices; otherwise, use standard += summation
subroutine magneticFieldPolygonFilament( &
        numVertices, vertices, current, &
        numEvalPos, evalPos, magneticField, &
        numProcessors, useCompensatedSummation)

    integer,  intent(in)                             :: numVertices
    real(wp), intent(in),  dimension(3, numVertices) :: vertices
    real(wp), intent(in)                             :: current
    integer,  intent(in)                             :: numEvalPos
    real(wp), intent(in),  dimension(3, numEvalPos)  :: evalPos
    real(wp), intent(out), dimension(3, numEvalPos)  :: magneticField
    integer,  intent(in), optional                   :: numProcessors
    logical,  intent(in), optional            :: useCompensatedSummation

    integer :: actNumProc
    logical :: actUseCompSum

    integer :: i, idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd, &
               nThreads, nSourcePerThread, nEvalPerThread, istat, idxThread
    real(wp), dimension(3) :: sumX, sumY, sumZ
    real(wp), dimension(:, :, :), allocatable :: magFldContribs

    ! handle optional arguments
    if (present(numProcessors)) then
        actNumProc = numProcessors
    else
        ! default: single-threaded
        actNumProc = 1
    end if ! present(numProcessors)

    if (present(useCompensatedSummation)) then
        actUseCompSum = useCompensatedSummation
    else
        ! default: use compensated summation
        actUseCompSum = .true.
    end if

    if (numVertices .lt. 2) then
        print *, "need at least 2 vertices, but got only ", numVertices
        return
    end if

    if (actNumProc .lt. 1) then
        print *, "need at least 1 processor, but got only ", actNumProc
        return
    end if

    if (current .eq. 0.0_wp) then
        magneticField(:, :) = 0.0_wp
        return
    end if

    if (actNumProc .eq. 1) then
        ! single-threaded call
        idxSourceStart = 1
        idxSourceEnd   = numVertices-1
        idxEvalStart   = 1
        idxEvalEnd     = numEvalPos
        call kernelMagneticFieldPolygonFilament( &
            vertices, current, &
            evalPos, &
            magneticField, &
            idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd, &
            actUseCompSum)
    else
        ! use multithreading

        if (numVertices-1 .gt. numEvalPos) then
            ! parallelize over nSource-1

            ! Note that each thread needs its own copy of the magneticField array,
            ! so this approach might need quite some memory in case the number of
            ! threads and the number of evaluation points is large.

            if (numVertices-1 .lt. actNumProc) then
                nThreads       = numVertices-1
                nSourcePerThread = 1
            else
                nThreads = actNumProc

                ! It is better that many threads do more than one thread needs to do more.
                nSourcePerThread = int(ceiling(real(numVertices-1, kind=wp) / &
                                          nThreads))
            end if

            allocate(magFldContribs(3, numEvalPos, nThreads), stat=istat)
            if (istat .ne. 0) then
                print *, "could not allocate magFldContribs: stat=", istat
                return
            end if

            ! parallelized evaluation
!$ifdef _OPENMP
!$omp parallel do private(idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd)
!$endif // _OPENMP
            do idxThread = 0, nThreads-1
                idxSourceStart =      idxThread    * nSourcePerThread + 1
                idxSourceEnd   = min((idxThread+1) * nSourcePerThread, &
                                     numVertices-1)
                idxEvalStart   = 1
                idxEvalEnd     = numEvalPos

                call kernelMagneticFieldPolygonFilament( &
                        vertices, current, &
                        evalPos, &
                        magFldContribs(:, :, idxThread+1), &
                        idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd, &
                        actUseCompSum)
            end do ! idxThread = 0, nThreads-1

            ! sum up contributions from source chunks
            if (actUseCompSum) then
                do i = 1, numEvalPos
                    ! TODO: bad memory access pattern here --> potential bottleneck !!!
                    sumX(:) = 0.0_wp
                    sumY(:) = 0.0_wp
                    sumZ(:) = 0.0_wp
                    do idxThread = 0, nThreads-1
                        call compAdd(magFldContribs(1, i, idxThread+1), sumX)
                        call compAdd(magFldContribs(2, i, idxThread+1), sumY)
                        call compAdd(magFldContribs(3, i, idxThread+1), sumZ)
                    end do ! idxThread = 0, nThreads-1
                    magneticField(1, i) = sum(sumX)
                    magneticField(2, i) = sum(sumY)
                    magneticField(3, i) = sum(sumZ)
                end do ! i = 1, numEvalPos

                deallocate(magFldContribs, stat=istat)
                if (istat .ne. 0) then
                    print *, "could not deallocate magFldContribs: stat=", istat
                    return
                end if
            else
                magneticField(:, :) = 0.0_wp
                do idxThread = 0, nThreads-1
                    do i = 1, numEvalPos
                        magneticField(:, i) = magneticField(:, i) &
                                            + magFldContribs(:, i, idxThread+1)
                    end do ! i = 1, numEvalPos
                end do ! idxThread = 0, nThreads-1
            end if ! actUseCompSum
        else ! numVertices-1 .gt. numEvalPos
            ! parallelize over nEval

            if (numEvalPos .lt. actNumProc) then
                nThreads       = numEvalPos
                nEvalPerThread = 1
            else
                nThreads = actNumProc
                nEvalPerThread = int(ceiling(real(numEvalPos, kind=wp) / &
                                          nThreads))
            end if

            ! parallelized evaluation
!$ifdef _OPENMP
!$omp parallel do private(idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd)
!$endif // _OPENMP
            do idxThread = 0, nThreads-1
                idxSourceStart = 1
                idxSourceEnd   = numVertices-1
                idxEvalStart   =      idxThread    * nEvalPerThread + 1
                idxEvalEnd     = min((idxThread+1) * nEvalPerThread, &
                                     numEvalPos)

                call kernelMagneticFieldPolygonFilament( &
                        vertices, current, &
                        evalPos, &
                        magneticField, &
                        idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd, &
                        actUseCompSum)
            end do ! idxThread = 0, nThreads-1
        end if ! numVertices-1 .gt. numEvalPos
    end if ! actNumProc .eq. 1
end subroutine ! magneticFieldPolygonFilament

!> Compute the magnetic field of a polygon filament
!> at a number of evaluation locations.
!>
!> @param void (*vertexSupplier)(int i, double *point): callback to put i-th current carrier polygon vertex into point as [3: x, y, z]; in m
!> @param current current along polygon; in A
!> @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
!> @param magneticField [3: x, y, z][numEvalPos] target array for magnetic field at evaluation locations; in T
!> @param numProcessors number of processors to use for parallelization
!> @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
!>                                of the contributions from the polygon vertices; otherwise, use standard += summation
subroutine magneticFieldPolygonFilamentVertexSupplier( &
        numVertices, vertexSupplier, current, &
        numEvalPos, evalPos, magneticField, &
        numProcessors, useCompensatedSummation)

    interface
        subroutine vertexSupplier(i, pointData)
            use mod_kinds, only: wp => dp
            implicit none
            integer,  intent(in)                :: i
            real(wp), intent(out), dimension(3) :: pointData
        end subroutine
    end interface

    integer,  intent(in)                             :: numVertices
    real(wp), intent(in)                             :: current
    integer,  intent(in)                             :: numEvalPos
    real(wp), intent(in),  dimension(3, numEvalPos)  :: evalPos
    real(wp), intent(out), dimension(3, numEvalPos)  :: magneticField
    integer,  intent(in), optional                   :: numProcessors
    logical,  intent(in), optional            :: useCompensatedSummation

    integer :: actNumProc
    logical :: actUseCompSum

    integer :: i, idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd, &
               nThreads, nSourcePerThread, nEvalPerThread, istat, idxThread
    real(wp), dimension(3) :: sumX, sumY, sumZ
    real(wp), dimension(:, :, :), allocatable :: magFldContribs

    ! handle optional arguments
    if (present(numProcessors)) then
        actNumProc = numProcessors
    else
        ! default: single-threaded
        actNumProc = 1
    end if ! present(numProcessors)

    if (present(useCompensatedSummation)) then
        actUseCompSum = useCompensatedSummation
    else
        ! default: use compensated summation
        actUseCompSum = .true.
    end if

    if (numVertices .lt. 2) then
        print *, "need at least 2 vertices, but got only ", numVertices
        return
    end if

    if (actNumProc .lt. 1) then
        print *, "need at least 1 processor, but got only ", actNumProc
        return
    end if

    if (current .eq. 0.0_wp) then
        magneticField(:, :) = 0.0_wp
        return
    end if

    if (actNumProc .eq. 1) then
        ! single-threaded call
        idxSourceStart = 1
        idxSourceEnd   = numVertices-1
        idxEvalStart   = 1
        idxEvalEnd     = numEvalPos
        call kernelMagneticFieldPolygonFilamentVertexSupplier( &
            vertexSupplier, current, &
            evalPos, &
            magneticField, &
            idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd, &
            actUseCompSum)
    else
        ! use multithreading

        if (numVertices-1 .gt. numEvalPos) then
            ! parallelize over nSource-1

            ! Note that each thread needs its own copy of the vectorPotential array,
            ! so this approach might need quite some memory in case the number of
            ! threads and the number of evaluation points is large.

            if (numVertices-1 .lt. actNumProc) then
                nThreads       = numVertices-1
                nSourcePerThread = 1
            else
                nThreads = actNumProc

                ! It is better that many threads do more than one thread needs to do more.
                nSourcePerThread = int(ceiling(real(numVertices-1, kind=wp) / &
                                          nThreads))
            end if

            allocate(magFldContribs(3, numEvalPos, nThreads), stat=istat)
            if (istat .ne. 0) then
                print *, "could not allocate magFldContribs: stat=", istat
                return
            end if

            ! parallelized evaluation
!$ifdef _OPENMP
!$omp parallel do private(idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd)
!$endif // _OPENMP
            do idxThread = 0, nThreads-1
                idxSourceStart =      idxThread    * nSourcePerThread + 1
                idxSourceEnd   = min((idxThread+1) * nSourcePerThread, &
                                     numVertices-1)
                idxEvalStart   = 1
                idxEvalEnd     = numEvalPos

                call kernelMagneticFieldPolygonFilamentVertexSupplier( &
                        vertexSupplier, current, &
                        evalPos, &
                        magFldContribs(:, :, idxThread+1), &
                        idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd, &
                        actUseCompSum)
            end do ! idxThread = 0, nThreads-1

            ! sum up contributions from source chunks
            if (actUseCompSum) then
                do i = 1, numEvalPos
                    ! TODO: bad memory access pattern here --> potential bottleneck !!!
                    sumX(:) = 0.0_wp
                    sumY(:) = 0.0_wp
                    sumZ(:) = 0.0_wp
                    do idxThread = 0, nThreads-1
                        call compAdd(magFldContribs(1, i, idxThread+1), sumX)
                        call compAdd(magFldContribs(2, i, idxThread+1), sumY)
                        call compAdd(magFldContribs(3, i, idxThread+1), sumZ)
                    end do ! idxThread = 0, nThreads-1
                    magneticField(1, i) = sum(sumX)
                    magneticField(2, i) = sum(sumY)
                    magneticField(3, i) = sum(sumZ)
                end do ! i = 1, numEvalPos

                deallocate(magFldContribs, stat=istat)
                if (istat .ne. 0) then
                    print *, "could not deallocate magFldContribs: stat=", istat
                    return
                end if
            else
                magneticField(:, :) = 0.0_wp
                do idxThread = 0, nThreads-1
                    do i = 1, numEvalPos
                        magneticField(:, i) = magneticField(:, i) &
                                            + magFldContribs(:, i, idxThread+1)
                    end do ! i = 1, numEvalPos
                end do ! idxThread = 0, nThreads-1
            end if ! actUseCompSum
        else ! numVertices-1 .gt. numEvalPos
            ! parallelize over nEval

            if (numEvalPos .lt. actNumProc) then
                nThreads       = numEvalPos
                nEvalPerThread = 1
            else
                nThreads = actNumProc
                nEvalPerThread = int(ceiling(real(numEvalPos, kind=wp) / &
                                          nThreads))
            end if

            ! parallelized evaluation
!$ifdef _OPENMP
!$omp parallel do private(idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd)
!$endif // _OPENMP
            do idxThread = 0, nThreads-1
                idxSourceStart = 1
                idxSourceEnd   = numVertices-1
                idxEvalStart   =      idxThread    * nEvalPerThread + 1
                idxEvalEnd     = min((idxThread+1) * nEvalPerThread, &
                                     numEvalPos)

                call kernelMagneticFieldPolygonFilamentVertexSupplier( &
                        vertexSupplier, current, &
                        evalPos, &
                        magneticField, &
                        idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd, &
                        actUseCompSum)
            end do ! idxThread = 0, nThreads-1
        end if ! numVertices-1 .gt. numEvalPos
    end if ! actNumProc .eq. 1
end subroutine ! magneticFieldPolygonFilamentVertexSupplier

end module ! abscab
