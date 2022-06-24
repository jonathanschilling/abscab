module abscab
use mod_cel
implicit none

!> vacuum magnetic permeability in Vs/Am (CODATA-2018)
real(wp), parameter :: MU_0 = 1.25663706212e-6_wp;

!> vacuum magnetic permeability, divided by pi
real(wp), parameter :: MU_0_BY_PI = MU_0 / PI;

!> vacuum magnetic permeability, divided by 2 pi
real(wp), parameter :: MU_0_BY_2_PI = MU_0 / (2.0_wp * PI);

!> vacuum magnetic permeability, divided by 4 pi
real(wp), parameter :: MU_0_BY_4_PI = MU_0 / (4.0_wp * PI);

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
    ! Must not use copysign here,
    ! since two negative signs must be able to cancel each other here!
    sws_A_z_ax_n = zP/abs(zP) * log(zP / (zP - 1.0_wp)) / 2.0_wp
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

!> Compute the normalized axial component of the magnetic vector potential of a straight wire segment,
!> evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
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

!> Compute the normalized axial component of the magnetic vector potential of a straight wire segment,
!> evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
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

    Ri_zP    = 2.0_wp * r_i * sinAlphaHalf * sinAlphaHalf ! r_i - z'
    Rf_p_zM1 = 2.0_wp * r_f * sinBetaHalf  * sinBetaHalf  ! r_f - (1 - z')

    n = Ri_zP + Rf_p_zM1

    sws_A_z_n = (log(2.0_wp + n) - log(n)) / 2.0_wp
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






end module ! abscab
