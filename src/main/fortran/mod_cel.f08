module mod_cel
use mod_kinds, only: wp => dp
use, intrinsic :: ieee_arithmetic
implicit none

! pi
real(wp), parameter :: PI = 4.0_wp * atan(1.0_wp)

! half of pi
real(wp), parameter :: PI_2 = PI / 2.0_wp

! sqrt of machine precision
real(wp), parameter :: sqrt_eps = sqrt(epsilon(1.0_wp))

contains

!> Compute the complete elliptic integral introduced in
!> "Numerical Calculation of Elliptic Integrals and Elliptic Functions. III"
!> by R. Bulirsch in "Numerische Mathematik" 13, 305-315 (1969):
!>
!> cel(k_c, p, a, b) =
!> \int_0^{\pi/2} \frac{a \cos^2{\varphi} + b \sin^2{\varphi}}
!>                     {  \cos^2{\varphi} + p \sin^2{\varphi}}
!>                \frac{\mathrm{d}\varphi}
!>                     {\sqrt{\cos^2{\varphi} + k_c^2 \sin^2{\varphi}}}
!>
!> @param k_c parameter k_c of cel(); absolute value must not be 0
!> @param p   parameter p of cel()
!> @param a   parameter a of cel()
!> @param b   parameter b of cel()
!> @return value of cel(k_c, p, a, b)
function cel(k_c_in, p_in, a_in, b_in)
    real(wp) :: cel
    real(wp) :: k_c_in
    real(wp) :: p_in
    real(wp) :: a_in
    real(wp) :: b_in

    real(wp) :: k_c, p, a, b
    real(wp) :: m, e, f, g, q

    if (k_c_in .eq. 0.0_wp) then
        if (b .eq. 0.0_wp) then
            ! when k_c is zero and b != 0, cel diverges (?)
            cel = ieee_value(1.0_wp, IEEE_POSITIVE_INF)
            return
        else
            k_c = sqrt_eps*sqrt_eps
        end if
    else
        k_c = abs(k_c_in)
    end if
    p = p_in
    a = a_in
    b = b_in

    m = 1.0_wp ! \mu
    e = k_c    ! \nu * \mu
    ! In the iterations, \nu_i is stored in k_c.

    ! initialization
    if (p .gt. 0.0_wp) then
        p = sqrt(p)
        b = b / p
    else ! q <= 0
        f = k_c*k_c        ! f = kc^2 (re-used here; later f = a_i)
        q = 1.0_wp-f       ! 1 - kc^2
        g = 1.0_wp-p
        f = f-p            ! kc^2 - p
        q = q*(b-a*p)      ! (1 - kc^2)*(b-a*p)
        p = sqrt(f/g)      ! sqrt((kc^2 - p)/(1-p))            --> p0
        a = (a-b)/g        ! (a-b)/(1-p)                       --> a0
        b = -q/(g*g*p)+a*p ! -(1-kc^2)*(b-a*p)/( (1-p)^2 * p ) --> b0
    end if

    ! iteration until convergence
    do
        f = a
        a = a + b/p
        g = e/p
        b = b + f*g
        b = b + b
        p = p + g
        g = m
        m = m + k_c
        if (abs(g - k_c) .gt. g * sqrt_eps) then
            k_c = sqrt(e)
            k_c = k_c + k_c
            e = k_c * m
        else
            exit
        end if
    end do

    ! final approximation:
    ! \pi/2 * (a * \mu + b)/(\mu * (\mu + p))
    cel = PI_2 * (a*m + b) / (m*(m+p))
end function ! cel

end module ! mod_cel
