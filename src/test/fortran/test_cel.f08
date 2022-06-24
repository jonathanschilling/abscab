program test_cel
use mod_cel
implicit none

    real(wp), parameter :: tolerance = 1.0e-15

    real(wp), parameter :: k_c = 0.1_wp
    real(wp), parameter :: p1 =  4.1_wp
    real(wp), parameter :: p2 = -4.1_wp
    real(wp), parameter :: a = 1.2_wp
    real(wp), parameter :: b = 1.1_wp

    real(wp), parameter :: cel1 =  1.5464442694017956_wp
    real(wp), parameter :: cel2 = -6.7687378198360556e-1_wp

    real(wp) :: c1,  c2
    real(wp) :: ra1, ra2

    c1 = cel(k_c, p1, a, b)
    c2 = cel(k_c, p2, a, b)

    ra1 = abs(cel1 - c1)/(1.0_wp + abs(cel1))
    ra2 = abs(cel2 - c2)/(1.0_wp + abs(cel2))

    if (ra1 .ge. tolerance) then
        print *, "case 1: rel/abs deviation = ", ra1
        stop 1
    end if

    if (ra2 .ge. tolerance) then
        print *, "case 2: rel/abs deviation = ", ra2
        stop 1
    end if

    print *, "test_cel: all test(s) passed :-)"

end program test_cel
