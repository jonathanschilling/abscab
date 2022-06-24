program test_cel
use mod_cel
implicit none

    real(wp) :: tolerance = 1.0e-15

    real(wp) :: k_c = 0.1_wp
    real(wp) :: p1 =  4.1_wp
    real(wp) :: p2 = -4.1_wp
    real(wp) :: a = 1.2_wp
    real(wp) :: b = 1.1_wp

    real(wp) :: cel1 =  1.5464442694017956_wp
    real(wp) :: cel2 = -6.7687378198360556e-1_wp

    real(wp) :: c1, c2
    real(wp) :: ra1, ra2

    c1 = cel(k_c, p1, a, b)
    c2 = cel(k_c, p2, a, b)

    ra1 = abs(cel1 - c1)/(1.0_wp + abs(cel1))
    ra2 = abs(cel2 - c2)/(1.0_wp + abs(cel2))

    print *, "case 1: rel/abs deviation = ", ra1
    print *, "case 2: rel/abs deviation = ", ra2

end program test_cel