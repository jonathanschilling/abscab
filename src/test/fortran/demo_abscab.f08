module mod_demo_abscab
use abscab
implicit none

real(wp) :: radius, rCorr

contains

subroutine demo_McGreivy()

    real(wp), parameter :: current = 17.0_wp ! A
    real(wp), dimension(3)   :: center, normal
    real(wp), dimension(3,1) :: evalPos, magneticField

    real(wp) :: bZRef

    radius  = 1.23_wp ! m
    center = (/ 0.0_wp, 0.0_wp, 0.0_wp /)
    normal = (/ 0.0_wp, 0.0_wp, 1.0_wp /)

    evalPos = reshape((/ 10.0_wp, 5.0_wp, 0.0_wp /), shape(evalPos))

    call magneticFieldCircularFilament(center, normal, radius, current, &
        1, evalPos, magneticField)
    bZRef = magneticField(3, 1)
    print *, "ref B_z = ", bZRef

    ! mimic circular wire loop as:
    ! a) Polygon with points on the circule to be mimiced
    ! b) Polygon with points slightly offset radially outward (McGreivy correction)
    ! --> a) should have 2nd-order convergence;
    !     b) should have 4th-order convergence wrt. number of Polygon points








end subroutine ! demo_McGreivy

end module ! mod_demo_abscab

program demo_abscab
    use mod_demo_abscab
    implicit none

    call demo_McGreivy()

end program ! demo_abscab
