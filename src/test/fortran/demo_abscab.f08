module mod_demo_abscab
use abscab
implicit none

real(wp) :: radius, rCorr, omega

contains

subroutine vertexSupplierStd(i, pointData)
    integer,  intent(in)                :: i
    real(wp), intent(out), dimension(3) :: pointData
    real(wp) :: phi
    phi = omega * real(i, kind=wp)
    pointData(1) = radius * cos(phi)
    pointData(2) = radius * sin(phi)
    pointData(3) = 0.0_wp
end subroutine ! vertexSupplierStd

subroutine vertexSupplierMcG(i, pointData)
    integer,  intent(in)                :: i
    real(wp), intent(out), dimension(3) :: pointData
    real(wp) :: phi
    phi = omega * real(i, kind=wp)
    pointData(1) = rCorr * cos(phi)
    pointData(2) = rCorr * sin(phi)
    pointData(3) = 0.0_wp
end subroutine ! vertexSupplierMcG

subroutine demo_McGreivy()
    use mod_testutil, only: errorMetric
!$ifdef _OPENMP
    use omp_lib
!$endif

    real(wp), parameter :: current = 17.0_wp ! A
    real(wp), dimension(3)   :: center, normal
    real(wp), dimension(3,1) :: evalPos, magneticField

    logical :: useCompensatedSummation
    integer  :: i, numCases, numPhi, numProcessors
    real(wp) :: bZRef, bZStd, bZMcG, dPhi

    integer, dimension(*), parameter :: allNumPhi = (/ &
        10, 30, 100, 300, 1000, 3000, &
        10000, 30000, 100000, 300000, &
        1000000, 3000000, &
        10000000, 30000000, &
        100000000, 300000000, 1000000000 /)
    real(wp), dimension(:), allocatable :: &
        allBzStdErr, allBzMcGErr
    real(wp), dimension(:,:), allocatable :: &
        resultTable

    useCompensatedSummation = .true.

!$ifdef _OPENMP
    numProcessors = omp_get_max_threads()
    print *, "using OpenMP"
!$else
!    numProcessors = 1
!$endif
    print *, "numProcessors = ", numProcessors

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
    numCases = size(allNumPhi)
    print *, "number of cases: ", numCases

    allocate(allBzStdErr(numCases), &
             allBzMcGErr(numCases), &
             resultTable(3, numCases))

    do i = 1, numCases
        numPhi = allNumPhi(i)
        print *, "case ",i,"/",numCases," numPhi=",numPhi

        omega = 2.0_wp * PI / real(numPhi - 1, kind=wp)

        call magneticFieldPolygonFilamentVertexSupplier( &
            numPhi, vertexSupplierStd, current, &
            1, evalPos, magneticField, &
            numProcessors, useCompensatedSummation)
        bZStd = magneticField(3,1)

        allBzStdErr(i) = errorMetric(bZRef, bZStd)
        print *, " on-poly B_z: ", bZStd, " (err ", allBzStdErr(i), ")"

        !> McGreivy radius correction
        !> spacing between points
        dPhi = 2.0 * PI / real(numPhi - 1, kind=wp)

        !> TODO: understand derivation of alpha for special case of closed circle
        !> |dr/ds| = 2*pi
        !> --> alpha = 1/R * (dr)^2 / 12
        !> == 4 pi^2 / (12 R)
        rCorr = radius * (1.0_wp + dPhi * dPhi/ 12.0_wp)

        call magneticFieldPolygonFilamentVertexSupplier( &
            numPhi, vertexSupplierMcG, current, &
            1, evalPos, magneticField, &
            numProcessors, useCompensatedSummation)
        bZMcG = magneticField(3,1)

        allBzMcGErr(i) = errorMetric(bZRef, bZMcG)
        print *, "McGreivy B_z: ", bZMcG, " (err ", allBzMcGErr(i), ")"

        resultTable(1, i) = real(numPhi, kind=wp)
        resultTable(2, i) = allBzStdErr(i)
        resultTable(3, i) = allBzMcGErr(i)
    end do ! i = 1, numCases

    ! TODO: dump resultTable to output file

    deallocate(allBzStdErr, allBzMcGErr, resultTable)

end subroutine ! demo_McGreivy

end module ! mod_demo_abscab

program demo_abscab
    use mod_demo_abscab
    implicit none

    call demo_McGreivy()

end program ! demo_abscab
