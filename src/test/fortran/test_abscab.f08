module mod_abscab_tests
use abscab
use mod_testutil
implicit none
contains

subroutine testStraightWireSegment(status)
    integer, intent(inout) :: status

    character(len=*), parameter :: &
        filename_rP    = "../resources/testPointsRpStraightWireSegment.dat", &
        filename_zP    = "../resources/testPointsZpStraightWireSegment.dat", &
        filename_A_z   = "../resources/StraightWireSegment_A_z_ref.dat", &
        filename_B_phi = "../resources/StraightWireSegment_B_phi_ref.dat"

    real(wp), parameter :: tolerance_A_z   = 1.0e-15_wp, &
                           tolerance_B_phi = 1.0e-15_wp

    integer :: rows_rP,    cols_rP, &
               rows_zP,    cols_zP, &
               rows_A_z,   cols_A_z, &
               rows_B_phi, cols_B_phi, &
               i, numCases, &
               aZStatus, bPhiStatus
    real(wp) :: rP, zP, aZ, bPhi, ref_A_z, ref_B_phi
    real(wp), dimension(:,:), allocatable :: &
        all_rP, all_zP, all_ref_A_z, all_ref_B_phi

    if (status .ne. 0) then
        ! skip if any previous test(s) failed
        return
    end if

    rows_rP = count_rows(filename_rP)
    cols_rP = count_cols(filename_rP)
    if (cols_rP .ne. 1) then
        print *, "error: expecting exactly 1 column in ", filename_rP
        status = 1
        return
    end if ! cols_rP .ne. 1
    numCases = rows_rP
    allocate(all_rP(cols_rP, rows_rP))
    call read_data(filename_rP, rows_rP, cols_rP, all_rP)

    rows_zP = count_rows(filename_zP)
    cols_zP = count_cols(filename_zP)
    if (cols_zP .ne. 1) then
        print *, "error: expecting exactly 1 column in ", filename_zP
        status = 1
        return
    end if ! cols_zP .ne. 1
    if (rows_zP .ne. numCases) then
        print *, "error: number of rows in ", filename_zP, &
                 "does not match number of rows in ", filename_rP
        status = 1
        return
    end if ! rows_zP .ne. numCases
    allocate(all_zP(cols_zP, rows_zP))
    call read_data(filename_zP, rows_zP, cols_zP, all_zP)

    rows_A_z = count_rows(filename_A_z)
    cols_A_z = count_cols(filename_A_z)
    if (cols_A_z .ne. 1) then
        print *, "error: expecting exactly 1 column in ", filename_A_z
        status = 1
        return
    end if ! cols_A_z .ne. 1
    if (rows_A_z .ne. numCases) then
        print *, "error: number of rows in ", filename_A_z, &
                 "does not match number of rows in ", filename_rP
        status = 1
        return
    end if ! rows_A_z .ne. numCases
    allocate(all_ref_A_z(cols_A_z, rows_A_z))
    call read_data(filename_A_z, rows_A_z, cols_A_z, all_ref_A_z)

    rows_B_phi = count_rows(filename_B_phi)
    cols_B_phi = count_cols(filename_B_phi)
    if (cols_B_phi .ne. 1) then
        print *, "error: expecting exactly 1 column in ", filename_B_phi
        status = 1
        return
    end if ! cols_B_phi .ne. 1
    if (rows_B_phi .ne. numCases) then
        print *, "error: number of rows in ", filename_B_phi, &
                 "does not match number of rows in ", filename_rP
        status = 1
        return
    end if ! rows_B_phi .ne. numCases
    allocate(all_ref_B_phi(cols_B_phi, rows_B_phi))
    call read_data(filename_B_phi, rows_B_phi, cols_B_phi, all_ref_B_phi)

    do i = 1, numCases
        rP        = all_rP(1, i)
        zP        = all_zP(1, i)
        ref_A_z   = all_ref_A_z(1, i)
        ref_B_phi = all_ref_B_phi(1, i)

        aZ   = straightWireSegment_A_z(rP, zP)
        bPhi = straightWireSegment_B_phi(rP, zP)

        aZStatus = assertRelAbsEquals(ref_A_z, aZ, tolerance_A_z)
        if (aZStatus .ne. 0) then
            print *, "error: mismatch at Straight Wire Segment A_z test case ", i
            print *, "  rho' = ", rP
            print *, "    z' = ", zP
            print *, "  ref A_z = ", ref_A_z
            print *, "  act A_z = ", aZ
        end if
        status = status + aZStatus

        bPhiStatus = assertRelAbsEquals(ref_B_phi, bPhi, tolerance_B_phi)
        if (bPhiStatus .ne. 0) then
            print *, "error: mismatch at Straight Wire Segment B_phi test case ", i
            print *, "  rho' = ", rP
            print *, "    z' = ", zP
            print *, "  ref B_phi = ", ref_B_phi
            print *, "  act B_phi = ", bPhi
        end if
        status = status + bPhiStatus

        if (status .ne. 0) then
            exit
        end if ! status .ne. 0
    end do ! i = 1, numCases

    deallocate(all_rP)
    deallocate(all_zP)
    deallocate(all_ref_A_z)
    deallocate(all_ref_B_phi)
end subroutine ! testStraightWireSegment

subroutine testCircularWireLoop(status)
    integer, intent(inout) :: status


end subroutine

end module ! mod_abscab_tests

program test_abscab
    use mod_abscab_tests
    implicit none
    integer :: status

    status = 0
    call testStraightWireSegment(status)
    call testCircularWireLoop(status)
    if (status .eq. 0) then
        print *, "success: all test(s) passed :-)"
    else
        print *, "error: some test(s) failed :-("
    end if

end program ! test_abscab
