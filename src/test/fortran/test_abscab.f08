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
               rows_B_phi, cols_B_phi

    real(wp), dimension(:,:), allocatable :: rP, zP, ref_A_z, ref_B_phi

    rows_rP = count_rows(filename_rP)
    cols_rP = count_cols(filename_rP)
    if (cols_rP .ne. 1) then
        print *, "error: more than 1 column in ", filename_rP
        return
    end if ! cols_rP .ne. 1
    allocate(rP(cols_rP, rows_rP))
    call read_data(filename_rP, rows_rP, cols_rP, rP)

    rows_zP = count_rows(filename_zP)
    cols_zP = count_cols(filename_zP)
    if (cols_zP .ne. 1) then
        print *, "error: more than 1 column in ", filename_zP
        return
    end if ! cols_zP .ne. 1
    allocate(zP(cols_zP, rows_zP))
    call read_data(filename_zP, rows_zP, cols_zP, zP)

    rows_A_z = count_rows(filename_A_z)
    cols_A_z = count_cols(filename_A_z)
    if (cols_A_z .ne. 1) then
        print *, "error: more than 1 column in ", filename_A_z
        return
    end if ! cols_A_z .ne. 1
    allocate(ref_A_z(cols_A_z, rows_A_z))
    call read_data(filename_A_z, rows_A_z, cols_A_z, ref_A_z)

    rows_B_phi = count_rows(filename_B_phi)
    cols_B_phi = count_cols(filename_B_phi)
    if (cols_B_phi .ne. 1) then
        print *, "error: more than 1 column in ", filename_B_phi
        return
    end if ! cols_B_phi .ne. 1
    allocate(ref_B_phi(cols_B_phi, rows_B_phi))
    call read_data(filename_B_phi, rows_B_phi, cols_B_phi, ref_B_phi)




    deallocate(rP)

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
