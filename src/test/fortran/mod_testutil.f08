module mod_testutil
use mod_kinds, only: wp => dp
implicit none
contains

function count_rows(filename)
    integer            :: count_rows
    character(len=*) :: filename

end function ! count_rows

function count_cols(filename)
    integer            :: count_cols
    character(len=*) :: filename

end function ! count_cols

subroutine read_data(filename, rows, cols, data)
    character(len=*), intent(in)  :: filename
    integer,            intent(in)  :: rows
    integer,            intent(in)  :: cols
    real(wp), dimension(:,:), intent(out) :: data


end subroutine ! read_data

function assertRelAbsEquals(expected, actual, tolerance)
    integer :: assertRelAbsEquals
    real(wp) :: expected, actual, tolerance

    real(wp) :: raErr

    raErr = abs(actual - expected) / (1.0_wp + abs(expected))
    if (raErr .ge. tolerance) then
        print *, "expected:", expected, "actual:", actual, &
          "(rel/abs error", raErr, "tolerance:", tolerance, ")"
        assertRelAbsEquals = 1
    else
        assertRelAbsEquals = 0
    end if

end function ! assertRelAbsEquals

end module ! mod_testutil
