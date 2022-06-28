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

end module ! mod_testutil
