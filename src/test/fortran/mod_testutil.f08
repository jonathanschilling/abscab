module mod_testutil
use mod_kinds, only: wp => dp
implicit none
contains

function count_rows(filename)
    integer            :: count_rows
    character(len=*) :: filename

    integer, parameter :: iunit = 42
    character(len=1000) :: line
    integer :: rows, istat

    rows = 0
    open(unit=iunit, file=trim(filename), status='old')
    do
        read(iunit, fmt='(e30.25)', iostat=istat)
        if (istat .lt. 0) then
            exit
        else if (istat .gt. 0) then
            print *, "error in reading '"//trim(filename)//"': iostat=", istat
            ! help to debug invalid inputs:
            ! re-read last line that lead to error and print it to screen
            backspace(iunit)
            read(iunit, fmt='(A)') line
            write(*, '(A)') 'Invalid line: '//trim(line)
            exit
        end if
        rows = rows + 1
    end do
    close(unit=iunit)

    count_rows = rows

end function ! count_rows

function count_cols(filename)
    integer            :: count_cols
    character(len=*) :: filename

    integer, parameter :: iunit = 42
    character(len=1000) :: line
    integer :: cols, istat, idxSp, colsThisLine

    cols = 0
    open(unit=iunit, file=trim(filename), status='old')
    do
        read(iunit, fmt='(A)', iostat=istat) line
        if (istat .lt. 0) then
            exit
        else if (istat .gt. 0) then
            print *, "error in reading '"//trim(filename)//"': iostat=", istat
            ! help to debug invalid inputs:
            ! re-read last line that lead to error and print it to screen
            backspace(iunit)
            read(iunit, fmt='(A)') line
            write(*, '(A)') 'Invalid line: '//trim(line)
            exit
        end if

        ! parse columns from line
        idxSp = index(trim(line), " ")
        if (idxSp .eq. 0) then
            ! no space in line --> single entry
            cols = max(cols, 1)
        else
            colsThisLine = 1
            do
                idxSp = index(trim(line(idxSp+1:)), " ")
                if (idxSp .ne. 0) then
                    colsThisLine = colsThisLine + 1
                else
                    cols = max(cols, colsThisLine)
                    exit
                end if
            end do
        end if
    end do
    close(unit=iunit)

    count_cols = cols

end function ! count_cols

subroutine read_data(filename, rows, cols, data)
    character(len=*), intent(in)  :: filename
    integer,            intent(in)  :: rows
    integer,            intent(in)  :: cols
    real(wp), dimension(:,:), intent(out) :: data

    integer, parameter :: iunit = 42
    character(len=1000) :: line
    integer :: r, istat

    open(unit=iunit, file=trim(filename), status='old', iostat=istat)
    if (istat .ne. 0) then
        print *, "error opening file:", istat
    end if
    do r = 1, rows
        read(iunit, fmt='(e30.25)', iostat=istat) data(1, r)
        if (istat .lt. 0) then
            exit
        else if (istat .gt. 0) then
            print *, "error in reading '"//trim(filename)//"': iostat=", istat
            ! help to debug invalid inputs:
            ! re-read last line that lead to error and print it to screen
            backspace(iunit)
            read(iunit, fmt='(A)') line
            write(*, '(A)') 'Invalid line: '//trim(line)
            exit
        end if
    end do
    close(unit=iunit)

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

function errorMetric(ref, act)
    real(wp) :: errorMetric
    real(wp) :: ref, act

    real(wp), parameter :: bad = 0.0_wp, good = -16.0_wp, &
                           tenToBad = 10**bad
    real(wp) :: relErr

    if (abs(ref) .gt. 0.0_wp) then
        if (act .ne. ref) then
            relErr = abs((act-ref)/ref)
            if (tenToBad .lt. relErr) then ! limit to max error of O(1)
                errorMetric = log10(tenToBad)
            else
                errorMetric = log10(relErr)
            end if
        else
            errorMetric = good
        end if
    else if (abs(act) .gt. 0.0_wp) then
        errorMetric = bad
    else
        errorMetric = good
    end if
end function ! errorMetric

end module ! mod_testutil
