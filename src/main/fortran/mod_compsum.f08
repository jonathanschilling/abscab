module mod_compsum
use mod_kinds, only: wp => dp
implicit none

contains

subroutine compAdd(contribution, compSum)
    real(wp), intent(in)                  :: contribution
    real(wp), intent(inout), dimension(3) :: compSum

    real(wp) :: s, cs, ccs, t, t2, c, cc

    s   = compSum(1)
    cs  = compSum(2)

    t = s + contribution
    if (abs(s) .ge. abs(contribution)) then
        c = (s - t) + contribution
    else
        c = (contribution - t) + s
    end if
    compSum(1) = t

    t2 = cs + c
    if (abs(cs) .ge. abs(c)) then
        cc = (cs - t2) + c
    else
        cc = (c - t2) + cs
    end if
    compSum(2) = t2
    compSum(3) = compSum(3) + cc
end subroutine ! compAdd

end module ! mod_compsum
