module abscab
use mod_kinds, only: wp => dp
implicit none

function sws_A_z_ax_f(zP)
    real(wp) :: sws_A_z_ax_f
    real(wp) :: zP

    sws_A_z_ax_f = atanh(1.0_wp / (abs(zP) + abs(1 - zP)))
end function sws_A_z_ax_f

end module abscab
