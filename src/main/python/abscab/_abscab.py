import numpy as np

from ._bulirsch_cel import cel
from ._compsum import compAdd

MU_0 = 1.25663706212e-6
"""vacuum magnetic permeability in Vs/Am (CODATA-2018)"""

MU_0_BY_PI = MU_0 / np.pi
r"""vacuum magnetic permeability, divided by :math:`\pi`"""

MU_0_BY_2_PI = MU_0 / (2.0 * np.pi)
r"""vacuum magnetic permeability, divided by :math:`2 \pi`"""

MU_0_BY_4_PI = MU_0 / (4.0 * np.pi)
r"""vacuum magnetic permeability, divided by :math:`4 \pi`"""

############## A_z of straight wire segment

def _sws_A_z_ax_f(zP):
    """Normalized A_z of Straight Wire Segment, along rhoP=0, far-field.

    Compute the normalized axial component of magnetic vector potential of straight wire segment,
    evaluated along axis of wire segment (rho = 0).
    This is a special case for points away from the wire ("far-field") for zP < -1 or zP >= 2.

    :param float zP: normalized vertical coordinate of evaluation location
    :return: normalized axial component of magnetic vector potential
    :rtype: float
    :meta private:
    """
    return np.arctanh(1 / (np.abs(zP) + np.abs(1 - zP)))

def _sws_A_z_ax_n(zP):
    """Normalized A_z of Straight Wire Segment, along rhoP=0, near-field.

    Compute the normalized axial component of magnetic vector potential of straight wire segment,
    evaluated along axis of wire segment (rhoP = 0).
    This is a special case for points close to the wire ("near-field") for -1 <= zP < 2.

    :param float zP: normalized axial coordinate of evaluation location; must not be in [0, 1] (on wire segment)
    :return: normalized axial component of magnetic vector potential
    :rtype: float
    :meta private:
    """

    # Two negative signs must be able to cancel each other here!
    return np.copysign(1.0, zP) * np.log(zP / (zP - 1)) / 2

def _sws_A_z_ax(zP):
    """Normalized A_z of Straight Wire Segment, along rhoP=0.

    Compute the normalized axial component of magnetic vector potential of straight wire segment,
    evaluated along axis of wire segment (rho = 0).

    :param float zP: normalized axial coordinate of evaluation location; must not be in [0, 1] (on wire segment)
    :return: normalized axial component of magnetic vector potential
    :rtype: float
    :meta private:
    """
    if zP < -1 or zP >= 2:
        return _sws_A_z_ax_f(zP)
    else:
        return _sws_A_z_ax_n(zP)

def _sws_A_z_rad_f(rhoP):
    """Normalized A_z of Straight Wire Segment, along zP=0 or zP=1, far-field.

    Compute the normalized axial component of the magnetic vector potential of a straight wire segment,
    evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
    This is a special case for points away from the wire ("far-field") for rhoP > 1.

    :param float rhoP: normalized radial coordinate of evaluation location; must not be zero (on wire segment)
    :return: normalized axial component of magnetic vector potential
    :rtype: float
    :meta private:
    """
    return np.arctanh(1 / (rhoP + np.hypot(rhoP, 1)))

def _sws_A_z_rad_n(rhoP):
    """Normalized A_z of Straight Wire Segment, along zP=0 or zP=1, near-field.

    Compute the normalized axial component of the magnetic vector potential of a straight wire segment,
    evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
    This is a special case for points close to the wire ("near-field") for rhoP <= 1.

    :param float rhoP: normalized radial coordinate of evaluation location; must not be zero (on wire segment)
    :return: normalized axial component of magnetic vector potential
    :rtype: float
    :meta private:
    """
    cat = 1 / np.hypot(rhoP, 1)  # cos(atan(...)  )
    sat = np.sin(np.arctan(rhoP) / 2) # sin(atan(...)/2)
    rc = rhoP * cat
    num = rc + 1 + cat
    den = rc + 2 * sat * sat
    return np.log(num / den) / 2

def _sws_A_z_rad(rhoP):
    """Normalized A_z of Straight Wire Segment, along zP=0 or zP=1.

    Compute the normalized axial component of the magnetic vector potential of a straight wire segment,
    evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).

    :param float rhoP: normalized radial coordinate of evaluation location; must not be zero (on wire segment)
    :return: normalized axial component of magnetic vector potential
    :rtype: float
    :meta private:
    """
    if rhoP > 1:
        return _sws_A_z_rad_f(rhoP)
    else:
        return _sws_A_z_rad_n(rhoP)

def _sws_A_z_f(rhoP, zP):
    """Normalized A_z of Straight Wire Segment, far-field.

    Compute the normalized axial component of the magnetic vector potential of a straight wire segment.
    This formulation is useful for points away from the wire ("far-field")
    at rhoP >= 1 or zP <= -1 or zP > 2.

    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized axial component of magnetic vector potential
    :rtype: float
    :meta private:
    """
    r_i = np.hypot(rhoP, zP)
    r_f = np.hypot(rhoP, 1 - zP)
    return np.arctanh(1 / (r_i + r_f))

def _sws_A_z_n(rhoP, zP):
    """Normalized A_z of Straight Wire Segment, near-field.

    Compute the normalized axial component of the magnetic vector potential of a straight wire segment.
    This formulation is useful for points close to the wire ("near-field")
    at rhoP < 1 and -1 < zP <= 2.

    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized axial component of magnetic vector potential
    :rtype: float
    :meta private:
    """
    omz = 1 - zP

    r_i = np.hypot(rhoP, zP)
    r_f = np.hypot(rhoP, omz)

    alpha = np.arctan2(rhoP, zP)
    sinAlphaHalf = np.sin(alpha / 2)

    beta = np.arctan2(rhoP, omz)
    sinBetaHalf = np.sin(beta / 2)

    Ri_zP    = r_i * sinAlphaHalf * sinAlphaHalf # 0.5 * (r_i - z')
    Rf_p_zM1 = r_f * sinBetaHalf  * sinBetaHalf  # 0.5 * (r_f - (1 - z'))

    n = Ri_zP + Rf_p_zM1

    return (np.log(1 + n) - np.log(n)) / 2

############## B_phi of straight wire segment

def _sws_B_phi_rad(rhoP):
    """Normalized B_phi of Straight Wire Segment, along zP=0 or zP=1.

    Compute the normalized tangential component of the magnetic field of a straight wire segment,
    evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).

    :param float rhoP: normalized radial coordinate of evaluation location
    :return: normalized tangential component of magnetic field
    :rtype: float
    :meta private:
    """
    return 1 / (rhoP * np.hypot(rhoP, 1))

def _sws_B_phi_f(rhoP, zP):
    """Normalized B_phi of Straight Wire Segment, far-field.

    Compute the normalized tangential component of the magnetic field of a straight wire segment.
    This formulation is useful for points away from the wire ("far-field")
    at rhoP >= 1 or zP <= 0 or zP >= 1 or rhoP/(1-zP) >= 1 or rhoP/zP >= 1.

    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized tangential component of magnetic field
    :rtype: float
    :meta private:
    """
    omz = 1 - zP

    r_i = np.hypot(rhoP, zP)
    r_f = np.hypot(rhoP, omz)

    num = rhoP * (1/r_i + 1/r_f)
    den = rhoP * rhoP - zP * omz + r_i * r_f

    return num / den

def _sws_B_phi_n(rhoP, zP):
    """Normalized B_phi of Straight Wire Segment, near-field.

    Compute the normalized tangential component of the magnetic field of a straight wire segment.
    This formulation is useful for points close to the wire ("near-field")
    at rhoP < 1 and 0 < zP < 1 and rhoP/(1-zP) < 1 and rhoP/zP < 1.

    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized tangential component of magnetic field
    :rtype: float
    :meta private:
    """
    omz = 1 - zP

    r_i = np.hypot(rhoP, zP)
    r_f = np.hypot(rhoP, omz)

    num = rhoP * (1/r_i + 1/r_f)

    alpha = np.arctan2(rhoP, zP)
    sinAlphaHalf = np.sin(alpha / 2)

    beta = np.arctan2(rhoP, omz)
    sinBetaHalf = np.sin(beta / 2)

    # r_f * sin^2(beta/2) + (1 - z') * sin^2(alpha/2)
    rfb_omza = r_f * sinBetaHalf * sinBetaHalf + omz * sinAlphaHalf * sinAlphaHalf

    #     r_i * r_f - z' * (1 - z')
    # =   r_i * r_f - r_i * (1 - z') + r_i * (1 - z') - z' * (1 - z')
    # =   r_i * r_f - r_i * r_f * cos(beta)
    #   + r_i * (1 - z') + (1 - z') * r_i * cos(alpha)
    # =   r_i *    r_f   * (1 - cos(beta))
    #   + r_i * (1 - z') * (1 - cos(alpha))
    # = 2 * r_i * [ r_f * sin^2(beta/2) + (1 - z') * sin^2(alpha/2) ]
    den = rhoP * rhoP + 2 * r_i * rfb_omza

    return num / den

############## A_phi of circular wire loop

def _cwl_A_phi_f(rhoP, zP):
    """Normalized A_phi of Circular Wire Loop, far-field.

    Compute the normalized tangential component of the magnetic vector potential of a circular wire loop.
    This formulation is useful for points away from the wire ("far-field")
    at rhoP < 1/2 or rhoP > 2 or :math:`|zP| >= 1`.

    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized tangential component of magnetic vector potential
    :rtype: float
    :meta private:
    """
    sqrt_kCSqNum = np.hypot(zP, 1 - rhoP)
    sqrt_kCSqDen = np.hypot(zP, 1 + rhoP)

    kC = sqrt_kCSqNum / sqrt_kCSqDen
    kSq = 4 * rhoP / (sqrt_kCSqDen * sqrt_kCSqDen)

    kCp1 = 1 + kC
    arg1 = 2 * np.sqrt(kC) / kCp1
    arg2 = 2 / (kCp1 * kCp1 * kCp1)
    C = cel(arg1, 1, 0, arg2)

    return kSq/sqrt_kCSqDen * C

def _cwl_A_phi_n(rhoP, zP):
    """Normalized A_phi of Circular Wire Loop, near-field.

    Compute the normalized tangential component of the magnetic vector potential of a circular wire loop.
    This formulation is useful for points close to the wire ("near-field")
    at 1/2 <= rhoP <= 2 and :math:`|zP| < 1`.

    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized tangential component of magnetic vector potential
    :rtype: float
    :meta private:
    """
    rhoP_m_1 = rhoP - 1

    n = zP / rhoP_m_1
    m = 1 + 2 / rhoP_m_1

    num = n * n + 1
    den = n * n + m * m

    kCSq = num / den

    prefac = 1 / (np.abs(rhoP - 1) * np.sqrt(den))
    celPart = cel(np.sqrt(kCSq), 1, -1, 1)
    return prefac * celPart

def _cwl_A_phi_v(zP):
    """Normalized A_phi of Circular Wire Loop, along rhoP=1, near-field.

    Compute the normalized tangential component of the magnetic vector potential of a circular wire loop.
    This formulation is useful for points along rhoP=1 with :math:`|zP| < 1`.

    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized tangential component of magnetic vector potential
    :rtype: float
    :meta private:
    """
    absZp = np.abs(zP)

    # 1/k_c
    kCInv = np.sqrt(4 + zP * zP) / absZp

    return cel(kCInv, 1, 1, -1) / absZp

############## B_rho of circular wire loop

def _cwl_B_rho_f(rhoP, zP):
    """Normalized B_rho of Circular Wire Loop, far-field.

    Compute the normalized radial component of the magnetic field of a circular wire loop.
    This formulation is useful for points away from the wire ("far-field")
    at rhoP < 1/2 or rhoP > 2 or :math:`|zP| >= 1`.

    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized radial component of magnetic field
    :rtype: float
    :meta private:
    """
    sqrt_kCSqNum = np.hypot(zP, 1 - rhoP)
    sqrt_kCSqDen = np.hypot(zP, 1 + rhoP)

    kCSqNum = sqrt_kCSqNum * sqrt_kCSqNum
    kCSqDen = sqrt_kCSqDen * sqrt_kCSqDen

    kCSq = kCSqNum / kCSqDen
    kC = np.sqrt(kCSq)

    D = cel(kC, 1, 0, 1)

    kCp1 = 1 + kC
    arg1 = 2 * np.sqrt(kC) / kCp1
    arg2 = 2 / (kCp1 * kCp1 * kCp1)
    C = cel(arg1, 1, 0, arg2)

    prefac = 4 * rhoP / (kCSqDen * sqrt_kCSqDen * kCSqNum)

    return prefac * zP * (D - C)

def _cwl_B_rho_n(rhoP, zP):
    """Normalized B_rho of Circular Wire Loop, near-field.

    Compute the normalized radial component of the magnetic field of a circular wire loop.
    This formulation is useful for points close to the wire ("near-field")
    at 1/2 <= rhoP <= 2 and :math:`|zP| < 1`.

    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized radial component of magnetic field
    :rtype: float
    :meta private:
    """
    rhoP_m_1 = rhoP - 1
    rd2 = rhoP_m_1 * rhoP_m_1

    n = zP / rhoP_m_1
    m = 1 + 2 / rhoP_m_1

    sqrt_kCSqNum = np.hypot(n, 1)
    sqrt_kCSqDen = np.hypot(n, m)

    kCSqNum = sqrt_kCSqNum * sqrt_kCSqNum
    kCSqDen = sqrt_kCSqDen * sqrt_kCSqDen

    kC = sqrt_kCSqNum / sqrt_kCSqDen

    D = cel(kC, 1, 0, 1)

    kCp1 = 1 + kC
    arg1 = 2 * np.sqrt(kC) / kCp1
    arg2 = 2 / (kCp1 * kCp1 * kCp1)
    C = arg2 * cel(arg1, 1, 0, 1)

    # z' / |rho' - 1|^5
    zP_rd5 = zP / (np.abs(rhoP_m_1) * rd2 * rd2)

    prefac = 4 * rhoP / (kCSqDen * sqrt_kCSqDen * kCSqNum)

    return prefac * zP_rd5 * (D - C)

def _cwl_B_rho_v(zP):
    """Normalized B_rho of Circular Wire Loop, along rhoP=1, near-field.

    Compute the normalized tangential component of the magnetic vector potential of a circular wire loop.
    This formulation is useful for points along rhoP=1 with :math:`|zP| < 1`.

    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized radial component of magnetic field
    :rtype: float
    :meta private:
    """
    zPSq = zP * zP

    kCSq = 1 / (1 + 4 / zPSq)
    kC = np.sqrt(kCSq)

    K = cel(kC, 1, 1, 1)
    E = cel(kC, 1, 1, kCSq)

    return np.copysign(kC / 2 * ((2 / zPSq + 1) * E - K), zP)

############## B_z of circular wire loop

def _cwl_B_z_f1(rhoP, zP):
    """Normalized B_z of Circular Wire Loop, far-field (1).

    Compute the normalized vertical component of the magnetic field of a circular wire loop.
    This formulation is useful for certain points away from the wire ("far-field")
    at rhoP < 1/2 or (rhoP <= 2 and :math:`|zP| >= 1`).

    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized vertical component of magnetic field
    :rtype: float
    :meta private:
    """
    sqrt_kCSqNum = np.hypot(zP, 1 - rhoP)
    sqrt_kCSqDen = np.hypot(zP, 1 + rhoP)

    kC = sqrt_kCSqNum / sqrt_kCSqDen

    K = cel(kC, 1, 1, 1)
    E = cel(kC, 1, 1, kC * kC)
    D = cel(kC, 1, 0, 1)

    prefac = 1 / (sqrt_kCSqDen * sqrt_kCSqNum * sqrt_kCSqNum)
    comb = (E - 2 * K + 2 * D)

    return prefac * (E + rhoP * comb)

def _cwl_B_z_f2(rhoP, zP):
    """Normalized B_z of Circular Wire Loop, far-field (2).

    Compute the normalized vertical component of the magnetic field of a circular wire loop.
    This formulation is useful for certain other points away from the wire ("far-field")
    at rhoP > 2.

    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized vertical component of magnetic field
    :rtype: float
    :meta private:
    """
    sqrt_kCSqNum = np.hypot(zP, 1 - rhoP)
    sqrt_kCSqDen = np.hypot(zP, 1 + rhoP)

    kC = sqrt_kCSqNum / sqrt_kCSqDen
    kCSq = kC * kC

    zPSqP1 = zP * zP + 1
    rhoPSq = rhoP * rhoP
    t1 = zPSqP1 / rhoPSq + 1
    t2 = 2 / rhoP

    # a is sqrt_kCSqDen normalized to rho'^2
    # b is sqrt_kCSqNum normalized to rho'^2
    # a == (z'^2 + (1 + rho')^2) / rho'^2 = (z'^2 + 1)/rho'^2 + 1  +  2/rho'
    # b == (z'^2 + (1 - rho')^2) / rho'^2 = (z'^2 + 1)/rho'^2 + 1  -  2/rho'
    a = t1 + t2
    b = t1 - t2

    # 1/prefac = sqrt( z'^2 + (1 + rho')^2)           * (z'^2 + (1 - rho')^2)
    #          = sqrt((z'^2 + (1 + rho')^2) / rho'^2) * (z'^2 + (1 - rho')^2) / rho'^2 * rho'^3
    #          = sqrt(a)                              * b                              * rho'^3
    prefac = 1 / (np.sqrt(a) * b * rhoPSq * rhoP)

    cdScale = 1 + (2 + zPSqP1 / rhoP) / rhoP

    E = cel(kC, 1, 1, kCSq)
    D = cel(kC, 1, 0, 1)

    kCP1 = 1 + kC
    arg1 = 2 * np.sqrt(kC) / kCP1
    arg2 = 2 / (kCP1 * kCP1 * kCP1)
    C = arg2 * cel(arg1, 1, 0, 1)

    # use C - D for (2 * D - E)/kSq
    return prefac * (E + 4 * (C - D) / cdScale)

def _cwl_B_z_n(rhoP, zP):
    """Normalized B_z of Circular Wire Loop, near-field.

    Compute the normalized vertical component of the magnetic field of a circular wire loop.
    This formulation is useful for points close to the wire ("near-field")
    at 1/2 <= rhoP <= 2, but not rhoP=1, and :math:`|zP| <= 1`.

    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized vertical component of magnetic field
    :rtype: float
    :meta private:
    """
    rp1 = rhoP - 1

    n = zP / rp1
    m = 1 + 2 / rp1

    sqrt_kCSqNum = np.hypot(n, 1)
    sqrt_kCSqDen = np.hypot(n, m)

    kCSqDen = sqrt_kCSqDen * sqrt_kCSqDen

    kC = sqrt_kCSqNum / sqrt_kCSqDen

    prefac = 1 / (np.abs(rp1) * rp1 * rp1 * kCSqDen * sqrt_kCSqDen)

    return prefac * cel(kC, kC * kC, 1 + rhoP, 1 - rhoP)

def _cwl_B_z_v(zP):
    """Normalized B_z of Circular Wire Loop, along rhoP=1, near-field.

    Compute the normalized tangential component of the magnetic vector potential of a circular wire loop.
    This formulation is useful for points along rhoP=1 with :math:`|zP| < 1`.

    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized radial component of magnetic field
    :rtype: float
    :meta private:
    """
    kCSq = zP * zP / (4 + zP * zP)
    kC = np.sqrt(kCSq)

    f = zP * zP + 4
    prefac = 1 / (f * np.sqrt(f))

    return prefac * cel(kC, kCSq, 2, 0)

# --------------------------------------------------

def straightWireSegment_A_z(rhoP, zP):
    """Compute the normalized axial component of the magnetic vector potential of a straight wire segment.

    This method selects the proper special case method to use
    for computing the normalized axial component of the magnetic vector potential
    of a straight wire segment for given evaluation location (rhoP, zP).

    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized axial component of magnetic vector potential
    :rtype: float
    :meta private:
    """
    if rhoP == 0.0:
        if zP < 0.0 or zP > 1.0:
            return _sws_A_z_ax(zP)
        else:
            raise RuntimeError("evaluation locations on the wire segment (rho'=%g z'=%g) are not allowed"%(rhoP, zP))
    elif zP == 0.0 or zP == 1.0:
        return _sws_A_z_rad(rhoP)
    elif rhoP >= 1.0 or zP <= -1.0 or zP > 2.0:
        return _sws_A_z_f(rhoP, zP)
    else:
        return _sws_A_z_n(rhoP, zP)

def straightWireSegment_B_phi(rhoP, zP):
    """Compute the normalized tangential component of the magnetic field of a straight wire segment.

    This method selects the proper special case method to use
    for computing the normalized tangential component of the magnetic field
    of a straight wire segment for given evaluation location (rhoP, zP).

    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized tangential component of magnetic field
    :rtype: float
    :meta private:
    """
    if rhoP == 0.0:
        if zP < 0.0 or zP > 1.0:
            return 0.0
        else:
            raise RuntimeError("evaluation locations on the wire segment (rho'=%g z'=%g) are not allowed"%(rhoP, zP))
    elif zP == 0.0 or zP == 1.0:
        return _sws_B_phi_rad(rhoP)
    elif rhoP >= zP or rhoP >= 1.0 - zP or zP < 0.0 or zP > 1.0:
        return _sws_B_phi_f(rhoP, zP)
    else:
        return _sws_B_phi_n(rhoP, zP)

def circularWireLoop_A_phi(rhoP, zP):
    """Compute the normalized tangential component of the magnetic vector potential of a circular wire loop.

    This method selects the proper special case method to use
    for computing the normalized tangential component of the magnetic vector potential
    of a circular wire loop for given evaluation location (rhoP, zP).

    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized tangential component of magnetic vector potential
    :rtype: float
    :meta private:
    """
    if rhoP == 0.0:
        return 0.0
    elif rhoP < 0.5 or rhoP > 2.0 or np.abs(zP) >= 1.0:
        return _cwl_A_phi_f(rhoP, zP)
    elif rhoP != 1.0:
        return _cwl_A_phi_n(rhoP, zP)
    else:
        if zP != 0.0:
            return _cwl_A_phi_v(zP)
        else:
            raise RuntimeError("evaluation at location of wire loop (rho' = 1, z' = 0) is not defined")

def circularWireLoop_B_rho(rhoP, zP):
    """Compute the normalized radial component of the magnetic field of a circular wire loop.

    This method selects the proper special case method to use
    for computing the normalized radial component of the magnetic field
    of a circular wire loop for given evaluation location (rhoP, zP).

    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized radial component of magnetic field
    :rtype: float
    :meta private:
    """
    if rhoP == 0.0 or zP == 0.0:
        if rhoP != 1.0:
            return 0.0
        else :
            raise RuntimeError("evaluation at location of wire loop (rho' = 1, z' = 0) is not defined")
    elif rhoP < 0.5 or rhoP > 2.0 or np.abs(zP) >= 1.0:
        return _cwl_B_rho_f(rhoP, zP)
    elif rhoP != 1.0:
        return _cwl_B_rho_n(rhoP, zP)
    else:
        return _cwl_B_rho_v(zP)

def circularWireLoop_B_z(rhoP, zP):
    """Compute the normalized vertical component of the magnetic field of a circular wire loop.

    This method selects the proper special case method to use
    for computing the normalized vertical component of the magnetic field
    of a circular wire loop for given evaluation location (rhoP, zP).

    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized radial component of magnetic field
    :rtype: float
    :meta private:
    """
    if rhoP < 0.5 or (rhoP <= 2 and np.abs(zP) > 1):
        return _cwl_B_z_f1(rhoP, zP)
    elif rhoP > 2:
        return _cwl_B_z_f2(rhoP, zP)
    elif rhoP != 1.0:
        return _cwl_B_z_n(rhoP, zP)
    else:
        if zP != 0.0:
            return _cwl_B_z_v(zP)
        else:
            raise RuntimeError("evaluation at location of wire loop (rho' = 1, z' = 0) is not defined")

# --------------------------------------------------

def vectorPotentialCircularFilament(center, normal, radius, current, evalPos):
    """Compute the magnetic vector potential of a circular wire loop.

    :param arr(float) center: [3: x, y, z] origin of loop (in meters)
    :param arr(float) normal: [3: x, y, z]  normal vector of loop (in meters); will be normalized internally
    :param float radius: radius of the wire loop (in meters)
    :param float current: loop current (in A)
    :param arr(float) evalPos: [nEvalPos][3: x, y, z] evaluation locations (in meters)
    :return: [nEvalPos][3: A_x, A_y, A_z] Cartesian components of the magnetic
             vector potential evaluated at the given locations (in Tm)
    :rtype: arr(float)
    """
    if not np.isfinite(radius) or radius <= 0.0:
        raise ValueError("radius must be finite and positive, but is %g"%(radius,))

    aPrefactor = MU_0_BY_PI * current

    # squared length of normal vector
    nLen2 = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]

    if nLen2 == 0.0:
        raise ValueError("length of normal vector must not be zero")

    # length of normal vector
    nLen = np.sqrt(nLen2)

    # unit normal vector of wire loop
    eX = normal[0] / nLen
    eY = normal[1] / nLen
    eZ = normal[2] / nLen

    nEvalPos = len(evalPos)

    vectorPotential = np.zeros((nEvalPos, 3))

    for idxEval in range(nEvalPos):

        # vector from center of wire loop to eval pos
        r0x = evalPos[idxEval, 0] - center[0]
        r0y = evalPos[idxEval, 1] - center[1]
        r0z = evalPos[idxEval, 2] - center[2]

        # z position along normal of wire loop
        alignedZ = eX * r0x + eY * r0y + eZ * r0z

        # normalized z component of evaluation location in coordinate system of wire loop
        zP = alignedZ / radius

        # r0 projected onto axis of wire loop
        rParallelX = alignedZ * eX
        rParallelY = alignedZ * eY
        rParallelZ = alignedZ * eZ

        # vector perpendicular to axis of wire loop, pointing at evaluation pos
        rPerpX = r0x - rParallelX
        rPerpY = r0y - rParallelY
        rPerpZ = r0z - rParallelZ

        # perpendicular distance squared between evalPos and axis of wire loop
        alignedRSq = rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ

        # prevent division-by-zero when computing radial unit vector
        # A_phi is zero anyway on-axis --> no contribution expected
        if alignedRSq > 0.0:

            # perpendicular distance between evalPos and axis of wire loop
            alignedR = np.sqrt(alignedRSq)

            # unit vector in radial direction
            eRX = rPerpX / alignedR
            eRY = rPerpY / alignedR
            eRZ = rPerpZ / alignedR

            # normalized rho component of evaluation location in coordinate system of wire loop
            rhoP = alignedR / radius

            # compute tangential component of magnetic vector potential, including current and mu_0
            aPhi = aPrefactor * circularWireLoop_A_phi(rhoP, zP)

            # compute cross product between e_z and e_rho to get e_phi
            ePhiX = eRY * eZ - eRZ * eY
            ePhiY = eRZ * eX - eRX * eZ
            ePhiZ = eRX * eY - eRY * eX

            # add contribution from wire loop to result
            vectorPotential[idxEval, 0] = aPhi * ePhiX
            vectorPotential[idxEval, 1] = aPhi * ePhiY
            vectorPotential[idxEval, 2] = aPhi * ePhiZ

    return vectorPotential

def magneticFieldCircularFilament(center, normal, radius, current, evalPos):
    """Compute the magnetic field of a circular wire loop.

    :param arr(float) center: [3: x, y, z] origin of loop (in meters)
    :param arr(float) normal: [3: x, y, z]  normal vector of loop (in meters); will be normalized internally
    :param float radius: radius of the wire loop (in meters)
    :param float current: loop current (in A)
    :param arr(float) evalPos: [nEvalPos][3: x, y, z] evaluation locations (in meters)
    :return: [nEvalPos][3: B_x, B_y, B_z] Cartesian components of the magnetic
             field evaluated at the given locations (in T)
    :rtype: arr(float)
    """
    if not np.isfinite(radius) or radius <= 0.0:
        raise ValueError("radius must be finite and positive, but is %g"%(radius,))

    bPrefactor = MU_0_BY_PI * current / radius

    # squared length of normal vector
    nLen2 = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]

    if nLen2 == 0.0:
        raise ValueError("length of normal vector must not be zero")

    # length of normal vector
    nLen = np.sqrt(nLen2)

    # unit normal vector of wire loop
    eX = normal[0] / nLen
    eY = normal[1] / nLen
    eZ = normal[2] / nLen

    nEvalPos = len(evalPos)

    magneticField = np.zeros((nEvalPos, 3))

    for idxEval in range(nEvalPos):

        # vector from center of wire loop to eval pos
        r0x = evalPos[idxEval, 0] - center[0]
        r0y = evalPos[idxEval, 1] - center[1]
        r0z = evalPos[idxEval, 2] - center[2]

        # z position along normal of wire loop
        alignedZ = eX * r0x + eY * r0y + eZ * r0z

        # normalized z component of evaluation location in coordinate system of wire loop
        zP = alignedZ / radius

        # r0 projected onto axis of wire loop
        rParallelX = alignedZ * eX
        rParallelY = alignedZ * eY
        rParallelZ = alignedZ * eZ

        # vector perpendicular to axis of wire loop, pointing at evaluation pos
        rPerpX = r0x - rParallelX
        rPerpY = r0y - rParallelY
        rPerpZ = r0z - rParallelZ

        # perpendicular distance squared between evalPos and axis of wire loop
        alignedRSq = rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ

        if alignedRSq > 0.0:
            # radial unit vector is only defined if evaluation pos is off-axis

            # perpendicular distance between evalPos and axis of wire loop
            alignedR = np.sqrt(alignedRSq)

            # unit vector in radial direction
            eRX = rPerpX / alignedR
            eRY = rPerpY / alignedR
            eRZ = rPerpZ / alignedR

            # normalized rho component of evaluation location in coordinate system of wire loop
            rhoP = alignedR / radius

            # compute radial component of normalized magnetic field
            # and scale by current and mu_0
            bRho = bPrefactor * circularWireLoop_B_rho(rhoP, zP)

            # add contribution from B_rho of wire loop to result
            magneticField[idxEval, 0] = bRho * eRX
            magneticField[idxEval, 1] = bRho * eRY
            magneticField[idxEval, 2] = bRho * eRZ
        else:
            rhoP = 0.0

        # compute vertical component of normalized magnetic field
        # and scale by current and mu_0
        bZ = bPrefactor * circularWireLoop_B_z(rhoP, zP)

        # add contribution from B_z of wire loop to result
        magneticField[idxEval, 0] += bZ * eX
        magneticField[idxEval, 1] += bZ * eY
        magneticField[idxEval, 2] += bZ * eZ

    return magneticField

# --------------------------------------------------

def vectorPotentialPolygonFilament(vertices, current, evalPos, useCompensatedSummation = True):
    """Compute the magnetic vector potential of a polygon filament at a number of evaluation locations.

    :param arr(float) vertices: [numVertices][3: x, y, z] points along polygon; in m
    :param float current: current along polygon; in A
    :param arr(float) evalPos: [numEvalPos][3: x, y, z] evaluation locations; in m
    :param bool useCompensatedSummation: If true, use Kahan-Babuska compensated summation to compute the superposition
                                         of the contributions from the polygon vertices; otherwise, use standard += summation.
    :return: [numEvalPos][3: x, y, z] magnetic vector potential at evaluation locations; in Tm
    :rtype: arr(float)
    """
    numEvalPos = len(evalPos)
    numVertices = len(vertices)
    if numVertices < 2:
        raise ValueError("must have at least 2 vertices, not %d"%(numVertices,))

    vectorPotential = np.zeros((numEvalPos, 3))

    if current == 0.0:
        return vectorPotential

    aPrefactor = MU_0_BY_2_PI * current

    # setup compensated summation objects
    if useCompensatedSummation:
        # need three doubles (s, cs, ccs) per eval pos --> see compsum.py
        aXSum = np.zeros((numEvalPos, 3))
        aYSum = np.zeros((numEvalPos, 3))
        aZSum = np.zeros((numEvalPos, 3))

    x_i = vertices[0, 0]
    y_i = vertices[0, 1]
    z_i = vertices[0, 2]

    for idxSource in range(numVertices-1):

        x_f = vertices[idxSource + 1, 0]
        y_f = vertices[idxSource + 1, 1]
        z_f = vertices[idxSource + 1, 2]

        # vector from start to end of i:th wire segment
        dx = x_f - x_i
        dy = y_f - y_i
        dz = z_f - z_i

        # squared length of wire segment
        l2 = dx * dx + dy * dy + dz * dz
        if l2 == 0.0:
            # skip zero-length segments: no contribution
            continue

        # length of wire segment
        l = np.sqrt(l2)

        # unit vector parallel to wire segment
        eX = dx / l
        eY = dy / l
        eZ = dz / l

        for idxEval in range(numEvalPos):

            # vector from start of wire segment to eval pos
            r0x = evalPos[idxEval, 0] - x_i
            r0y = evalPos[idxEval, 1] - y_i
            r0z = evalPos[idxEval, 2] - z_i

            # z position along axis of wire segment
            alignedZ = eX * r0x + eY * r0y + eZ * r0z

            # normalized z component of evaluation location in coordinate system of wire segment
            zP = alignedZ / l

            # vector perpendicular to axis of wire segment, pointing at evaluation pos
            rPerpX = r0x - alignedZ * eX
            rPerpY = r0y - alignedZ * eY
            rPerpZ = r0z - alignedZ * eZ

            # perpendicular distance between evalPos and axis of wire segment
            alignedR = np.sqrt(rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ)

            # normalized rho component of evaluation location in coordinate system of wire segment
            rhoP = alignedR / l

            # compute parallel component of magnetic vector potential, including current and mu_0
            aParallel = aPrefactor * straightWireSegment_A_z(rhoP, zP)

            # add contribution from wire segment to result
            if useCompensatedSummation:
                compAdd(aParallel * eX, aXSum[idxEval, :])
                compAdd(aParallel * eY, aYSum[idxEval, :])
                compAdd(aParallel * eZ, aZSum[idxEval, :])
            else:
                vectorPotential[idxEval, 0] += aParallel * eX
                vectorPotential[idxEval, 1] += aParallel * eY
                vectorPotential[idxEval, 2] += aParallel * eZ

        # shift to next point
        x_i = x_f
        y_i = y_f
        z_i = z_f

    if useCompensatedSummation:
        # obtain compensated sums from summation objects
        for idxEval in range(numEvalPos):
            vectorPotential[idxEval, 0] = aXSum[idxEval, 0] + aXSum[idxEval, 1] + aXSum[idxEval, 2]
            vectorPotential[idxEval, 1] = aYSum[idxEval, 0] + aYSum[idxEval, 1] + aYSum[idxEval, 2]
            vectorPotential[idxEval, 2] = aZSum[idxEval, 0] + aZSum[idxEval, 1] + aZSum[idxEval, 2]

    return vectorPotential

def vectorPotentialPolygonFilamentVertexSupplier(numVertices, vertexSupplier, current, evalPos, useCompensatedSummation = True):
    """Compute the magnetic vector potential of a polygon filament at a number of evaluation locations.

    :param int numVertices: number of polygon vertices to take into account
    :param callable(int i), arr(float) vertexSupplier: should return points along polygon as [3: x, y, z] for i=0,1,...,(numVertices-1); in m
    :param float current: current along polygon; in A
    :param arr(float) evalPos: [numEvalPos][3: x, y, z] evaluation locations; in m
    :param bool useCompensatedSummation: If true, use Kahan-Babuska compensated summation to compute the superposition
                                         of the contributions from the polygon vertices; otherwise, use standard += summation.
    :return: [numEvalPos][3: x, y, z] magnetic vector potential at evaluation locations; in Tm
    :rtype: arr(float)
    """
    numEvalPos = len(evalPos)
    if numVertices < 2:
        raise ValueError("must have at least 2 vertices, not %d"%(numVertices,))

    vectorPotential = np.zeros((numEvalPos, 3))

    if current == 0.0:
        return vectorPotential

    aPrefactor = MU_0_BY_2_PI * current

    # setup compensated summation objects
    if useCompensatedSummation:
        # need three doubles (s, cs, ccs) per eval pos --> see compsum.py
        aXSum = np.zeros((numEvalPos, 3))
        aYSum = np.zeros((numEvalPos, 3))
        aZSum = np.zeros((numEvalPos, 3))

    pointData = vertexSupplier(0)
    x_i = pointData[0]
    y_i = pointData[1]
    z_i = pointData[2]

    for idxSource in range(numVertices-1):

        pointData = vertexSupplier(idxSource + 1)
        x_f = pointData[0]
        y_f = pointData[1]
        z_f = pointData[2]

        # vector from start to end of i:th wire segment
        dx = x_f - x_i
        dy = y_f - y_i
        dz = z_f - z_i

        # squared length of wire segment
        l2 = dx * dx + dy * dy + dz * dz
        if l2 == 0.0:
            # skip zero-length segments: no contribution
            continue

        # length of wire segment
        l = np.sqrt(l2)

        # unit vector parallel to wire segment
        eX = dx / l
        eY = dy / l
        eZ = dz / l

        for idxEval in range(numEvalPos):

            # vector from start of wire segment to eval pos
            r0x = evalPos[idxEval, 0] - x_i
            r0y = evalPos[idxEval, 1] - y_i
            r0z = evalPos[idxEval, 2] - z_i

            # z position along axis of wire segment
            alignedZ = eX * r0x + eY * r0y + eZ * r0z

            # normalized z component of evaluation location in coordinate system of wire segment
            zP = alignedZ / l

            # vector perpendicular to axis of wire segment, pointing at evaluation pos
            rPerpX = r0x - alignedZ * eX
            rPerpY = r0y - alignedZ * eY
            rPerpZ = r0z - alignedZ * eZ

            # perpendicular distance between evalPos and axis of wire segment
            alignedR = np.sqrt(rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ)

            # normalized rho component of evaluation location in coordinate system of wire segment
            rhoP = alignedR / l

            # compute parallel component of magnetic vector potential, including current and mu_0
            aParallel = aPrefactor * straightWireSegment_A_z(rhoP, zP)

            # add contribution from wire segment to result
            if useCompensatedSummation:
                compAdd(aParallel * eX, aXSum[idxEval, :])
                compAdd(aParallel * eY, aYSum[idxEval, :])
                compAdd(aParallel * eZ, aZSum[idxEval, :])
            else:
                vectorPotential[idxEval, 0] += aParallel * eX
                vectorPotential[idxEval, 1] += aParallel * eY
                vectorPotential[idxEval, 2] += aParallel * eZ

        # shift to next point
        x_i = x_f
        y_i = y_f
        z_i = z_f

    if useCompensatedSummation:
        # obtain compensated sums from summation objects
        for idxEval in range(numEvalPos):
            vectorPotential[idxEval, 0] = aXSum[idxEval, 0] + aXSum[idxEval, 1] + aXSum[idxEval, 2]
            vectorPotential[idxEval, 1] = aYSum[idxEval, 0] + aYSum[idxEval, 1] + aYSum[idxEval, 2]
            vectorPotential[idxEval, 2] = aZSum[idxEval, 0] + aZSum[idxEval, 1] + aZSum[idxEval, 2]

    return vectorPotential

def magneticFieldPolygonFilament(vertices, current, evalPos, useCompensatedSummation = True):
    """Compute the magnetic field of a polygon filament at a number of evaluation locations.

    :param arr(float) vertices: [numVertices][3: x, y, z] points along polygon; in m
    :param float current: current along polygon; in A
    :param arr(float) evalPos: [numEvalPos][3: x, y, z] evaluation locations; in m
    :param bool useCompensatedSummation: If true, use Kahan-Babuska compensated summation to compute the superposition
                                         of the contributions from the polygon vertices; otherwise, use standard += summation.
    :return: [numEvalPos][3: x, y, z] magnetic field at evaluation locations; in T
    :rtype: arr(float)
    """
    numEvalPos = len(evalPos)
    numVertices = len(vertices)
    if numVertices < 2:
        raise ValueError("must have at least 2 vertices, not %d"%(numVertices,))

    magneticField = np.zeros((numEvalPos, 3))

    if current == 0.0:
        return magneticField

    # needs additional division by length of wire segment!
    bPrefactorL = MU_0_BY_4_PI * current

    # setup compensated summation objects
    if useCompensatedSummation:
        # need three doubles (s, cs, ccs) per eval pos --> see compsum.py
        bXSum = np.zeros((numEvalPos, 3))
        bYSum = np.zeros((numEvalPos, 3))
        bZSum = np.zeros((numEvalPos, 3))

    x_i = vertices[0, 0]
    y_i = vertices[0, 1]
    z_i = vertices[0, 2]

    for idxSource in range(numVertices-1):

        x_f = vertices[idxSource + 1, 0]
        y_f = vertices[idxSource + 1, 1]
        z_f = vertices[idxSource + 1, 2]

        # vector from start to end of i:th wire segment
        dx = x_f - x_i
        dy = y_f - y_i
        dz = z_f - z_i

        # squared length of wire segment
        l2 = dx * dx + dy * dy + dz * dz
        if l2 == 0.0:
            # skip zero-length segments: no contribution
            continue

        # length of wire segment
        l = np.sqrt(l2)

        # assemble full prefactor for B_phi
        bPrefactor = bPrefactorL / l

        # unit vector parallel to wire segment
        eX = dx / l
        eY = dy / l
        eZ = dz / l

        for idxEval in range(numEvalPos):

            # vector from start of wire segment to eval pos
            r0x = evalPos[idxEval, 0] - x_i
            r0y = evalPos[idxEval, 1] - y_i
            r0z = evalPos[idxEval, 2] - z_i

            # z position along axis of wire segment
            alignedZ = eX * r0x + eY * r0y + eZ * r0z

            # normalized z component of evaluation location in coordinate system of wire segment
            zP = alignedZ / l

            # vector perpendicular to axis of wire segment, pointing at evaluation pos
            rPerpX = r0x - alignedZ * eX
            rPerpY = r0y - alignedZ * eY
            rPerpZ = r0z - alignedZ * eZ

            # perpendicular distance squared between evalPos and axis of wire segment
            alignedRSq = rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ

            # B_phi is zero along axis of filament
            if alignedRSq > 0.0:

                # perpendicular distance between evalPos and axis of wire segment
                alignedR = np.sqrt(alignedRSq)

                # normalized rho component of evaluation location in coordinate system of wire segment
                rhoP = alignedR / l

                # compute tangential component of magnetic vector potential, including current and mu_0
                bPhi = bPrefactor * straightWireSegment_B_phi(rhoP, zP)

                # unit vector in radial direction
                eRX = rPerpX / alignedR
                eRY = rPerpY / alignedR
                eRZ = rPerpZ / alignedR

                # compute cross product between e_z and e_rho to get e_phi
                ePhiX = eY * eRZ - eZ * eRY
                ePhiY = eZ * eRX - eX * eRZ
                ePhiZ = eX * eRY - eY * eRX

                # add contribution from wire segment to result
                if useCompensatedSummation:
                    compAdd(bPhi * ePhiX, bXSum[idxEval, :])
                    compAdd(bPhi * ePhiY, bYSum[idxEval, :])
                    compAdd(bPhi * ePhiZ, bZSum[idxEval, :])
                else:
                    magneticField[idxEval, 0] += bPhi * ePhiX
                    magneticField[idxEval, 1] += bPhi * ePhiY
                    magneticField[idxEval, 2] += bPhi * ePhiZ

        # shift to next point
        x_i = x_f
        y_i = y_f
        z_i = z_f

    if useCompensatedSummation:
        # obtain compensated sums from summation objects
        for idxEval in range(numEvalPos):
            magneticField[idxEval, 0] = bXSum[idxEval, 0] + bXSum[idxEval, 1] + bXSum[idxEval, 2]
            magneticField[idxEval, 1] = bYSum[idxEval, 0] + bYSum[idxEval, 1] + bYSum[idxEval, 2]
            magneticField[idxEval, 2] = bZSum[idxEval, 0] + bZSum[idxEval, 1] + bZSum[idxEval, 2]

    return magneticField

def magneticFieldPolygonFilamentVertexSupplier(numVertices, vertexSupplier, current, evalPos, useCompensatedSummation = True):
    """Compute the magnetic field of a polygon filament at a number of evaluation locations.

    :param int numVertices: number of polygon vertices to take into account
    :param callable(int i), arr(float) vertexSupplier: should return points along polygon as [3: x, y, z] for i=0,1,...,(numVertices-1); in m
    :param float current: current along polygon; in A
    :param arr(float) evalPos: [numEvalPos][3: x, y, z] evaluation locations; in m
    :param bool useCompensatedSummation: If true, use Kahan-Babuska compensated summation to compute the superposition
                                         of the contributions from the polygon vertices; otherwise, use standard += summation.
    :return: [numEvalPos][3: x, y, z] magnetic field at evaluation locations; in T
    :rtype: arr(float)
    """
    numEvalPos = len(evalPos)
    if numVertices < 2:
        raise ValueError("must have at least 2 vertices, not %d"%(numVertices,))

    magneticField = np.zeros((numEvalPos, 3))

    if current == 0.0:
        return magneticField

    # needs additional division by length of wire segment!
    bPrefactorL = MU_0_BY_4_PI * current

    # setup compensated summation objects
    if useCompensatedSummation:
        # need three doubles (s, cs, ccs) per eval pos --> see compsum.py
        bXSum = np.zeros((numEvalPos, 3))
        bYSum = np.zeros((numEvalPos, 3))
        bZSum = np.zeros((numEvalPos, 3))

    pointData = vertexSupplier(0)
    x_i = pointData[0]
    y_i = pointData[1]
    z_i = pointData[2]

    for idxSource in range(numVertices-1):

        pointData = vertexSupplier(idxSource + 1)
        x_f = pointData[0]
        y_f = pointData[1]
        z_f = pointData[2]

        # vector from start to end of i:th wire segment
        dx = x_f - x_i
        dy = y_f - y_i
        dz = z_f - z_i

        # squared length of wire segment
        l2 = dx * dx + dy * dy + dz * dz
        if l2 == 0.0:
            # skip zero-length segments: no contribution
            continue

        # length of wire segment
        l = np.sqrt(l2)

        # assemble full prefactor for B_phi
        bPrefactor = bPrefactorL / l

        # unit vector parallel to wire segment
        eX = dx / l
        eY = dy / l
        eZ = dz / l

        for idxEval in range(numEvalPos):

            # vector from start of wire segment to eval pos
            r0x = evalPos[idxEval, 0] - x_i
            r0y = evalPos[idxEval, 1] - y_i
            r0z = evalPos[idxEval, 2] - z_i

            # z position along axis of wire segment
            alignedZ = eX * r0x + eY * r0y + eZ * r0z

            # normalized z component of evaluation location in coordinate system of wire segment
            zP = alignedZ / l

            # vector perpendicular to axis of wire segment, pointing at evaluation pos
            rPerpX = r0x - alignedZ * eX
            rPerpY = r0y - alignedZ * eY
            rPerpZ = r0z - alignedZ * eZ

            # perpendicular distance squared between evalPos and axis of wire segment
            alignedRSq = rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ

            # B_phi is zero along axis of filament
            if alignedRSq > 0.0:

                # perpendicular distance between evalPos and axis of wire segment
                alignedR = np.sqrt(alignedRSq)

                # normalized rho component of evaluation location in coordinate system of wire segment
                rhoP = alignedR / l

                # compute tangential component of magnetic vector potential, including current and mu_0
                bPhi = bPrefactor * straightWireSegment_B_phi(rhoP, zP)

                # unit vector in radial direction
                eRX = rPerpX / alignedR
                eRY = rPerpY / alignedR
                eRZ = rPerpZ / alignedR

                # compute cross product between e_z and e_rho to get e_phi
                ePhiX = eY * eRZ - eZ * eRY
                ePhiY = eZ * eRX - eX * eRZ
                ePhiZ = eX * eRY - eY * eRX

                # add contribution from wire segment to result
                if useCompensatedSummation:
                    compAdd(bPhi * ePhiX, bXSum[idxEval, :])
                    compAdd(bPhi * ePhiY, bYSum[idxEval, :])
                    compAdd(bPhi * ePhiZ, bZSum[idxEval, :])
                else:
                    magneticField[idxEval, 0] += bPhi * ePhiX
                    magneticField[idxEval, 1] += bPhi * ePhiY
                    magneticField[idxEval, 2] += bPhi * ePhiZ

        # shift to next point
        x_i = x_f
        y_i = y_f
        z_i = z_f

    if useCompensatedSummation:
        # obtain compensated sums from summation objects
        for idxEval in range(numEvalPos):
            magneticField[idxEval, 0] = bXSum[idxEval, 0] + bXSum[idxEval, 1] + bXSum[idxEval, 2]
            magneticField[idxEval, 1] = bYSum[idxEval, 0] + bYSum[idxEval, 1] + bYSum[idxEval, 2]
            magneticField[idxEval, 2] = bZSum[idxEval, 0] + bZSum[idxEval, 1] + bZSum[idxEval, 2]

    return magneticField
