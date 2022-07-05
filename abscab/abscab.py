import numpy as np
from bulirsch_cel import cel

"""vacuum magnetic permeability in Vs/Am (CODATA-2018)"""
MU_0 = 1.25663706212e-6

"""vacuum magnetic permeability, divided by pi"""
MU_0_BY_PI = MU_0 / np.pi

"""vacuum magnetic permeability, divided by 2 pi"""
MU_0_BY_2_PI = MU_0 / (2.0 * np.pi)

"""vacuum magnetic permeability, divided by 4 pi"""
MU_0_BY_4_PI = MU_0 / (4.0 * np.pi)

############## A_z of straight wire segment

def sws_A_z_ax_f(zP):
    """Normalized A_z of Straight Wire Segment, along rhoP=0, far-field.
    
    Compute the normalized axial component of magnetic vector potential of straight wire segment,
    evaluated along axis of wire segment (rho = 0).
    This is a special case for points away from the wire ("far-field") for zP < -1 or zP >= 2.
    
    :param float zP: normalized vertical coordinate of evaluation location
    :return: normalized axial component of magnetic vector potential
    :rtype: float
    """
    return np.arctanh(1 / (np.abs(zP) + np.abs(1 - zP)))

def sws_A_z_ax_n(zP):
    """Normalized A_z of Straight Wire Segment, along rhoP=0, near-field.
    
    Compute the normalized axial component of magnetic vector potential of straight wire segment,
    evaluated along axis of wire segment (rhoP = 0).
    This is a special case for points close to the wire ("near-field") for -1 <= zP < 2.
    
    :param float zP: normalized axial coordinate of evaluation location; must not be in [0, 1] (on wire segment)
    :return: normalized axial component of magnetic vector potential
    :rtype: float
    """
    
    # Two negative signs must be able to cancel each other here!
    return np.copysign(1.0, zP) * np.log(zP / (zP - 1)) / 2

def sws_A_z_ax(zP):
    """Normalized A_z of Straight Wire Segment, along rhoP=0.
    
    Compute the normalized axial component of magnetic vector potential of straight wire segment,
    evaluated along axis of wire segment (rho = 0).
    
    :param float zP: normalized axial coordinate of evaluation location; must not be in [0, 1] (on wire segment)
    :return: normalized axial component of magnetic vector potential
    :rtype: float
    """
    if zP < -1 or zP >= 2:
        return sws_A_z_ax_f(zP)
    else:
        return sws_A_z_ax_n(zP)

def sws_A_z_rad_f(rhoP):
    """Normalized A_z of Straight Wire Segment, along zP=0 or zP=1, far-field.
    
    Compute the normalized axial component of the magnetic vector potential of a straight wire segment,
    evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
    This is a special case for points away from the wire ("far-field") for rhoP > 1.
  
    :param float rhoP: normalized radial coordinate of evaluation location; must not be zero (on wire segment)
    :return: normalized axial component of magnetic vector potential
    :rtype: float
    """
    return np.arctanh(1 / (rhoP + np.hypot(rhoP, 1)))

def sws_A_z_rad_n(rhoP):
    """Normalized A_z of Straight Wire Segment, along zP=0 or zP=1, near-field.
    
    Compute the normalized axial component of the magnetic vector potential of a straight wire segment,
    evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
    This is a special case for points close to the wire ("near-field") for rhoP <= 1.
    
    :param float rhoP: normalized radial coordinate of evaluation location; must not be zero (on wire segment)
    :return: normalized axial component of magnetic vector potential
    :rtype: float
    """
    cat = 1 / np.hypot(rhoP, 1)  # cos(atan(...)  )
    sat = np.sin(np.arctan(rhoP) / 2) # sin(atan(...)/2)
    rc = rhoP * cat
    num = rc + 1 + cat
    den = rc + 2 * sat * sat
    return np.log(num / den) / 2

def sws_A_z_rad(rhoP):
    """Normalized A_z of Straight Wire Segment, along zP=0 or zP=1.
    
    Compute the normalized axial component of the magnetic vector potential of a straight wire segment,
    evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
    
    :param float rhoP: normalized radial coordinate of evaluation location; must not be zero (on wire segment)
    :return: normalized axial component of magnetic vector potential
    :rtype: float
    """
    if rhoP > 1:
        return sws_A_z_rad_f(rhoP)
    else:
        return sws_A_z_rad_n(rhoP)

def sws_A_z_f(rhoP, zP):
    """Normalized A_z of Straight Wire Segment, far-field.
    
    Compute the normalized axial component of the magnetic vector potential of a straight wire segment.
    This formulation is useful for points away from the wire ("far-field")
    at rhoP >= 1 or zP <= -1 or zP > 2.
 
    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized axial component of magnetic vector potential
    :rtype: float
    """
    r_i = np.hypot(rhoP, zP)
    r_f = np.hypot(rhoP, 1 - zP)
    return np.arctanh(1 / (r_i + r_f))

def sws_A_z_n(rhoP, zP):
    """Normalized A_z of Straight Wire Segment, near-field.
    
    Compute the normalized axial component of the magnetic vector potential of a straight wire segment.
    This formulation is useful for points close to the wire ("near-field")
    at rhoP < 1 and -1 < zP <= 2.
 
    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized axial component of magnetic vector potential
    :rtype: float
    """
    omz = 1 - zP

    r_i = np.hypot(rhoP, zP)
    r_f = np.hypot(rhoP, omz)

    alpha = np.atan2(rhoP, zP)
    sinAlphaHalf = np.sin(alpha / 2)

    beta = np.atan2(rhoP, omz)
    sinBetaHalf = np.sin(beta / 2)

    Ri_zP    = 2 * r_i * sinAlphaHalf * sinAlphaHalf # r_i - z'
    Rf_p_zM1 = 2 * r_f * sinBetaHalf  * sinBetaHalf  # r_f - (1 - z')

    n = Ri_zP + Rf_p_zM1

    return (np.log(2 + n) - np.log(n)) / 2

############## B_phi of straight wire segment

def sws_B_phi_rad(rhoP):
    """Normalized B_phi of Straight Wire Segment, along zP=0 or zP=1.
    
    Compute the normalized tangential component of the magnetic field of a straight wire segment,
    evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
    
    :param float rhoP: normalized radial coordinate of evaluation location
    :return: normalized tangential component of magnetic field
    :rtype: float
    """
    return 1 / (rhoP * np.hypot(rhoP, 1))

def sws_B_phi_f(rhoP, zP):
    """Normalized B_phi of Straight Wire Segment, far-field.
    
    Compute the normalized tangential component of the magnetic field of a straight wire segment.
    This formulation is useful for points away from the wire ("far-field")
    at rhoP >= 1 or zP <= 0 or zP >= 1 or rhoP/(1-zP) >= 1 or rhoP/zP >= 1.
 
    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized tangential component of magnetic field
    :rtype: float
    """
    omz = 1 - zP

    r_i = np.hypot(rhoP, zP)
    r_f = np.hypot(rhoP, omz)

    num = rhoP * (1/r_i + 1/r_f)
    den = rhoP * rhoP - zP * omz + r_i * r_f

    return num / den

def sws_B_phi_n(rhoP, zP):
    """Normalized B_phi of Straight Wire Segment, near-field.
    
    Compute the normalized tangential component of the magnetic field of a straight wire segment.
    This formulation is useful for points close to the wire ("near-field")
    at rhoP < 1 and 0 < zP < 1 and rhoP/(1-zP) < 1 and rhoP/zP < 1.
 
    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized tangential component of magnetic field
    :rtype: float
    """
    omz = 1 - zP

    r_i = np.hypot(rhoP, zP)
    r_f = np.hypot(rhoP, omz)

    num = rhoP * (1/r_i + 1/r_f)

    alpha = np.atan2(rhoP, zP)
    sinAlphaHalf = np.sin(alpha / 2)

    beta = np.atan2(rhoP, omz)
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

def cwl_A_phi_f(rhoP, zP):
    """Normalized A_phi of Circular Wire Loop, far-field.
    
    Compute the normalized tangential component of the magnetic vector potential of a circular wire loop.
    This formulation is useful for points away from the wire ("far-field")
    at rhoP < 1/2 or rhoP > 2 or |zP| >= 1.
 
    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized tangential component of magnetic vector potential
    :rtype: float
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

def cwl_A_phi_n(rhoP, zP):
    """Normalized A_phi of Circular Wire Loop, near-field.
    
    Compute the normalized tangential component of the magnetic vector potential of a circular wire loop.
    This formulation is useful for points close to the wire ("near-field")
    at 1/2 <= rhoP <= 2 and |zP| < 1.
 
    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized tangential component of magnetic vector potential
    :rtype: float
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

def cwl_A_phi_v(zP):
    """Normalized A_phi of Circular Wire Loop, along rhoP=1, near-field.
    
    Compute the normalized tangential component of the magnetic vector potential of a circular wire loop.
    This formulation is useful for points along rhoP=1 with |zP| < 1.
 
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized tangential component of magnetic vector potential
    :rtype: float
    """
    absZp = np.abs(zP)

    # 1/k_c
    kCInv = np.sqrt(4 + zP * zP) / absZp

    return cel(kCInv, 1, 1, -1) / absZp

############## B_rho of circular wire loop

def cwl_B_rho_f(rhoP, zP):
    """Normalized B_rho of Circular Wire Loop, far-field.
    
    Compute the normalized radial component of the magnetic field of a circular wire loop.
    This formulation is useful for points away from the wire ("far-field")
    at rhoP < 1/2 or rhoP > 2 or |zP| >= 1.
 
    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized radial component of magnetic field
    :rtype: float
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

def cwl_B_rho_n(rhoP, zP):
    """Normalized B_rho of Circular Wire Loop, near-field.
    
    Compute the normalized radial component of the magnetic field of a circular wire loop.
    This formulation is useful for points close to the wire ("near-field")
    at 1/2 <= rhoP <= 2 and |zP| < 1.
 
    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized radial component of magnetic field
    :rtype: float
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

def cwl_B_rho_v(zP):
    """Normalized B_rho of Circular Wire Loop, along rhoP=1, near-field.
    
    Compute the normalized tangential component of the magnetic vector potential of a circular wire loop.
    This formulation is useful for points along rhoP=1 with |zP| < 1.
 
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized radial component of magnetic field
    :rtype: float
    """
    zPSq = zP * zP

    kCSq = 1 / (1 + 4 / zPSq)
    kC = np.sqrt(kCSq)

    K = cel(kC, 1, 1, 1)
    E = cel(kC, 1, 1, kCSq)

    return np.copysign(kC / 2 * ((2 / zPSq + 1) * E - K), zP)

############## B_z of circular wire loop

def cwl_B_z_f1(rhoP, zP):
    """Normalized B_z of Circular Wire Loop, far-field (1).
    
    Compute the normalized vertical component of the magnetic field of a circular wire loop.
    This formulation is useful for certain points away from the wire ("far-field")
    at rhoP < 1/2 or (rhoP <= 2 and |zP| >= 1).
 
    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized vertical component of magnetic field
    :rtype: float
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

def cwl_B_z_f2(rhoP, zP):
    """Normalized B_z of Circular Wire Loop, far-field (2).
    
    Compute the normalized vertical component of the magnetic field of a circular wire loop.
    This formulation is useful for certain other points away from the wire ("far-field")
    at rhoP > 2.
 
    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized vertical component of magnetic field
    :rtype: float
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

def cwl_B_z_n(rhoP, zP):
    """Normalized B_z of Circular Wire Loop, near-field.
    
    Compute the normalized vertical component of the magnetic field of a circular wire loop.
    This formulation is useful for points close to the wire ("near-field")
    at 1/2 <= rhoP <= 2, but not rhoP=1, and |zP| <= 1.
 
    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized vertical component of magnetic field
    :rtype: float
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

def cwl_B_z_v(zP):
    """Normalized B_z of Circular Wire Loop, along rhoP=1, near-field.
    
    Compute the normalized tangential component of the magnetic vector potential of a circular wire loop.
    This formulation is useful for points along rhoP=1 with |zP| < 1.
 
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized radial component of magnetic field
    :rtype: float
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
    """
    if rhoP == 0.0:
        return sws_A_z_ax(zP)
    elif zP == 0.0 or zP == 1.0:
        return sws_A_z_rad(rhoP)
    elif rhoP >= 1.0 or zP <= -1.0 or zP > 2.0:
        return sws_A_z_f(rhoP, zP)
    else:
        return sws_A_z_n(rhoP, zP)

def straightWireSegment_B_phi(rhoP, zP):
    """Compute the normalized tangential component of the magnetic field of a straight wire segment.
    
    This method selects the proper special case method to use
    for computing the normalized tangential component of the magnetic field
    of a straight wire segment for given evaluation location (rhoP, zP).
    
    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized tangential component of magnetic field
    :rtype: float
    """
    if rhoP == 0.0:
        return 0.0
    elif zP == 0.0 or zP == 1.0:
        return sws_B_phi_rad(rhoP)
    elif rhoP >= 1.0 or zP <= 0.0 or zP >= 1.0 or rhoP / (1 - zP) >= 1.0 or rhoP / zP >= 1.0:
        return sws_B_phi_f(rhoP, zP)
    else:
        return sws_B_phi_n(rhoP, zP)

def circularWireLoop_A_phi(rhoP, zP):
    """Compute the normalized tangential component of the magnetic vector potential of a circular wire loop.
    
    This method selects the proper special case method to use
    for computing the normalized tangential component of the magnetic vector potential
    of a circular wire loop for given evaluation location (rhoP, zP).
    
    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized tangential component of magnetic vector potential
    :rtype: float
    """
    if rhoP == 0.0:
        return 0.0
    elif rhoP < 0.5 or rhoP > 2.0 or np.abs(zP) >= 1.0:
        return cwl_A_phi_f(rhoP, zP)
    elif rhoP != 1.0:
        return cwl_A_phi_n(rhoP, zP)
    else:
        return cwl_A_phi_v(zP)

def circularWireLoop_B_rho(rhoP, zP):
    """Compute the normalized radial component of the magnetic field of a circular wire loop.
    
    This method selects the proper special case method to use
    for computing the normalized radial component of the magnetic field
    of a circular wire loop for given evaluation location (rhoP, zP).
    
    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized radial component of magnetic field
    :rtype: float
    """
    if rhoP == 0.0 or zP == 0.0:
        return 0.0
    elif rhoP < 0.5 or rhoP > 2.0 or np.abs(zP) >= 1.0:
        return cwl_B_rho_f(rhoP, zP)
    elif rhoP != 1.0:
        return cwl_B_rho_n(rhoP, zP)
    else:
        return cwl_B_rho_v(zP)

def circularWireLoop_B_z(rhoP, zP):
    """Compute the normalized vertical component of the magnetic field of a circular wire loop.
    
    This method selects the proper special case method to use
    for computing the normalized vertical component of the magnetic field
    of a circular wire loop for given evaluation location (rhoP, zP).
    
    :param float rhoP: normalized radial coordinate of evaluation location
    :param float zP: normalized axial coordinate of evaluation location
    :return: normalized radial component of magnetic field
    :rtype: float
    """
    if rhoP < 0.5 or (rhoP <= 2 and np.abs(zP) > 1):
        return cwl_B_z_f1(rhoP, zP)
    elif rhoP > 2:
        return cwl_B_z_f2(rhoP, zP)
    elif rhoP != 1.0:
        return cwl_B_z_n(rhoP, zP)
    else:
        return cwl_B_z_v(zP)

# --------------------------------------------------


