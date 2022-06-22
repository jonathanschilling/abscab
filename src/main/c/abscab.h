#ifndef ABSCAB_H
#define ABSCAB_H

#include "cel.h"

/////// A_z of straight wire segment

/**
 * Compute the normalized axial component of magnetic vector potential of straight wire segment,
 * evaluated along axis of wire segment (rho = 0).
 * This is a special case for points away from the wire ("far-field") for zP < -1 or zP >= 2.
 *
 * @param zP normalized axial coordinate of evaluation location; must not be in [0, 1] (on wire segment)
 * @return normalized axial component of magnetic vector potential
 */
double sws_A_z_ax_f(double zP) {
	return atanh(1 / (fabs(zP) + fabs(1 - zP)));
}

/**
 * Compute the normalized axial component of magnetic vector potential of straight wire segment,
 * evaluated along axis of wire segment (rhoP = 0).
 * This is a special case for points close to the wire ("near-field") for -1 <= zP < 2.
 *
 * @param zP normalized axial coordinate of evaluation location; must not be in [0, 1] (on wire segment)
 * @return normalized axial component of magnetic vector potential
 */
double sws_A_z_ax_n(double zP) {
	return copysign(log(zP / (zP - 1)) / 2, zP);
}

/**
 * Compute the normalized axial component of magnetic vector potential of straight wire segment,
 * evaluated along axis of wire segment (rho = 0).
 *
 * @param zP normalized axial coordinate of evaluation location; must not be in [0, 1] (on wire segment)
 * @return normalized axial component of magnetic vector potential
 */
double sws_A_z_ax(double zP) {
	if (zP < -1 || zP >= 2) {
		return sws_A_z_ax_f(zP);
	} else {
		return sws_A_z_ax_n(zP);
	}
}

/**
 * Compute the normalized axial component of the magnetic vector potential of a straight wire segment,
 * evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
 * This is a special case for points away from the wire ("far-field") for rhoP > 1.
 *
 * @param rhoP normalized radial coordinate of evaluation location; must not be zero (on wire segment)
 * @return normalized axial component of magnetic vector potential
 */
double sws_A_z_rad_f(double rhoP) {
	return atanh(1 / (rhoP + hypot(rhoP, 1)));
}

/**
 * Compute the normalized axial component of the magnetic vector potential of a straight wire segment,
 * evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
 * This is a special case for points close to the wire ("near-field") for rhoP <= 1.
 *
 * @param rhoP normalized radial coordinate of evaluation location; must not be zero (on wire segment)
 * @return normalized axial component of magnetic vector potential
 */
double sws_A_z_rad_n(double rhoP) {
	double cat = 1 / hypot(rhoP, 1);  // cos(atan(...)  )
	double sat = sin(atan(rhoP) / 2); // sin(atan(...)/2)
	double rc = rhoP * cat;
	double num = rc + 1 + cat;
	double den = rc + 2 * sat * sat;
	return log(num / den) / 2;
}

/**
 * Compute the normalized axial component of the magnetic vector potential of a straight wire segment,
 * evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
 *
 * @param rhoP normalized radial coordinate of evaluation location; must not be zero (on wire segment)
 * @return normalized axial component of magnetic vector potential
 */
double sws_A_z_rad(double rhoP) {
	if (rhoP > 1) {
		return sws_A_z_rad_f(rhoP);
	} else {
		return sws_A_z_rad_n(rhoP);
	}
}

/**
 * Compute the normalized axial component of the magnetic vector potential of a straight wire segment,
 * evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
 * This formulation is useful for points away from the wire ("far-field")
 * at rhoP >= 1 or zP <= -1 or zP > 2.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized axial component of magnetic vector potential
 */
double sws_A_z_f(double rhoP, double zP) {
	double r_i = hypot(rhoP, zP);
	double r_f = hypot(rhoP, 1 - zP);
	return atanh(1 / (r_i + r_f));
}

/**
 * Compute the normalized axial component of the magnetic vector potential of a straight wire segment,
 * evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
 * This formulation is useful for points close to the wire ("near-field")
 * at rhoP < 1 and -1 < zP <= 2.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized axial component of magnetic vector potential
 */
double sws_A_z_n(double rhoP, double zP) {
	double omz = 1 - zP;

	double r_i = hypot(rhoP, zP);
	double r_f = hypot(rhoP, omz);

	double alpha = atan2(rhoP, zP);
	double sinAlphaHalf = sin(alpha / 2);

	double beta = atan2(rhoP, omz);
	double sinBetaHalf = sin(beta / 2);

	double Ri_zP    = 2 * r_i * sinAlphaHalf * sinAlphaHalf; // r_i - z'
	double Rf_p_zM1 = 2 * r_f * sinBetaHalf  * sinBetaHalf;  // r_f - (1 - z')

	double n = Ri_zP + Rf_p_zM1;

	return (log(2 + n) - log(n)) / 2;
}

/////// B_phi of straight wire segment

/**
 * Compute the normalized tangential component of the magnetic field of a straight wire segment,
 * evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @return normalized tangential component of magnetic field
 */
double sws_B_phi_rad(double rhoP) {
	return 1 / (rhoP * hypot(rhoP, 1));
}

/**
 * Compute the normalized tangential component of the magnetic field of a straight wire segment.
 * This formulation is useful for points away from the wire ("far-field")
 * at rhoP >= 1 or zP <= 0 or zP >= 1 or rhoP/(1-zP) >= 1 or rhoP/zP >= 1.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized tangential component of magnetic field
 */
double sws_B_phi_f(double rhoP, double zP) {
	double omz = 1 - zP;

	double r_i = hypot(rhoP, zP);
	double r_f = hypot(rhoP, omz);

	double num = rhoP * (1/r_i + 1/r_f);
	double den = rhoP * rhoP - zP * omz + r_i * r_f;

	return num / den;
}

/**
 * Compute the normalized tangential component of the magnetic field of a straight wire segment.
 * This formulation is useful for points close to the wire ("near-field")
 * at rhoP < 1 and 0 < zP < 1 and rhoP/(1-zP) < 1 and rhoP/zP < 1.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized tangential component of magnetic field
 */
double sws_B_phi_n(double rhoP, double zP) {
	double omz = 1 - zP;

	double r_i = hypot(rhoP, zP);
	double r_f = hypot(rhoP, omz);

	double num = rhoP * (1/r_i + 1/r_f);

	double alpha = atan2(rhoP, zP);
	double sinAlphaHalf = sin(alpha / 2);

	double beta = atan2(rhoP, omz);
	double sinBetaHalf = sin(beta / 2);

	// r_f * sin^2(beta/2) + (1 - z') * sin^2(alpha/2)
	double rfb_omza = r_f * sinBetaHalf * sinBetaHalf + omz * sinAlphaHalf * sinAlphaHalf;

	//     r_i * r_f - z' * (1 - z')
	// =   r_i * r_f - r_i * (1 - z') + r_i * (1 - z') - z' * (1 - z')
	// =   r_i * r_f - r_i * r_f * cos(beta)
	//   + r_i * (1 - z') + (1 - z') * r_i * cos(alpha)
	// =   r_i *    r_f   * (1 - cos(beta))
	//   + r_i * (1 - z') * (1 - cos(alpha))
	// = 2 * r_i * [ r_f * sin^2(beta/2) + (1 - z') * sin^2(alpha/2) ]
	double den = rhoP * rhoP + 2 * r_i * rfb_omza;

	return num / den;
}

///// A_phi of circular wire loop

/**
 * Compute the normalized tangential component of the magnetic vector potential of a circular wire loop.
 * This formulation is useful for points away from the wire ("far-field")
 * at rhoP < 1/2 or rhoP > 2 or |zP| >= 1.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized tangential component of magnetic vector potential
 */
double cwl_A_phi_f(double rhoP, double zP) {
	double sqrt_kCSqNum = hypot(zP, 1 - rhoP);
	double sqrt_kCSqDen = hypot(zP, 1 + rhoP);

	double kC = sqrt_kCSqNum / sqrt_kCSqDen;
	double kSq = 4 * rhoP / (sqrt_kCSqDen * sqrt_kCSqDen);

	double kCp1 = 1 + kC;
	double arg1 = 2 * sqrt(kC) / kCp1;
	double arg2 = 2 / (kCp1 * kCp1 * kCp1);
	double C = cel(arg1, 1, 0, arg2);

	return kSq/sqrt_kCSqDen * C;
}

/**
 * Compute the normalized tangential component of the magnetic vector potential of a circular wire loop.
 * This formulation is useful for points close to the wire ("near-field")
 * at 1/2 <= rhoP <= 2 and |zP| < 1.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized tangential component of magnetic vector potential
 */
double cwl_A_phi_n(double rhoP, double zP) {
	double rhoP_m_1 = rhoP - 1;

	double n = zP / rhoP_m_1;
	double m = 1 + 2 / rhoP_m_1;

	double num = n * n + 1;
	double den = n * n + m * m;

	double kCSq = num / den;

	double prefac = 1 / (fabs(rhoP - 1) * sqrt(den));
	double celPart = cel(sqrt(kCSq), 1, -1, 1);
	return prefac * celPart;
}

/**
 * Compute the normalized tangential component of the magnetic vector potential of a circular wire loop.
 * This formulation is useful for points along rhoP=1 with |zP| < 1.
 *
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized tangential component of magnetic vector potential
 */
double cwl_A_phi_v(double zP) {
	double absZp = fabs(zP);

	// 1/k_c
	double kCInv = sqrt(4 + zP * zP) / absZp;

	return cel(kCInv, 1, 1, -1) / absZp;
}

//////// B_rho of circular wire loop

/**
 * Compute the normalized radial component of the magnetic field of a circular wire loop.
 * This formulation is useful for points away from the wire ("far-field")
 * at rhoP < 1/2 or rhoP > 2 or |zP| >= 1.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized radial component of magnetic field
 */
double cwl_B_rho_f(double rhoP, double zP) {
	double sqrt_kCSqNum = hypot(zP, 1 - rhoP);
	double sqrt_kCSqDen = hypot(zP, 1 + rhoP);

	double kCSqNum = sqrt_kCSqNum * sqrt_kCSqNum;
	double kCSqDen = sqrt_kCSqDen * sqrt_kCSqDen;

	double kCSq = kCSqNum / kCSqDen;
	double kC = sqrt(kCSq);

	double D = cel(kC, 1, 0, 1);

	double kCp1 = 1 + kC;
	double arg1 = 2 * sqrt(kC) / kCp1;
	double arg2 = 2 / (kCp1 * kCp1 * kCp1);
	double C = cel(arg1, 1, 0, arg2);

	double prefac = 4 * rhoP / (kCSqDen * sqrt_kCSqDen * kCSqNum);

	return prefac * zP * (D - C);
}

/**
 * Compute the normalized radial component of the magnetic field of a circular wire loop.
 * This formulation is useful for points close to the wire ("near-field")
 * at 1/2 <= rhoP <= 2 and |zP| < 1.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized radial component of magnetic field
 */
double cwl_B_rho_n(double rhoP, double zP) {
	double rhoP_m_1 = rhoP - 1;
	double rd2 = rhoP_m_1 * rhoP_m_1;

	double n = zP / rhoP_m_1;
	double m = 1 + 2 / rhoP_m_1;

	double sqrt_kCSqNum = hypot(n, 1);
	double sqrt_kCSqDen = hypot(n, m);

	double kCSqNum = sqrt_kCSqNum * sqrt_kCSqNum;
	double kCSqDen = sqrt_kCSqDen * sqrt_kCSqDen;

	double kC = sqrt_kCSqNum / sqrt_kCSqDen;

	double D = cel(kC, 1, 0, 1);

	double kCp1 = 1 + kC;
	double arg1 = 2 * sqrt(kC) / kCp1;
	double arg2 = 2 / (kCp1 * kCp1 * kCp1);
	double C = arg2 * cel(arg1, 1, 0, 1);

	// z' / |rho' - 1|^5
	double zP_rd5 = zP / (fabs(rhoP_m_1) * rd2 * rd2);

	double prefac = 4 * rhoP / (kCSqDen * sqrt_kCSqDen * kCSqNum);

	return prefac * zP_rd5 * (D - C);
}

/**
 * Compute the normalized radial component of the magnetic field of a circular wire loop.
 * This formulation is useful for points along rhoP=1 with |zP| < 1.
 *
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized radial component of magnetic field
 */
double cwl_B_rho_v(double zP) {
	double zPSq = zP * zP;

	double kCSq = 1 / (1 + 4 / zPSq);
	double kC = sqrt(kCSq);

	double K = cel(kC, 1, 1, 1);
	double E = cel(kC, 1, 1, kCSq);

	return copysign(kC / 2 * ((2 / zPSq + 1) * E - K), zP);
}

////// B_z of circular wire loop

/**
 * Compute the normalized vertical component of the magnetic field of a circular wire loop.
 * This formulation is useful for certain points away from the wire ("far-field")
 * at rhoP < 1/2 or (rhoP <= 2 and |zP| >= 1).
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized vertical component of magnetic field
 */
double cwl_B_z_f1(double rhoP, double zP) {
	double sqrt_kCSqNum = hypot(zP, 1 - rhoP);
	double sqrt_kCSqDen = hypot(zP, 1 + rhoP);

	double kC = sqrt_kCSqNum / sqrt_kCSqDen;

	double K = cel(kC, 1, 1, 1);
	double E = cel(kC, 1, 1, kC * kC);
	double D = cel(kC, 1, 0, 1);

	double prefac = 1 / (sqrt_kCSqDen * sqrt_kCSqNum * sqrt_kCSqNum);
	double comb = (E - 2 * K + 2 * D);

	return prefac * (E + rhoP * comb);
}

/**
 * Compute the normalized vertical component of the magnetic field of a circular wire loop.
 * This formulation is useful for certain other points away from the wire ("far-field")
 * at rhoP > 2.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized vertical component of magnetic field
 */
double cwl_B_z_f2(double rhoP, double zP) {
	double sqrt_kCSqNum = hypot(zP, 1 - rhoP);
	double sqrt_kCSqDen = hypot(zP, 1 + rhoP);

	double kC = sqrt_kCSqNum / sqrt_kCSqDen;
	double kCSq = kC * kC;

	double zPSqP1 = zP * zP + 1;
	double rhoPSq = rhoP * rhoP;
	double t1 = zPSqP1 / rhoPSq + 1;
	double t2 = 2 / rhoP;

	// a is sqrt_kCSqDen normalized to rho'^2
	// b is sqrt_kCSqNum normalized to rho'^2
	// a == (z'^2 + (1 + rho')^2) / rho'^2 = (z'^2 + 1)/rho'^2 + 1  +  2/rho'
	// b == (z'^2 + (1 - rho')^2) / rho'^2 = (z'^2 + 1)/rho'^2 + 1  -  2/rho'
	double a = t1 + t2;
	double b = t1 - t2;

	// 1/prefac = sqrt( z'^2 + (1 + rho')^2)           * (z'^2 + (1 - rho')^2)
	//          = sqrt((z'^2 + (1 + rho')^2) / rho'^2) * (z'^2 + (1 - rho')^2) / rho'^2 * rho'^3
	//          = sqrt(a)                              * b                              * rho'^3
	double prefac = 1 / (sqrt(a) * b * rhoPSq * rhoP);

	double cdScale = 1 + (2 + zPSqP1 / rhoP) / rhoP;

	double E = cel(kC, 1, 1, kCSq);
	double D = cel(kC, 1, 0, 1);

	double kCP1 = 1 + kC;
	double arg1 = 2 * sqrt(kC) / kCP1;
	double arg2 = 2 / (kCP1 * kCP1 * kCP1);
	double C = arg2 * cel(arg1, 1, 0, 1);

	// use C - D for (2 * D - E)/kSq
	return prefac * (E + 4 * (C - D) / cdScale);
}

/**
 * Compute the normalized vertical component of the magnetic field of a circular wire loop.
 * This formulation is useful for points close to the wire ("near-field")
 * at 1/2 <= rhoP <= 2, but not rhoP=1, and |zP| <= 1.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized vertical component of magnetic field
 */
double cwl_B_z_n(double rhoP, double zP) {
	double rp1 = rhoP - 1;

	double n = zP / rp1;
	double m = 1 + 2 / rp1;

	double sqrt_kCSqNum = hypot(n, 1);
	double sqrt_kCSqDen = hypot(n, m);

	double kCSqDen = sqrt_kCSqDen * sqrt_kCSqDen;

	double kC = sqrt_kCSqNum / sqrt_kCSqDen;

	double prefac = 1 / (fabs(rp1) * rp1 * rp1 * kCSqDen * sqrt_kCSqDen);

	return prefac * cel(kC, kC * kC, 1 + rhoP, 1 - rhoP);
}

/**
 * Compute the normalized vertical component of the magnetic field of a circular wire loop.
 * This formulation is useful for points along rhoP=1 with |zP| <= 1.
 *
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized vertical component of magnetic field
 */
double cwl_B_z_v(double zP) {
	double kCSq = zP * zP / (4 + zP * zP);
	double kC = sqrt(kCSq);

	double f = zP * zP + 4;
	double prefac = 1 / (f * sqrt(f));

	return prefac * cel(kC, kCSq, 2, 0);
}

// --------------------------------------------------

/**
 * Compute the normalized axial component of the magnetic vector potential of a straight wire segment.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized axial component of magnetic vector potential
 */
double straightWireSegment_A_z(double rhoP, double zP) {
	if (rhoP == 0.0) {
		return sws_A_z_ax(zP);
	} else if (zP == 0.0 || zP == 1.0) {
		return sws_A_z_rad(rhoP);
	} else if (rhoP >= 1.0 || zP <= -1.0 || zP > 2.0) {
		return sws_A_z_f(rhoP, zP);
	} else {
		return sws_A_z_n(rhoP, zP);
	}
}

/**
 * Compute the normalized tangential component of the magnetic field of a straight wire segment.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized tangential component of magnetic field
 */
double straightWireSegment_B_phi(double rhoP, double zP) {
	if (rhoP == 0.0) {
		return 0.0;
	} else if (zP == 0.0 || zP == 1.0) {
		return sws_B_phi_rad(rhoP);
	} else if (rhoP >= 1.0 || zP <= 0.0 || zP >= 1.0 || rhoP / (1 - zP) >= 1.0 || rhoP / zP >= 1.0) {
		return sws_B_phi_f(rhoP, zP);
	} else {
		return sws_B_phi_n(rhoP, zP);
	}
}

/**
 * Geometric part of magnetic vector potential computation for circular wire
 * loop at rho'=1, z'=0 (normalized coordinates). This routine selects special
 * case routines to get the most accurate formulation for given evaluation
 * coordinates.
 *
 * @param rhoP normalized radial evaluation position
 * @param zP   normalized vertical evaluation position
 * @return A_phi: toroidal component of magnetic vector potential: geometric
 *         part (no mu0*I/pi factor included)
 */
double circularWireLoop_A_phi(double rhoP, double zP) {
	if (rhoP == 0.0) {
		return 0.0;
	} else if (rhoP < 0.5 || rhoP > 2.0 || fabs(zP) >= 1.0) {
		return cwl_A_phi_f(rhoP, zP);
	} else if (rhoP != 1.0) {
		return cwl_A_phi_n(rhoP, zP);
	} else {
		return cwl_A_phi_v(zP);
	}
}

/**
 * Geometric part of radial magnetic field computation for circular wire loop at
 * rho'=1, z'=0 (normalized coordinates). This routine selects special case
 * routines to get the most accurate formulation for given evaluation
 * coordinates.
 *
 * @param rhoP normalized radial evaluation position
 * @param zP   normalized vertical evaluation position
 * @return B_rho: radial component of magnetic field: geometric part (no
 *         mu0*I/(pi*a) factor included)
 */
double circularWireLoop_B_rho(double rhoP, double zP) {
	if (rhoP == 0.0 || zP == 0.0) {
		return 0.0;
	} else if (rhoP < 0.5 || rhoP > 2.0 || fabs(zP) >= 1.0) {
		return cwl_B_rho_f(rhoP, zP);
	} else if (rhoP != 1.0) {
		return cwl_B_rho_n(rhoP, zP);
	} else {
		return cwl_B_rho_v(zP);
	}
}

/**
 * Geometric part of vertical magnetic field computation for circular wire loop
 * at rho'=1, z'=0 (normalized coordinates). This routine selects special case
 * routines to get the most accurate formulation for given evaluation
 * coordinates.
 *
 * @param rhoP normalized radial evaluation position
 * @param zP   normalized vertical evaluation position
 * @return B_z: vertical component of magnetic field: geometric part (no
 *         mu0*I/(pi*a) factor included)
 */
double circularWireLoop_B_z(double rhoP, double zP) {
	if (rhoP < 0.5 || (rhoP <= 2 && fabs(zP) > 1)) {
		return cwl_B_z_f1(rhoP, zP);
	} else if (rhoP > 2) {
		return cwl_B_z_f2(rhoP, zP);
	} else if (rhoP != 1.0) {
		return cwl_B_z_n(rhoP, zP);
	} else {
		return cwl_B_z_v(zP);
	}
}

// --------------------------------------------------







#endif // ABSCAB_H
