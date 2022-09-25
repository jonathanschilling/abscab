#ifndef ABSCAB_H
#define ABSCAB_H

// for memset()
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#else // _OPENMP
#define omp_get_max_threads() 1
#endif // _OPENMP

#include "cel.h"
#include "compsum.h"

#define min(x,y) ((x) < (y) ? (x) : (y))

/** vacuum magnetic permeability in Vs/Am (CODATA-2018) */
const double MU_0 = 1.25663706212e-6;

/** vacuum magnetic permeability, divided by pi */
const double MU_0_BY_PI = MU_0 / M_PI;

/** vacuum magnetic permeability, divided by 2 pi */
const double MU_0_BY_2_PI = MU_0 / (2.0 * M_PI);

/** vacuum magnetic permeability, divided by 4 pi */
const double MU_0_BY_4_PI = MU_0 / (4.0 * M_PI);

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
	// Two negative signs must be able to cancel each other here!
	return copysign(1.0, zP) * log(zP / (zP - 1)) / 2;
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
 * Compute the normalized axial component of the magnetic vector potential of a straight wire segment.
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
 * Compute the normalized axial component of the magnetic vector potential of a straight wire segment.
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

	double Ri_zP    = r_i * sinAlphaHalf * sinAlphaHalf; // 0.5 * (r_i - z')
	double Rf_p_zM1 = r_f * sinBetaHalf  * sinBetaHalf;  // 0.5 * (r_f - (1 - z'))

	double n = Ri_zP + Rf_p_zM1;

	return (log(1 + n) - log(n)) / 2;
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
		if (zP < 0 || zP > 1) {
			return sws_A_z_ax(zP);
		} else {
			fprintf(stderr, "evaluation locations on the wire segment (rho'=%g z'=%g) are not allowed\n", rhoP, zP);
		}
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
		if (zP < 0 || zP > 1) {
			return 0.0;
		} else {
			fprintf(stderr, "evaluation locations on the wire segment (rho'=%g z'=%g) are not allowed\n", rhoP, zP);
		}
	} else if (zP == 0.0 || zP == 1.0) {
		return sws_B_phi_rad(rhoP);
	} else if (rhoP >= zP || rhoP >= 1 - zP || zP < 0.0 || zP > 1.0) {
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
		if (zP != 0) {
			return cwl_A_phi_v(zP);
		} else {
			fprintf(stderr, "evaluation at location of wire loop (rho' = 1, z' = 0) is not defined\n");
		}
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
		if (rhoP != 1.0) {
			return 0.0;
		} else {
			fprintf(stderr, "evaluation at location of wire loop (rho' = 1, z' = 0) is not defined\n");
		}
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
		if (zP != 0) {
			return cwl_B_z_v(zP);
		} else {
			fprintf(stderr, "evaluation at location of wire loop (rho' = 1, z' = 0) is not defined\n");
		}
	}
}

// --------------------------------------------------

/**
 * Compute the magnetic vector potential of a circular wire loop.
 *
 * @param center  [3: x, y, z] origin of loop (in meters)
 * @param normal  [3: x, y, z] normal vector of loop (in meters); will be
 *                normalized internally
 * @param radius  radius of the wire loop (in meters)
 * @param current loop current (in A)
 * @param evalPos [3: x, y, z][nEvalPos] evaluation locations (in meters)
 * @param vectorPotential [3: A_x, A_y, A_z][nEvalPos] Cartesian components of the magnetic
 *         vector potential evaluated at the given locations (in Tm); has to be allocated on entry
 */
void vectorPotentialCircularFilament(double *center, double *normal, double radius,
		double current, int nEvalPos, double *evalPos, double *vectorPotential) {

	if (!isfinite(radius) || radius <= 0.0) {
		printf("radius must be finite and positive, but is %g\n", radius);
		return;
	}

	double aPrefactor = MU_0_BY_PI * current;

	// squared length of normal vector
	double nLen2 = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];

	if (nLen2 == 0.0) {
		printf("length of normal vector must not be zero");
		return;
	}

	// length of normal vector
	double nLen = sqrt(nLen2);

	// unit normal vector of wire loop
	double eX = normal[0] / nLen;
	double eY = normal[1] / nLen;
	double eZ = normal[2] / nLen;

	for (int idxEval = 0; idxEval < nEvalPos; ++idxEval) {

		// vector from center of wire loop to eval pos
		double r0x = evalPos[3 * idxEval + 0] - center[0];
		double r0y = evalPos[3 * idxEval + 1] - center[1];
		double r0z = evalPos[3 * idxEval + 2] - center[2];

		// z position along normal of wire loop
		double alignedZ = eX * r0x + eY * r0y + eZ * r0z;

		// normalized z component of evaluation location in coordinate system of wire loop
		double zP = alignedZ / radius;

		// r0 projected onto axis of wire loop
		double rParallelX = alignedZ * eX;
		double rParallelY = alignedZ * eY;
		double rParallelZ = alignedZ * eZ;

		// vector perpendicular to axis of wire loop, pointing at evaluation pos
		double rPerpX = r0x - rParallelX;
		double rPerpY = r0y - rParallelY;
		double rPerpZ = r0z - rParallelZ;

		// perpendicular distance squared between evalPos and axis of wire loop
		double alignedRSq = rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ;

		// prevent division-by-zero when computing radial unit vector
		// A_phi is zero anyway on-axis --> no contribution expected
		if (alignedRSq > 0.0) {

			// perpendicular distance between evalPos and axis of wire loop
			double alignedR = sqrt(alignedRSq);

			// unit vector in radial direction
			double eRX = rPerpX / alignedR;
			double eRY = rPerpY / alignedR;
			double eRZ = rPerpZ / alignedR;

			// normalized rho component of evaluation location in coordinate system of wire loop
			double rhoP = alignedR / radius;

			// compute tangential component of magnetic vector potential, including current and mu_0
			double aPhi = aPrefactor * circularWireLoop_A_phi(rhoP, zP);

			// compute cross product between e_z and e_rho to get e_phi
			double ePhiX = eRY * eZ - eRZ * eY;
			double ePhiY = eRZ * eX - eRX * eZ;
			double ePhiZ = eRX * eY - eRY * eX;

			// add contribution from wire loop to result
			vectorPotential[3 * idxEval + 0] = aPhi * ePhiX;
			vectorPotential[3 * idxEval + 1] = aPhi * ePhiY;
			vectorPotential[3 * idxEval + 2] = aPhi * ePhiZ;
		}
	}
}

/**
 * Compute the magnetic field of a circular wire loop.
 *
 * @param center  [3: x, y, z] origin of loop (in meters)
 * @param normal  [3: x, y, z] normal vector of loop (in meters); will be
 *                normalized internally
 * @param radius  radius of the wire loop (in meters)
 * @param current loop current (in A)
 * @param evalPos [3: x, y, z][nEvalPos] evaluation locations (in meters)
 * @param magneticField [3: B_x, B_y, B_z][nEvalPos] Cartesian components of the magnetic
 *         field evaluated at the given locations (in T); has to be allocated on entry
 */
void magneticFieldCircularFilament(double *center, double *normal, double radius,
		double current, int nEvalPos, double *evalPos, double *magneticField) {

	if (!isfinite(radius) || radius <= 0.0) {
		printf("radius must be finite and positive, but is %g\n", radius);
		return;
	}

	double bPrefactor = MU_0_BY_PI * current / radius;

	// squared length of normal vector
	double nLen2 = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];

	if (nLen2 == 0.0) {
		printf("length of normal vector must not be zero");
		return;
	}

	// length of normal vector
	double nLen = sqrt(nLen2);

	// unit normal vector of wire loop
	double eX = normal[0] / nLen;
	double eY = normal[1] / nLen;
	double eZ = normal[2] / nLen;

	for (int idxEval = 0; idxEval < nEvalPos; ++idxEval) {

		// vector from center of wire loop to eval pos
		double r0x = evalPos[3 * idxEval + 0] - center[0];
		double r0y = evalPos[3 * idxEval + 1] - center[1];
		double r0z = evalPos[3 * idxEval + 2] - center[2];

		// z position along normal of wire loop
		double alignedZ = eX * r0x + eY * r0y + eZ * r0z;

		// normalized z component of evaluation location in coordinate system of wire loop
		double zP = alignedZ / radius;

		// r0 projected onto axis of wire loop
		double rParallelX = alignedZ * eX;
		double rParallelY = alignedZ * eY;
		double rParallelZ = alignedZ * eZ;

		// vector perpendicular to axis of wire loop, pointing at evaluation pos
		double rPerpX = r0x - rParallelX;
		double rPerpY = r0y - rParallelY;
		double rPerpZ = r0z - rParallelZ;

		// perpendicular distance squared between evalPos and axis of wire loop
		double alignedRSq = rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ;

		double rhoP;
		if (alignedRSq > 0.0) {
			// radial unit vector is only defined if evaluation pos is off-axis

			// perpendicular distance between evalPos and axis of wire loop
			double alignedR = sqrt(alignedRSq);

			// unit vector in radial direction
			double eRX = rPerpX / alignedR;
			double eRY = rPerpY / alignedR;
			double eRZ = rPerpZ / alignedR;

			// normalized rho component of evaluation location in coordinate system of wire loop
			rhoP = alignedR / radius;

			// compute radial component of normalized magnetic field
			// and scale by current and mu_0
			double bRho = bPrefactor * circularWireLoop_B_rho(rhoP, zP);

			// add contribution from B_rho of wire loop to result
			magneticField[3 * idxEval + 0] = bRho * eRX;
			magneticField[3 * idxEval + 1] = bRho * eRY;
			magneticField[3 * idxEval + 2] = bRho * eRZ;
		} else {
			rhoP = 0.0;
		}

		// compute vertical component of normalized magnetic field
		// and scale by current and mu_0
		double bZ = bPrefactor * circularWireLoop_B_z(rhoP, zP);

		// add contribution from B_z of wire loop to result
		magneticField[3 * idxEval + 0] += bZ * eX;
		magneticField[3 * idxEval + 1] += bZ * eY;
		magneticField[3 * idxEval + 2] += bZ * eZ;
	}
}

// --------------------------------------------------

/**
 * Compute the magnetic vector potential of a polygon filament
 * at a number of evaluation locations.
 *
 * @param vertices [3: x, y, z][numVertices] points along polygon; in m
 * @param current current along polygon; in A
 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
 * @param vectorPotential [3: x, y, z][numEvalPos] target array for magnetic vector potential at evaluation locations; in Tm
 * @param idxSourceStart first index in {@code vertices} to take into account
 * @param idxSourceEnd (last+1) index in {@code vertices} to take into account
 * @param idxEvalStart first index in {@code evalPos} to take into account
 * @param idxEvalEnd (last+1) index in {@code evalPos} to take into account
 * @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
 *                                of the contributions from the polygon vertices; otherwise, use standard += summation
 */
void kernelVectorPotentialPolygonFilament(
		double *vertices,
		double current,
		double *evalPos,
		double *vectorPotential,
		int idxSourceStart,
		int idxSourceEnd,
		int idxEvalStart,
		int idxEvalEnd,
		bool useCompensatedSummation) {

	double aPrefactor = MU_0_BY_2_PI * current;

	// setup compensated summation objects
	double *aXSum;
	double *aYSum;
	double *aZSum;
	if (useCompensatedSummation) {
		int numEvalPos = idxEvalEnd - idxEvalStart;

		// need three doubles (s, cs, ccs) per eval pos --> see compsum.h
		int numBytesToAllocate = 3 * numEvalPos * sizeof(double);

		aXSum = (double *) malloc (numBytesToAllocate);
		aYSum = (double *) malloc (numBytesToAllocate);
		aZSum = (double *) malloc (numBytesToAllocate);

		// initialize target array to zero
		memset(aXSum, 0, numBytesToAllocate);
		memset(aYSum, 0, numBytesToAllocate);
		memset(aZSum, 0, numBytesToAllocate);
	} else {
		aXSum = NULL;
		aYSum = NULL;
		aZSum = NULL;

		// initialize target array to zero
		for (int idxEval = idxEvalStart; idxEval < idxEvalEnd; ++idxEval) {
			memset(vectorPotential + 3 * idxEval, 0, 3 * sizeof(double));
		}
	}

	double x_i = vertices[3 * idxSourceStart + 0];
	double y_i = vertices[3 * idxSourceStart + 1];
	double z_i = vertices[3 * idxSourceStart + 2];

	for (int idxSource = idxSourceStart; idxSource < idxSourceEnd; ++idxSource) {

		double x_f = vertices[3 * (idxSource + 1) + 0];
		double y_f = vertices[3 * (idxSource + 1) + 1];
		double z_f = vertices[3 * (idxSource + 1) + 2];

		// vector from start to end of i:th wire segment
		double dx = x_f - x_i;
		double dy = y_f - y_i;
		double dz = z_f - z_i;

		// squared length of wire segment
		double l2 = dx * dx + dy * dy + dz * dz;
		if (l2 == 0.0) {
			// skip zero-length segments: no contribution
			continue;
		}

		// length of wire segment
		double l = sqrt(l2);

		// unit vector parallel to wire segment
		double eX = dx / l;
		double eY = dy / l;
		double eZ = dz / l;

		for (int idxEval = idxEvalStart; idxEval < idxEvalEnd; ++idxEval) {

			// vector from start of wire segment to eval pos
			double r0x = evalPos[3 * idxEval + 0] - x_i;
			double r0y = evalPos[3 * idxEval + 1] - y_i;
			double r0z = evalPos[3 * idxEval + 2] - z_i;

			// z position along axis of wire segment
			double alignedZ = eX * r0x + eY * r0y + eZ * r0z;

			// normalized z component of evaluation location in coordinate system of wire segment
			double zP = alignedZ / l;

			// vector perpendicular to axis of wire segment, pointing at evaluation pos
			double rPerpX = r0x - alignedZ * eX;
			double rPerpY = r0y - alignedZ * eY;
			double rPerpZ = r0z - alignedZ * eZ;

			// perpendicular distance between evalPos and axis of wire segment
			double alignedR = sqrt(rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ);

			// normalized rho component of evaluation location in coordinate system of wire segment
			double rhoP = alignedR / l;

			// compute parallel component of magnetic vector potential, including current and mu_0
			double aParallel = aPrefactor * straightWireSegment_A_z(rhoP, zP);

			// add contribution from wire segment to result
			if (useCompensatedSummation) {
				compAdd(aParallel * eX, aXSum + 3 * (idxEval - idxEvalStart));
				compAdd(aParallel * eY, aYSum + 3 * (idxEval - idxEvalStart));
				compAdd(aParallel * eZ, aZSum + 3 * (idxEval - idxEvalStart));
			} else {
				vectorPotential[3 * idxEval + 0] += aParallel * eX;
				vectorPotential[3 * idxEval + 1] += aParallel * eY;
				vectorPotential[3 * idxEval + 2] += aParallel * eZ;
			}
		}

		// shift to next point
		x_i = x_f;
		y_i = y_f;
		z_i = z_f;
	}

	if (useCompensatedSummation) {
		// obtain compensated sums from summation objects
		for (int idxEval = idxEvalStart; idxEval < idxEvalEnd; ++idxEval) {
			int relIdx = 3 * (idxEval - idxEvalStart);
			vectorPotential[3 * idxEval + 0] = aXSum[relIdx + 0] + aXSum[relIdx + 1] + aXSum[relIdx + 2];
			vectorPotential[3 * idxEval + 1] = aYSum[relIdx + 0] + aYSum[relIdx + 1] + aYSum[relIdx + 2];
			vectorPotential[3 * idxEval + 2] = aZSum[relIdx + 0] + aZSum[relIdx + 1] + aZSum[relIdx + 2];
		}

		free(aXSum);
		free(aYSum);
		free(aZSum);
	}
}

/**
 * Compute the magnetic vector potential of a polygon filament
 * at a number of evaluation locations.
 *
 * @param void (*vertexSupplier)(int i, double *point): callback to put i-th current carrier polygon vertex into point as [3: x, y, z]; in m
 * @param current current along polygon; in A
 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
 * @param vectorPotential [3: x, y, z][numEvalPos] target array for magnetic vector potential at evaluation locations; in Tm
 * @param idxSourceStart first index in {@code vertices} to take into account
 * @param idxSourceEnd (last+1) index in {@code vertices} to take into account
 * @param idxEvalStart first index in {@code evalPos} to take into account
 * @param idxEvalEnd (last+1) index in {@code evalPos} to take into account
 * @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
 *                                of the contributions from the polygon vertices; otherwise, use standard += summation
 */
void kernelVectorPotentialPolygonFilamentVertexSupplier(
		void (*vertexSupplier)(int i, double *point),
		double current,
		double *evalPos,
		double *vectorPotential,
		int idxSourceStart,
		int idxSourceEnd,
		int idxEvalStart,
		int idxEvalEnd,
		bool useCompensatedSummation) {

	double aPrefactor = MU_0_BY_2_PI * current;

	// setup compensated summation objects
	double *aXSum;
	double *aYSum;
	double *aZSum;
	if (useCompensatedSummation) {
		int numEvalPos = idxEvalEnd - idxEvalStart;

		// need three doubles (s, cs, ccs) per eval pos --> see compsum.h
		int numBytesToAllocate = 3 * numEvalPos * sizeof(double);

		aXSum = (double *) malloc (numBytesToAllocate);
		aYSum = (double *) malloc (numBytesToAllocate);
		aZSum = (double *) malloc (numBytesToAllocate);

		// initialize target array to zero
		memset(aXSum, 0, numBytesToAllocate);
		memset(aYSum, 0, numBytesToAllocate);
		memset(aZSum, 0, numBytesToAllocate);
	} else {
		aXSum = NULL;
		aYSum = NULL;
		aZSum = NULL;

		// initialize target array to zero
		for (int idxEval = idxEvalStart; idxEval < idxEvalEnd; ++idxEval) {
			memset(vectorPotential + 3 * idxEval, 0, 3 * sizeof(double));
		}
	}

	// get first point from vertexSupplier
	double *pointData = (double *) malloc(3 * sizeof(double));
	vertexSupplier(idxSourceStart, pointData);
	double x_i = pointData[0];
	double y_i = pointData[1];
	double z_i = pointData[2];

	for (int idxSource = idxSourceStart; idxSource < idxSourceEnd; ++idxSource) {

		// get next point from vertexSupplier
		vertexSupplier(idxSource+1, pointData);
		double x_f = pointData[0];
		double y_f = pointData[1];
		double z_f = pointData[2];

		// vector from start to end of i:th wire segment
		double dx = x_f - x_i;
		double dy = y_f - y_i;
		double dz = z_f - z_i;

		// squared length of wire segment
		double l2 = dx * dx + dy * dy + dz * dz;
		if (l2 == 0.0) {
			// skip zero-length segments: no contribution
			continue;
		}

		// length of wire segment
		double l = sqrt(l2);

		// unit vector parallel to wire segment
		double eX = dx / l;
		double eY = dy / l;
		double eZ = dz / l;

		for (int idxEval = idxEvalStart; idxEval < idxEvalEnd; ++idxEval) {

			// vector from start of wire segment to eval pos
			double r0x = evalPos[3 * idxEval + 0] - x_i;
			double r0y = evalPos[3 * idxEval + 1] - y_i;
			double r0z = evalPos[3 * idxEval + 2] - z_i;

			// z position along axis of wire segment
			double alignedZ = eX * r0x + eY * r0y + eZ * r0z;

			// normalized z component of evaluation location in coordinate system of wire segment
			double zP = alignedZ / l;

			// vector perpendicular to axis of wire segment, pointing at evaluation pos
			double rPerpX = r0x - alignedZ * eX;
			double rPerpY = r0y - alignedZ * eY;
			double rPerpZ = r0z - alignedZ * eZ;

			// perpendicular distance between evalPos and axis of wire segment
			double alignedR = sqrt(rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ);

			// normalized rho component of evaluation location in coordinate system of wire segment
			double rhoP = alignedR / l;

			// compute parallel component of magnetic vector potential, including current and mu_0
			double aParallel = aPrefactor * straightWireSegment_A_z(rhoP, zP);

			// add contribution from wire segment to result
			if (useCompensatedSummation) {
				compAdd(aParallel * eX, aXSum + 3 * (idxEval - idxEvalStart));
				compAdd(aParallel * eY, aYSum + 3 * (idxEval - idxEvalStart));
				compAdd(aParallel * eZ, aZSum + 3 * (idxEval - idxEvalStart));
			} else {
				vectorPotential[3 * idxEval + 0] += aParallel * eX;
				vectorPotential[3 * idxEval + 1] += aParallel * eY;
				vectorPotential[3 * idxEval + 2] += aParallel * eZ;
			}
		}

		// shift to next point
		x_i = x_f;
		y_i = y_f;
		z_i = z_f;
	}

	// return target for vertexSupplier not needed anymore
	free(pointData);

	if (useCompensatedSummation) {
		// obtain compensated sums from summation objects
		for (int idxEval = idxEvalStart; idxEval < idxEvalEnd; ++idxEval) {
			int relIdx = 3 * (idxEval - idxEvalStart);
			vectorPotential[3 * idxEval + 0] = aXSum[relIdx + 0] + aXSum[relIdx + 1] + aXSum[relIdx + 2];
			vectorPotential[3 * idxEval + 1] = aYSum[relIdx + 0] + aYSum[relIdx + 1] + aYSum[relIdx + 2];
			vectorPotential[3 * idxEval + 2] = aZSum[relIdx + 0] + aZSum[relIdx + 1] + aZSum[relIdx + 2];
		}

		free(aXSum);
		free(aYSum);
		free(aZSum);
	}
}

/**
 * Compute the magnetic field of a polygon filament
 * at a number of evaluation locations.
 *
 * @param vertices [3: x, y, z][numVertices] points along polygon; in m
 * @param current current along polygon; in A
 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
 * @param magneticField [3: x, y, z][numEvalPos] target array for magnetic field at evaluation locations; in T
 * @param idxSourceStart first index in {@code vertices} to take into account
 * @param idxSourceEnd (last+1) index in {@code vertices} to take into account
 * @param idxEvalStart first index in {@code evalPos} to take into account
 * @param idxEvalEnd (last+1) index in {@code evalPos} to take into account
 * @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
 *                                of the contributions from the polygon vertices; otherwise, use standard += summation
 */
void kernelMagneticFieldPolygonFilament(
		double *vertices,
		double current,
		double *evalPos,
		double *magneticField,
		int idxSourceStart,
		int idxSourceEnd,
		int idxEvalStart,
		int idxEvalEnd,
		bool useCompensatedSummation) {

	// needs additional division by length of wire segment!
	double bPrefactorL = MU_0_BY_4_PI * current;

	// setup compensated summation targets
	double *bXSum;
	double *bYSum;
	double *bZSum;
	if (useCompensatedSummation) {
		int numEvalPos = idxEvalEnd - idxEvalStart;

		// need three doubles (s, cs, ccs) per eval pos --> see compsum.h
		int numBytesToAllocate = 3 * numEvalPos * sizeof(double);

		bXSum = (double *) malloc (numBytesToAllocate);
		bYSum = (double *) malloc (numBytesToAllocate);
		bZSum = (double *) malloc (numBytesToAllocate);

		// initialize target array to zero
		memset(bXSum, 0, numBytesToAllocate);
		memset(bYSum, 0, numBytesToAllocate);
		memset(bZSum, 0, numBytesToAllocate);
	} else {
		bXSum = NULL;
		bYSum = NULL;
		bZSum = NULL;

		// initialize target array to zero
		for (int idxEval = idxEvalStart; idxEval < idxEvalEnd; ++idxEval) {
			memset(magneticField + 3 * idxEval, 0, 3 * sizeof(double));
		}
	}

	double x_i = vertices[3 * idxSourceStart + 0];
	double y_i = vertices[3 * idxSourceStart + 1];
	double z_i = vertices[3 * idxSourceStart + 2];

	for (int idxSource = idxSourceStart; idxSource < idxSourceEnd; ++idxSource) {

		double x_f = vertices[3 * (idxSource + 1) + 0];
		double y_f = vertices[3 * (idxSource + 1) + 1];
		double z_f = vertices[3 * (idxSource + 1) + 2];

		// vector from start to end of i:th wire segment
		double dx = x_f - x_i;
		double dy = y_f - y_i;
		double dz = z_f - z_i;

		// squared length of wire segment
		double l2 = dx * dx + dy * dy + dz * dz;
		if (l2 == 0.0) {
			// skip zero-length segments: no contribution
			continue;
		}

		// length of wire segment
		double l = sqrt(l2);

		// assemble full prefactor for B_phi
		double bPrefactor = bPrefactorL / l;

		// unit vector parallel to wire segment
		double eX = dx / l;
		double eY = dy / l;
		double eZ = dz / l;

		for (int idxEval = idxEvalStart; idxEval < idxEvalEnd; ++idxEval) {

			// vector from start of wire segment to eval pos
			double r0x = evalPos[3 * idxEval + 0] - x_i;
			double r0y = evalPos[3 * idxEval + 1] - y_i;
			double r0z = evalPos[3 * idxEval + 2] - z_i;

			// z position along axis of wire segment
			double alignedZ = eX * r0x + eY * r0y + eZ * r0z;

			// vector perpendicular to axis of wire segment, pointing at evaluation pos
			double rPerpX = r0x - alignedZ * eX;
			double rPerpY = r0y - alignedZ * eY;
			double rPerpZ = r0z - alignedZ * eZ;

			// perpendicular distance squared between evalPos and axis of wire segment
			double alignedRSq = rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ;

			// B_phi is zero along axis of filament
			if (alignedRSq > 0.0) {

				// perpendicular distance between evalPos and axis of wire segment
				double alignedR = sqrt(alignedRSq);

				// normalized rho component of evaluation location in coordinate system of wire segment
				double rhoP = alignedR / l;

				// normalized z component of evaluation location in coordinate system of wire segment
				double zP = alignedZ / l;

				// compute tangential component of magnetic vector potential, including current and mu_0
				double bPhi = bPrefactor * straightWireSegment_B_phi(rhoP, zP);

				// unit vector in radial direction
				double eRX = rPerpX / alignedR;
				double eRY = rPerpY / alignedR;
				double eRZ = rPerpZ / alignedR;

				// compute cross product between e_z and e_rho to get e_phi
				double ePhiX = eY * eRZ - eZ * eRY;
				double ePhiY = eZ * eRX - eX * eRZ;
				double ePhiZ = eX * eRY - eY * eRX;

				// add contribution from wire segment to result
				if (useCompensatedSummation) {
					compAdd(bPhi * ePhiX, bXSum + 3 * (idxEval - idxEvalStart));
					compAdd(bPhi * ePhiY, bYSum + 3 * (idxEval - idxEvalStart));
					compAdd(bPhi * ePhiZ, bZSum + 3 * (idxEval - idxEvalStart));
				} else {
					magneticField[3 * idxEval + 0] += bPhi * ePhiX;
					magneticField[3 * idxEval + 1] += bPhi * ePhiY;
					magneticField[3 * idxEval + 2] += bPhi * ePhiZ;
				}
			}
		}

		// shift to next point
		x_i = x_f;
		y_i = y_f;
		z_i = z_f;
	}

	if (useCompensatedSummation) {
		// obtain compensated sums from summation objects
		for (int idxEval = idxEvalStart; idxEval < idxEvalEnd; ++idxEval) {
			int relIdx = 3 * (idxEval - idxEvalStart);
			magneticField[3 * idxEval + 0] = bXSum[relIdx + 0] + bXSum[relIdx + 1] + bXSum[relIdx + 2];
			magneticField[3 * idxEval + 1] = bYSum[relIdx + 0] + bYSum[relIdx + 1] + bYSum[relIdx + 2];
			magneticField[3 * idxEval + 2] = bZSum[relIdx + 0] + bZSum[relIdx + 1] + bZSum[relIdx + 2];
		}

		free(bXSum);
		free(bYSum);
		free(bZSum);
	}
}

/**
 * Compute the magnetic field of a polygon filament
 * at a number of evaluation locations.
 *
 * @param void (*vertexSupplier)(int i, double *point): callback to put i-th current carrier polygon vertex into point as [3: x, y, z]; in m
 * @param current current along polygon; in A
 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
 * @param magneticField [3: x, y, z][numEvalPos] target array for magnetic field at evaluation locations; in T
 * @param idxSourceStart first index in {@code vertices} to take into account
 * @param idxSourceEnd (last+1) index in {@code vertices} to take into account
 * @param idxEvalStart first index in {@code evalPos} to take into account
 * @param idxEvalEnd (last+1) index in {@code evalPos} to take into account
 * @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
 *                                of the contributions from the polygon vertices; otherwise, use standard += summation
 */
void kernelMagneticFieldPolygonFilamentVertexSupplier(
		void (*vertexSupplier)(int i, double *point),
		double current,
		double *evalPos,
		double *magneticField,
		int idxSourceStart,
		int idxSourceEnd,
		int idxEvalStart,
		int idxEvalEnd,
		bool useCompensatedSummation) {

	// needs additional division by length of wire segment!
	double bPrefactorL = MU_0_BY_4_PI * current;

	// setup compensated summation targets
	double *bXSum;
	double *bYSum;
	double *bZSum;
	if (useCompensatedSummation) {
		int numEvalPos = idxEvalEnd - idxEvalStart;

		// need three doubles (s, cs, ccs) per eval pos --> see compsum.h
		int numBytesToAllocate = 3 * numEvalPos * sizeof(double);

		bXSum = (double *) malloc (numBytesToAllocate);
		bYSum = (double *) malloc (numBytesToAllocate);
		bZSum = (double *) malloc (numBytesToAllocate);

		// initialize target array to zero
		memset(bXSum, 0, numBytesToAllocate);
		memset(bYSum, 0, numBytesToAllocate);
		memset(bZSum, 0, numBytesToAllocate);
	} else {
		bXSum = NULL;
		bYSum = NULL;
		bZSum = NULL;

		// initialize target array to zero
		for (int idxEval = idxEvalStart; idxEval < idxEvalEnd; ++idxEval) {
			memset(magneticField + 3 * idxEval, 0, 3 * sizeof(double));
		}
	}

	// get first point from vertexSupplier
	double *pointData = (double *) malloc(3 * sizeof(double));
	vertexSupplier(idxSourceStart, pointData);
	double x_i = pointData[0];
	double y_i = pointData[1];
	double z_i = pointData[2];

	for (int idxSource = idxSourceStart; idxSource < idxSourceEnd; ++idxSource) {

		// get next point from vertexSupplier
		vertexSupplier(idxSource+1, pointData);
		double x_f = pointData[0];
		double y_f = pointData[1];
		double z_f = pointData[2];

		// vector from start to end of i:th wire segment
		double dx = x_f - x_i;
		double dy = y_f - y_i;
		double dz = z_f - z_i;

		// squared length of wire segment
		double l2 = dx * dx + dy * dy + dz * dz;
		if (l2 == 0.0) {
			// skip zero-length segments: no contribution
			continue;
		}

		// length of wire segment
		double l = sqrt(l2);

		// assemble full prefactor for B_phi
		double bPrefactor = bPrefactorL / l;

		// unit vector parallel to wire segment
		double eX = dx / l;
		double eY = dy / l;
		double eZ = dz / l;

		for (int idxEval = idxEvalStart; idxEval < idxEvalEnd; ++idxEval) {

			// vector from start of wire segment to eval pos
			double r0x = evalPos[3 * idxEval + 0] - x_i;
			double r0y = evalPos[3 * idxEval + 1] - y_i;
			double r0z = evalPos[3 * idxEval + 2] - z_i;

			// z position along axis of wire segment
			double alignedZ = eX * r0x + eY * r0y + eZ * r0z;

			// vector perpendicular to axis of wire segment, pointing at evaluation pos
			double rPerpX = r0x - alignedZ * eX;
			double rPerpY = r0y - alignedZ * eY;
			double rPerpZ = r0z - alignedZ * eZ;

			// perpendicular distance squared between evalPos and axis of wire segment
			double alignedRSq = rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ;

			// B_phi is zero along axis of filament
			if (alignedRSq > 0.0) {

				// perpendicular distance between evalPos and axis of wire segment
				double alignedR = sqrt(alignedRSq);

				// normalized rho component of evaluation location in coordinate system of wire segment
				double rhoP = alignedR / l;

				// normalized z component of evaluation location in coordinate system of wire segment
				double zP = alignedZ / l;

				// compute tangential component of magnetic vector potential, including current and mu_0
				double bPhi = bPrefactor * straightWireSegment_B_phi(rhoP, zP);

				// unit vector in radial direction
				double eRX = rPerpX / alignedR;
				double eRY = rPerpY / alignedR;
				double eRZ = rPerpZ / alignedR;

				// compute cross product between e_z and e_rho to get e_phi
				double ePhiX = eY * eRZ - eZ * eRY;
				double ePhiY = eZ * eRX - eX * eRZ;
				double ePhiZ = eX * eRY - eY * eRX;

				// add contribution from wire segment to result
				if (useCompensatedSummation) {
					compAdd(bPhi * ePhiX, bXSum + 3 * (idxEval - idxEvalStart));
					compAdd(bPhi * ePhiY, bYSum + 3 * (idxEval - idxEvalStart));
					compAdd(bPhi * ePhiZ, bZSum + 3 * (idxEval - idxEvalStart));
				} else {
					magneticField[3 * idxEval + 0] += bPhi * ePhiX;
					magneticField[3 * idxEval + 1] += bPhi * ePhiY;
					magneticField[3 * idxEval + 2] += bPhi * ePhiZ;
				}
			}
		}

		// shift to next point
		x_i = x_f;
		y_i = y_f;
		z_i = z_f;
	}

	// return target for vertexSupplier not needed anymore
	free(pointData);

	if (useCompensatedSummation) {
		// obtain compensated sums from summation objects
		for (int idxEval = idxEvalStart; idxEval < idxEvalEnd; ++idxEval) {
			int relIdx = 3 * (idxEval - idxEvalStart);
			magneticField[3 * idxEval + 0] = bXSum[relIdx + 0] + bXSum[relIdx + 1] + bXSum[relIdx + 2];
			magneticField[3 * idxEval + 1] = bYSum[relIdx + 0] + bYSum[relIdx + 1] + bYSum[relIdx + 2];
			magneticField[3 * idxEval + 2] = bZSum[relIdx + 0] + bZSum[relIdx + 1] + bZSum[relIdx + 2];
		}

		free(bXSum);
		free(bYSum);
		free(bZSum);
	}
}

// --------------------------------------------------

/**
 * Compute the magnetic vector potential of a polygon filament
 * at a number of evaluation locations.
 *
 * @param vertices [3: x, y, z][numVertices] points along polygon; in m
 * @param current current along polygon; in A
 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
 * @param vectorPotential [3: x, y, z][numEvalPos] target array for magnetic vector potential at evaluation locations; in Tm
 * @param numProcessors number of processors to use for parallelization
 * @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
 *                                of the contributions from the polygon vertices; otherwise, use standard += summation
 */
void vectorPotentialPolygonFilament_specPar_specSum(
		int numVertices,
		double *vertices,
		double current,
		int numEvalPos,
		double *evalPos,
		double *vectorPotential,
		int numProcessors,
		bool useCompensatedSummation) {

	if (numVertices < 2) {
		printf("need at least 2 vertices, but only got %d\n", numVertices);
		return;
	}

	if (numProcessors < 1) {
		printf("need at least 1 processor, but only got %d\n", numProcessors);
		return;
	}

	if (current == 0.0) {
		// TODO: memset(0, vectorPotential)
		return;
	}

	if (numProcessors == 1) {
		// single-threaded call
		int idxSourceStart = 0;
		int idxSourceEnd   = numVertices-1;
		int idxEvalStart   = 0;
		int idxEvalEnd     = numEvalPos;
		kernelVectorPotentialPolygonFilament(
				vertices, current,
				evalPos,
				vectorPotential,
				idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
				useCompensatedSummation);
	} else {
		// use multithreading

		if (numVertices-1 > numEvalPos) {
			// parallelize over nSource-1

			// Note that each thread needs its own copy of the vectorPotential array,
			// so this approach might need quite some memory in case the number of
			// threads and the number of evaluation points is large.

			int nThreads;
			int nSourcePerThread;
			if (numVertices-1 < numProcessors) {
				nThreads = numVertices-1;
				nSourcePerThread = 1;
			} else {
				nThreads = numProcessors;

				// It is better that many threads do more
				// than one thread needs to do more.
				nSourcePerThread = (int) ceil( (numVertices-1.0) / nThreads);
			}

			double *vectorPotentialContributions = (double *) malloc(nThreads * 3 * numEvalPos * sizeof(double));
			if (vectorPotentialContributions == NULL) {
				printf("failed to allocate temporary array for vector potential contributions\n");
				return;
			}

			// parallelized evaluation
			int idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd;
#ifdef _OPENMP
#pragma omp parallel for private(idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd)
#endif // _OPENMP
			for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
				idxSourceStart =      idxThread    * nSourcePerThread;
				idxSourceEnd   = min((idxThread+1) * nSourcePerThread, numVertices-1);
				idxEvalStart   = 0;
				idxEvalEnd     = numEvalPos;

				kernelVectorPotentialPolygonFilament(
						vertices, current,
						evalPos,
						vectorPotentialContributions + idxThread * 3 * numEvalPos,
						idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
						useCompensatedSummation);
			}

			// sum up contributions from source chunks
			if (useCompensatedSummation) {
				double *sumX = (double *) malloc(3 * sizeof(double));
				double *sumY = (double *) malloc(3 * sizeof(double));
				double *sumZ = (double *) malloc(3 * sizeof(double));
				for (int i=0; i<numEvalPos; ++i) {
					memset(sumX, 0, 3 * sizeof(*sumX));
					memset(sumY, 0, 3 * sizeof(*sumY));
					memset(sumZ, 0, 3 * sizeof(*sumZ));
					// TODO: bad memory access pattern here --> potential bottleneck !!!
					for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
						compAdd(vectorPotentialContributions[idxThread * 3 * numEvalPos + 3 * i + 0], sumX);
						compAdd(vectorPotentialContributions[idxThread * 3 * numEvalPos + 3 * i + 1], sumY);
						compAdd(vectorPotentialContributions[idxThread * 3 * numEvalPos + 3 * i + 2], sumZ);
					}
					vectorPotential[3 * i + 0] = sumX[0] + sumX[1] + sumX[2];
					vectorPotential[3 * i + 1] = sumY[0] + sumY[1] + sumY[2];
					vectorPotential[3 * i + 2] = sumZ[0] + sumZ[1] + sumZ[2];
				}
				free(sumX);
				free(sumY);
				free(sumZ);
			} else {
				// TODO: memset(vectorPotential, 0)
				for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
					for (int i=0; i<numEvalPos; ++i) {
						vectorPotential[3 * i + 0] += vectorPotentialContributions[idxThread * 3 * numEvalPos + 3 * i + 0];
						vectorPotential[3 * i + 1] += vectorPotentialContributions[idxThread * 3 * numEvalPos + 3 * i + 1];
						vectorPotential[3 * i + 2] += vectorPotentialContributions[idxThread * 3 * numEvalPos + 3 * i + 2];
					}
				}
			}

			free(vectorPotentialContributions);
		} else { // nEval > nSource
			// parallelize over nEval

			int nThreads;
			int nEvalPerThread;
			if (numEvalPos < numProcessors) {
				nThreads = numEvalPos;
				nEvalPerThread = 1;
			} else {
				nThreads = numProcessors;
				nEvalPerThread = (int) ceil( ((double) numEvalPos) / nThreads );
			}

			// parallelized evaluation
			int idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd;
#ifdef _OPENMP
#pragma omp parallel for private(idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd)
#endif // _OPENMP
			for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
				idxSourceStart = 0;
				idxSourceEnd   = numVertices-1;
				idxEvalStart   =      idxThread    * nEvalPerThread;
				idxEvalEnd     = min((idxThread+1) * nEvalPerThread, numEvalPos);

				kernelVectorPotentialPolygonFilament(
						vertices, current,
						evalPos,
						vectorPotential,
						idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
						useCompensatedSummation);
			}
		} // parallelize over nSource or nEval
	} // parallelization
}

/**
 * Compute the magnetic vector potential of a polygon filament
 * at a number of evaluation locations.
 *
 * @param void (*vertexSupplier)(int i, double *point): callback to put i-th current carrier polygon vertex into point as [3: x, y, z]; in m
 * @param current current along polygon; in A
 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
 * @param vectorPotential [3: x, y, z][numEvalPos] target array for magnetic vector potential at evaluation locations; in Tm
 * @param numProcessors number of processors to use for parallelization
 * @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
 *                                of the contributions from the polygon vertices; otherwise, use standard += summation
 */
void vectorPotentialPolygonFilamentVertexSupplier_specPar_specSum(
		int numVertices,
		void (*vertexSupplier)(int i, double *point),
		double current,
		int numEvalPos,
		double *evalPos,
		double *vectorPotential,
		int numProcessors,
		bool useCompensatedSummation) {

	if (numVertices < 2) {
		printf("need at least 2 vertices, but only got %d\n", numVertices);
		return;
	}

	if (numProcessors < 1) {
		printf("need at least 1 processor, but only got %d\n", numProcessors);
		return;
	}

	if (current == 0.0) {
		// TODO: memset(0, vectorPotential)
		return;
	}

	if (numProcessors == 1) {
		// single-threaded call
		int idxSourceStart = 0;
		int idxSourceEnd   = numVertices-1;
		int idxEvalStart   = 0;
		int idxEvalEnd     = numEvalPos;
		kernelVectorPotentialPolygonFilamentVertexSupplier(
				vertexSupplier, current,
				evalPos,
				vectorPotential,
				idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
				useCompensatedSummation);
	} else {
		// use multithreading

		if (numVertices-1 > numEvalPos) {
			// parallelize over nSource-1

			// Note that each thread needs its own copy of the vectorPotential array,
			// so this approach might need quite some memory in case the number of
			// threads and the number of evaluation points is large.

			int nThreads;
			int nSourcePerThread;
			if (numVertices-1 < numProcessors) {
				nThreads = numVertices-1;
				nSourcePerThread = 1;
			} else {
				nThreads = numProcessors;

				// It is better that many threads do more
				// than one thread needs to do more.
				nSourcePerThread = (int) ceil( (numVertices-1.0) / nThreads);
			}

			double *vectorPotentialContributions = (double *) malloc(nThreads * 3 * numEvalPos * sizeof(double));
			if (vectorPotentialContributions == NULL) {
				printf("failed to allocate temporary array for vector potential contributions\n");
				return;
			}

			// parallelized evaluation
			int idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd;
#ifdef _OPENMP
#pragma omp parallel for private(idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd)
#endif // _OPENMP
			for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
				idxSourceStart =      idxThread    * nSourcePerThread;
				idxSourceEnd   = min((idxThread+1) * nSourcePerThread, numVertices-1);
				idxEvalStart   = 0;
				idxEvalEnd     = numEvalPos;

				kernelVectorPotentialPolygonFilamentVertexSupplier(
						vertexSupplier, current,
						evalPos,
						vectorPotentialContributions + idxThread * 3 * numEvalPos,
						idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
						useCompensatedSummation);
			}

			// sum up contributions from source chunks
			if (useCompensatedSummation) {
				double *sumX = (double *) malloc(3 * sizeof(double));
				double *sumY = (double *) malloc(3 * sizeof(double));
				double *sumZ = (double *) malloc(3 * sizeof(double));
				for (int i=0; i<numEvalPos; ++i) {
					memset(sumX, 0, 3 * sizeof(*sumX));
					memset(sumY, 0, 3 * sizeof(*sumY));
					memset(sumZ, 0, 3 * sizeof(*sumZ));
					// TODO: bad memory access pattern here --> potential bottleneck !!!
					for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
						compAdd(vectorPotentialContributions[idxThread * 3 * numEvalPos + 3 * i + 0], sumX);
						compAdd(vectorPotentialContributions[idxThread * 3 * numEvalPos + 3 * i + 1], sumY);
						compAdd(vectorPotentialContributions[idxThread * 3 * numEvalPos + 3 * i + 2], sumZ);
					}
					vectorPotential[3 * i + 0] = sumX[0] + sumX[1] + sumX[2];
					vectorPotential[3 * i + 1] = sumY[0] + sumY[1] + sumY[2];
					vectorPotential[3 * i + 2] = sumZ[0] + sumZ[1] + sumZ[2];
				}
				free(sumX);
				free(sumY);
				free(sumZ);
			} else {
				// TODO: memset(vectorPotential, 0)
				for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
					for (int i=0; i<numEvalPos; ++i) {
						vectorPotential[3 * i + 0] += vectorPotentialContributions[idxThread * 3 * numEvalPos + 3 * i + 0];
						vectorPotential[3 * i + 1] += vectorPotentialContributions[idxThread * 3 * numEvalPos + 3 * i + 1];
						vectorPotential[3 * i + 2] += vectorPotentialContributions[idxThread * 3 * numEvalPos + 3 * i + 2];
					}
				}
			}

			free(vectorPotentialContributions);
		} else { // nEval > nSource
			// parallelize over nEval

			int nThreads;
			int nEvalPerThread;
			if (numEvalPos < numProcessors) {
				nThreads = numEvalPos;
				nEvalPerThread = 1;
			} else {
				nThreads = numProcessors;
				nEvalPerThread = (int) ceil( ((double) numEvalPos) / nThreads );
			}

			// parallelized evaluation
			int idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd;
#ifdef _OPENMP
#pragma omp parallel for private(idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd)
#endif // _OPENMP
			for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
				idxSourceStart = 0;
				idxSourceEnd   = numVertices-1;
				idxEvalStart   =      idxThread    * nEvalPerThread;
				idxEvalEnd     = min((idxThread+1) * nEvalPerThread, numEvalPos);

				kernelVectorPotentialPolygonFilamentVertexSupplier(
						vertexSupplier, current,
						evalPos,
						vectorPotential,
						idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
						useCompensatedSummation);
			}
		} // parallelize over nSource or nEval
	} // parallelization
}

/**
 * Compute the magnetic field of a polygon filament
 * at a number of evaluation locations.
 *
 * @param vertices [3: x, y, z][numVertices] points along polygon; in m
 * @param current current along polygon; in A
 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
 * @param magneticField [3: x, y, z][numEvalPos] target array for magnetic field at evaluation locations; inT
 * @param numProcessors number of processors to use for parallelization
 * @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
 *                                of the contributions from the polygon vertices; otherwise, use standard += summation
 */
void magneticFieldPolygonFilament_specPar_specSum(
		int numVertices,
		double *vertices,
		double current,
		int numEvalPos,
		double *evalPos,
		double *magneticField,
		int numProcessors,
		bool useCompensatedSummation) {

	if (numVertices < 2) {
		printf("need at least 2 vertices, but only got %d\n", numVertices);
		return;
	}

	if (numProcessors < 1) {
		printf("need at least 1 processor, but only got %d\n", numProcessors);
		return;
	}

	if (current == 0.0) {
		// TODO: memset(0, magneticField)
		return;
	}

	if (numProcessors == 1) {
		// single-threaded call
		int idxSourceStart = 0;
		int idxSourceEnd   = numVertices-1;
		int idxEvalStart   = 0;
		int idxEvalEnd     = numEvalPos;
		kernelMagneticFieldPolygonFilament(
				vertices, current,
				evalPos,
				magneticField,
				idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
				useCompensatedSummation);
	} else {
		// use multithreading

		if (numVertices-1 > numEvalPos) {
			// parallelize over nSource-1

			// Note that each thread needs its own copy of the vectorPotential array,
			// so this approach might need quite some memory in case the number of
			// threads and the number of evaluation points is large.

			int nThreads;
			int nSourcePerThread;
			if (numVertices-1 < numProcessors) {
				nThreads = numVertices-1;
				nSourcePerThread = 1;
			} else {
				nThreads = numProcessors;

				// It is better that many threads do more
				// than one thread needs to do more.
				nSourcePerThread = (int) ceil( (numVertices-1.0) / nThreads);
			}

			double *magneticFieldContributions = (double *) malloc(nThreads * 3 * numEvalPos * sizeof(double));
			if (magneticFieldContributions == NULL) {
				printf("failed to allocate temporary array for magnetic field contributions\n");
				return;
			}

			// parallelized evaluation
			int idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd;
#ifdef _OPENMP
#pragma omp parallel for private(idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd)
#endif // _OPENMP
			for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
				idxSourceStart =      idxThread    * nSourcePerThread;
				idxSourceEnd   = min((idxThread+1) * nSourcePerThread, numVertices-1);
				idxEvalStart   = 0;
				idxEvalEnd     = numEvalPos;

				kernelMagneticFieldPolygonFilament(
						vertices, current,
						evalPos,
						magneticFieldContributions + idxThread * 3 * numEvalPos,
						idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
						useCompensatedSummation);
			}

			// sum up contributions from source chunks
			if (useCompensatedSummation) {
				double *sumX = (double *) malloc(3 * sizeof(double));
				double *sumY = (double *) malloc(3 * sizeof(double));
				double *sumZ = (double *) malloc(3 * sizeof(double));
				for (int i=0; i<numEvalPos; ++i) {
					memset(sumX, 0, 3 * sizeof(*sumX));
					memset(sumY, 0, 3 * sizeof(*sumY));
					memset(sumZ, 0, 3 * sizeof(*sumZ));
					// TODO: bad memory access pattern here --> potential bottleneck !!!
					for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
						compAdd(magneticFieldContributions[idxThread * 3 * numEvalPos + 3 * i + 0], sumX);
						compAdd(magneticFieldContributions[idxThread * 3 * numEvalPos + 3 * i + 1], sumY);
						compAdd(magneticFieldContributions[idxThread * 3 * numEvalPos + 3 * i + 2], sumZ);
					}
					magneticField[3 * i + 0] = sumX[0] + sumX[1] + sumX[2];
					magneticField[3 * i + 1] = sumY[0] + sumY[1] + sumY[2];
					magneticField[3 * i + 2] = sumZ[0] + sumZ[1] + sumZ[2];
				}
				free(sumX);
				free(sumY);
				free(sumZ);
			} else {
				// TODO: memset(magneticField, 0)
				for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
					for (int i=0; i<numEvalPos; ++i) {
						magneticField[3 * i + 0] += magneticFieldContributions[idxThread * 3 * numEvalPos + 3 * i + 0];
						magneticField[3 * i + 1] += magneticFieldContributions[idxThread * 3 * numEvalPos + 3 * i + 1];
						magneticField[3 * i + 2] += magneticFieldContributions[idxThread * 3 * numEvalPos + 3 * i + 2];
					}
				}
			}

			free(magneticFieldContributions);
		} else { // nEval > nSource
			// parallelize over nEval

			int nThreads;
			int nEvalPerThread;
			if (numEvalPos < numProcessors) {
				nThreads = numEvalPos;
				nEvalPerThread = 1;
			} else {
				nThreads = numProcessors;
				nEvalPerThread = (int) ceil( ((double) numEvalPos) / nThreads );
			}

			// parallelized evaluation
			int idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd;
#ifdef _OPENMP
#pragma omp parallel for private(idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd)
#endif // _OPENMP
			for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
				idxSourceStart = 0;
				idxSourceEnd   = numVertices-1;
				idxEvalStart   =      idxThread    * nEvalPerThread;
				idxEvalEnd     = min((idxThread+1) * nEvalPerThread, numEvalPos);

				kernelMagneticFieldPolygonFilament(
						vertices, current,
						evalPos,
						magneticField,
						idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
						useCompensatedSummation);
			}
		} // parallelize over nSource or nEval
	} // parallelization
}


/**
 * Compute the magnetic field of a polygon filament
 * at a number of evaluation locations.
 *
 * @param void (*vertexSupplier)(int i, double *point): callback to put i-th current carrier polygon vertex into point as [3: x, y, z]; in m
 * @param current current along polygon; in A
 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
 * @param magneticField [3: x, y, z][numEvalPos] target array for magnetic field at evaluation locations; in T
 * @param numProcessors number of processors to use for parallelization
 * @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
 *                                of the contributions from the polygon vertices; otherwise, use standard += summation
 */
void magneticFieldPolygonFilamentVertexSupplier_specPar_specSum(
		int numVertices,
		void (*vertexSupplier)(int i, double *point),
		double current,
		int numEvalPos,
		double *evalPos,
		double *magneticField,
		int numProcessors,
		bool useCompensatedSummation) {

	if (numVertices < 2) {
		printf("need at least 2 vertices, but only got %d\n", numVertices);
		return;
	}

	if (numProcessors < 1) {
		printf("need at least 1 processor, but only got %d\n", numProcessors);
		return;
	}

	if (current == 0.0) {
		// TODO: memset(0, magneticField)
		return;
	}

	if (numProcessors == 1) {
		// single-threaded call
		int idxSourceStart = 0;
		int idxSourceEnd   = numVertices-1;
		int idxEvalStart   = 0;
		int idxEvalEnd     = numEvalPos;
		kernelMagneticFieldPolygonFilamentVertexSupplier(
				vertexSupplier, current,
				evalPos,
				magneticField,
				idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
				useCompensatedSummation);
	} else {
		// use multithreading

		if (numVertices-1 > numEvalPos) {
			// parallelize over nSource-1

			// Note that each thread needs its own copy of the magneticField array,
			// so this approach might need quite some memory in case the number of
			// threads and the number of evaluation points is large.

			int nThreads;
			int nSourcePerThread;
			if (numVertices-1 < numProcessors) {
				nThreads = numVertices-1;
				nSourcePerThread = 1;
			} else {
				nThreads = numProcessors;

				// It is better that many threads do more
				// than one thread needs to do more.
				nSourcePerThread = (int) ceil( (numVertices-1.0) / nThreads);
			}

			double *magneticFieldContributions = (double *) malloc(nThreads * 3 * numEvalPos * sizeof(double));
			if (magneticFieldContributions == NULL) {
				printf("failed to allocate temporary array for magnetic field contributions\n");
				return;
			}

			// parallelized evaluation
			int idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd;
#ifdef _OPENMP
#pragma omp parallel for private(idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd)
#endif // _OPENMP
			for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
				idxSourceStart =      idxThread    * nSourcePerThread;
				idxSourceEnd   = min((idxThread+1) * nSourcePerThread, numVertices-1);
				idxEvalStart   = 0;
				idxEvalEnd     = numEvalPos;

				kernelMagneticFieldPolygonFilamentVertexSupplier(
						vertexSupplier, current,
						evalPos,
						magneticFieldContributions + idxThread * 3 * numEvalPos,
						idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
						useCompensatedSummation);
			}

			// sum up contributions from source chunks
			if (useCompensatedSummation) {
				double *sumX = (double *) malloc(3 * sizeof(double));
				double *sumY = (double *) malloc(3 * sizeof(double));
				double *sumZ = (double *) malloc(3 * sizeof(double));
				for (int i=0; i<numEvalPos; ++i) {
					memset(sumX, 0, 3 * sizeof(*sumX));
					memset(sumY, 0, 3 * sizeof(*sumY));
					memset(sumZ, 0, 3 * sizeof(*sumZ));
					// TODO: bad memory access pattern here --> potential bottleneck !!!
					for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
						compAdd(magneticFieldContributions[idxThread * 3 * numEvalPos + 3 * i + 0], sumX);
						compAdd(magneticFieldContributions[idxThread * 3 * numEvalPos + 3 * i + 1], sumY);
						compAdd(magneticFieldContributions[idxThread * 3 * numEvalPos + 3 * i + 2], sumZ);
					}
					magneticField[3 * i + 0] = sumX[0] + sumX[1] + sumX[2];
					magneticField[3 * i + 1] = sumY[0] + sumY[1] + sumY[2];
					magneticField[3 * i + 2] = sumZ[0] + sumZ[1] + sumZ[2];
				}
				free(sumX);
				free(sumY);
				free(sumZ);
			} else {
				// TODO: memset(magneticField, 0)
				for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
					for (int i=0; i<numEvalPos; ++i) {
						magneticField[3 * i + 0] += magneticFieldContributions[3 * (numEvalPos * idxThread + i) + 0];
						magneticField[3 * i + 1] += magneticFieldContributions[3 * (numEvalPos * idxThread + i) + 1];
						magneticField[3 * i + 2] += magneticFieldContributions[3 * (numEvalPos * idxThread + i) + 2];
					}
				}
			}

			free(magneticFieldContributions);
		} else { // nEval > nSource
			// parallelize over nEval

			int nThreads;
			int nEvalPerThread;
			if (numEvalPos < numProcessors) {
				nThreads = numEvalPos;
				nEvalPerThread = 1;
			} else {
				nThreads = numProcessors;
				nEvalPerThread = (int) ceil( ((double) numEvalPos) / nThreads );
			}

			// parallelized evaluation
			int idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd;
#ifdef _OPENMP
#pragma omp parallel for private(idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd)
#endif // _OPENMP
			for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
				idxSourceStart = 0;
				idxSourceEnd   = numVertices-1;
				idxEvalStart   =      idxThread    * nEvalPerThread;
				idxEvalEnd     = min((idxThread+1) * nEvalPerThread, numEvalPos);

				kernelMagneticFieldPolygonFilamentVertexSupplier(
						vertexSupplier, current,
						evalPos,
						magneticField,
						idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
						useCompensatedSummation);
			}
		} // parallelize over nSource or nEval
	} // parallelization
}

// --------------------------------------------------

/**
 * Compute the magnetic vector potential of a polygon filament
 * at a number of evaluation locations.
 * Kahan-Babuska compensated summation is used to compute the superposition
 * of the contributions from the polygon vertices.
 *
 * @param vertices [3: x, y, z][numVertices] points along polygon; in m
 * @param current current along polygon; in A
 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
 * @param vectorPotential [3: x, y, z][numEvalPos] target array for magnetic vector potential at evaluation locations; in Tm
 * @param numProcessors number of processors to use for parallelization
 */
void vectorPotentialPolygonFilament_specPar(
		int numVertices,
		double *vertices,
		double current,
		int numEvalPos,
		double *evalPos,
		double *vectorPotential,
		int numProcessors) {

	bool useCompensatedSummation = true;
	vectorPotentialPolygonFilament_specPar_specSum(
			numVertices, vertices, current,
			numEvalPos, evalPos,
			vectorPotential,
			numProcessors,
			useCompensatedSummation);
}

/**
 * Compute the magnetic vector potential of a polygon filament
 * at a number of evaluation locations.
 * Kahan-Babuska compensated summation is used to compute the superposition
 * of the contributions from the polygon vertices.
 *
 * @param void (*vertexSupplier)(int i, double *point): callback to put i-th current carrier polygon vertex into point as [3: x, y, z]; in m
 * @param current current along polygon; in A
 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
 * @param vectorPotential [3: x, y, z][numEvalPos] target array for magnetic vector potential at evaluation locations; in Tm
 * @param numProcessors number of processors to use for parallelization
 */
void vectorPotentialPolygonFilamentVertexSupplier_specPar(
		int numVertices,
		void (*vertexSupplier)(int i, double *point),
		double current,
		int numEvalPos,
		double *evalPos,
		double *vectorPotential,
		int numProcessors) {

	bool useCompensatedSummation = true;
	vectorPotentialPolygonFilamentVertexSupplier_specPar_specSum(
			numVertices, vertexSupplier, current,
			numEvalPos, evalPos,
			vectorPotential,
			numProcessors,
			useCompensatedSummation);
}

/**
 * Compute the magnetic field of a polygon filament
 * at a number of evaluation locations.
 * Kahan-Babuska compensated summation is used to compute the superposition
 * of the contributions from the polygon vertices.
 *
 * @param vertices [3: x, y, z][numVertices] points along polygon; in m
 * @param current current along polygon; in A
 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
 * @param magneticField [3: x, y, z][numEvalPos] target array for magnetic field at evaluation locations; inT
 * @param numProcessors number of processors to use for parallelization
 */
void magneticFieldPolygonFilament_specPar(
		int numVertices,
		double *vertices,
		double current,
		int numEvalPos,
		double *evalPos,
		double *magneticField,
		int numProcessors) {

	bool useCompensatedSummation = true;
	magneticFieldPolygonFilament_specPar_specSum(
			numVertices, vertices, current,
			numEvalPos, evalPos,
			magneticField,
			numProcessors,
			useCompensatedSummation);
}

/**
 * Compute the magnetic field of a polygon filament
 * at a number of evaluation locations.
 * Kahan-Babuska compensated summation is used to compute the superposition
 * of the contributions from the polygon vertices.
 *
 * @param void (*vertexSupplier)(int i, double *point): callback to put i-th current carrier polygon vertex into point as [3: x, y, z]; in m
 * @param current current along polygon; in A
 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
 * @param magneticField [3: x, y, z][numEvalPos] target array for magnetic field at evaluation locations; inT
 * @param numProcessors number of processors to use for parallelization
 */
void magneticFieldPolygonFilamentVertexSupplier_specPar(
		int numVertices,
		void (*vertexSupplier)(int i, double *point),
		double current,
		int numEvalPos,
		double *evalPos,
		double *magneticField,
		int numProcessors) {

	bool useCompensatedSummation = true;
	magneticFieldPolygonFilamentVertexSupplier_specPar_specSum(
			numVertices, vertexSupplier, current,
			numEvalPos, evalPos,
			magneticField,
			numProcessors,
			useCompensatedSummation);
}

// --------------------------------------------------

/**
 * Compute the magnetic vector potential of a polygon filament
 * at a number of evaluation locations.
 * Kahan-Babuska compensated summation is used to compute the superposition
 * of the contributions from the polygon vertices.
 * The computation is parallelized over all available processors.
 *
 * @param vertices [3: x, y, z][numVertices] points along polygon; in m
 * @param current current along polygon; in A
 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
 * @param vectorPotential [3: x, y, z][numEvalPos] target array for magnetic vector potential at evaluation locations; in Tm
 */
void vectorPotentialPolygonFilament(
		int numVertices,
		double *vertices,
		double current,
		int numEvalPos,
		double *evalPos,
		double *vectorPotential) {

	int numProcessors = omp_get_max_threads();
	vectorPotentialPolygonFilament_specPar(
			numVertices, vertices, current,
			numEvalPos, evalPos,
			vectorPotential,
			numProcessors);
}

/**
 * Compute the magnetic vector potential of a polygon filament
 * at a number of evaluation locations.
 * Kahan-Babuska compensated summation is used to compute the superposition
 * of the contributions from the polygon vertices.
 * The computation is parallelized over all available processors.
 *
 * @param void (*vertexSupplier)(int i, double *point): callback to put i-th current carrier polygon vertex into point as [3: x, y, z]; in m
 * @param current current along polygon; in A
 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
 * @param vectorPotential [3: x, y, z][numEvalPos] target array for magnetic vector potential at evaluation locations; in Tm
 */
void vectorPotentialPolygonFilamentVertexSupplier(
		int numVertices,
		void (*vertexSupplier)(int i, double *point),
		double current,
		int numEvalPos,
		double *evalPos,
		double *vectorPotential) {

	int numProcessors = omp_get_max_threads();
	vectorPotentialPolygonFilamentVertexSupplier_specPar(
			numVertices, vertexSupplier, current,
			numEvalPos, evalPos,
			vectorPotential,
			numProcessors);
}

/**
 * Compute the magnetic field of a polygon filament
 * at a number of evaluation locations.
 * Kahan-Babuska compensated summation is used to compute the superposition
 * of the contributions from the polygon vertices.
 * The computation is parallelized over all available processors.
 *
 * @param vertices [3: x, y, z][numVertices] points along polygon; in m
 * @param current current along polygon; in A
 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
 * @param magneticField [3: x, y, z][numEvalPos] target array for magnetic field at evaluation locations; inT
 */
void magneticFieldPolygonFilament(
		int numVertices,
		double *vertices,
		double current,
		int numEvalPos,
		double *evalPos,
		double *magneticField) {

	int numProcessors = omp_get_max_threads();
	magneticFieldPolygonFilament_specPar(
			numVertices, vertices, current,
			numEvalPos, evalPos,
			magneticField,
			numProcessors);
}

/**
 * Compute the magnetic field of a polygon filament
 * at a number of evaluation locations.
 * Kahan-Babuska compensated summation is used to compute the superposition
 * of the contributions from the polygon vertices.
 * The computation is parallelized over all available processors.
 *
 * @param void (*vertexSupplier)(int i, double *point): callback to put i-th current carrier polygon vertex into point as [3: x, y, z]; in m
 * @param current current along polygon; in A
 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
 * @param magneticField [3: x, y, z][numEvalPos] target array for magnetic field at evaluation locations; inT
 */
void magneticFieldPolygonFilamentVertexSupplier(
		int numVertices,
		void (*vertexSupplier)(int i, double *point),
		double current,
		int numEvalPos,
		double *evalPos,
		double *magneticField) {

	int numProcessors = omp_get_max_threads();
	magneticFieldPolygonFilamentVertexSupplier_specPar(
			numVertices, vertexSupplier, current,
			numEvalPos, evalPos,
			magneticField,
			numProcessors);
}

#endif // ABSCAB_H
