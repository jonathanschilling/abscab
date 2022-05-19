package de.labathome.abscab;

import org.apache.commons.math3.util.FastMath;

import de.labathome.CompleteEllipticIntegral;

/** Accurate Biot-Savart routines with Correct Asymptotic Behavior */
public class ABSCAB {

	/**
	 * Full-field magnetic vector potential of straight wire segment; only A_z component is present.
	 * @param rhoP
	 * @param zP
	 * @return
	 */
	public static double straightWireSegment_A_z(double rhoP, double zP) {
		if (rhoP == 0.0) {
			return A_z_along_rhoP_0(rhoP, zP);
		} else if (zP == 0.0 || zP == 1.0) {
			return A_z_along_zP_0_or_1(rhoP, zP);
		} else if (-1 < zP && zP <= 2.0 && rhoP < 1.0) {
			// near-field
			if (zP >= 1.0 || rhoP/(1-zP) >= 1.0) {
				return A_z_6a(rhoP, zP);
			} else if (zP >= 0.0 && rhoP/zP <= 1.0) {
				return A_z_6b(rhoP, zP);
			} else {
				return A_z_6c(rhoP, zP);
			}
		} else {
			return A_z_1(rhoP, zP);
		}
	}

	/**
	 * Full-field magnetic field of straight wire segment; only B_phi component is present.
	 * @param rhoP
	 * @param zP
	 * @return
	 */
	public static double straightWireSegment_B_phi(double rhoP, double zP) {
		if (rhoP == 0.0) {
			return B_phi_2(rhoP, zP);
		} else if (zP == 0.0 || zP == 1.0) {
			return B_phi_3(rhoP, zP);
		} else if (zP >= 1.0 || zP <= 0.0 || rhoP >= 1.0 ||
				rhoP/(1-zP) >= 1.0 || rhoP/zP >= 1.0) {
			return B_phi_4(rhoP, zP);
		} else {
			return B_phi_5(rhoP, zP);
		}
	}

	/**
	 * Geometric part of magnetic vector potential computation
	 * for circular wire loop at rho'=1, z'=0 (normalized coordinates).
	 * This routine selects special case routines
	 * to get the most accurate formulation for given evaluation coordinates.
	 *
	 * @param rhoP normalized radial evaluation position
	 * @param zP normalized vertical evaluation position
	 * @return A_phi: toroidal component of magnetic vector potential: geometric part (no mu0*I/pi factor included)
	 */
	public static double circularWireLoop_A_phi(double rhoP, double zP) {
		if (rhoP == 0.0) {
			return 0.0;
		} else if (zP >= 1.0 || rhoP < 0.5 || rhoP > 2.0) {
			return A_phi_1(rhoP, zP);
		} else if (rhoP != 1.0) {
			return A_phi_6(rhoP, zP);
		} else {
			// rhoP == 1, zP < 1
			return A_phi_5(rhoP, zP);
		}
	}

	/**
	 * Geometric part of radial magnetic field computation
	 * for circular wire loop at rho'=1, z'=0 (normalized coordinates).
	 * This routine selects special case routines
	 * to get the most accurate formulation for given evaluation coordinates.
	 *
	 * @param rhoP normalized radial evaluation position
	 * @param zP normalized vertical evaluation position
	 * @return B_rho: radial component of magnetic field: geometric part (no mu0*I/(pi*a) factor included)
	 */
	public static double circularWireLoop_B_rho(double rhoP, double zP) {
		if (rhoP == 0.0 || zP == 0.0) {
			return 0.0;
		} else if (zP >= 1.0 || rhoP < 0.5 || rhoP > 2.0) {
			return B_rho_3(rhoP, zP);
		} else if (rhoP != 1.0) {
			return B_rho_1(rhoP, zP);
		} else {
			return B_rho_4(rhoP, zP);
		}
	}

	/**
	 * Geometric part of vertical magnetic field computation
	 * for circular wire loop at rho'=1, z'=0 (normalized coordinates).
	 * This routine selects special case routines
	 * to get the most accurate formulation for given evaluation coordinates.
	 *
	 * @param rhoP normalized radial evaluation position
	 * @param zP normalized vertical evaluation position
	 * @return B_z: vertical component of magnetic field: geometric part (no mu0*I/(pi*a) factor included)
	 */
	public static double circularWireLoop_B_z(double rhoP, double zP) {
		if (rhoP < 0.5 || (rhoP < 2 && zP >= 1)) {
			return B_z_1(rhoP, zP);
		} else if (rhoP >= 2) {
			return B_z_2(rhoP, zP);
		} else if (rhoP == 1.0) {
			return B_z_4(rhoP, zP);
		} else if (zP != 0.0){
			return B_z_5(rhoP, zP);
		} else {
			return B_z_6(rhoP, zP);
		}
	}
















	static double A_z_1(double rhoP, double zP) {

		double Ri = Math.hypot(rhoP, zP);
		double Rf = Math.hypot(rhoP, 1.0 - zP);

		return FastMath.atanh(1.0 / (Ri + Rf));
	}






	/**
	 * combined solution for rhoP=0
	 * @param rhoP
	 * @param zP
	 * @return
	 */
	static double A_z_along_rhoP_0(double rhoP, double zP) {

		if (zP < -1 || zP > 2) {
			return A_z_2(rhoP, zP);
		} else {
			return A_z_2b(rhoP, zP);
		}
	}

	/**
	 * special case for rho'=0
	 * @param rhoP
	 * @param zP
	 * @return
	 */
	static double A_z_2(double rhoP, double zP) {
		return FastMath.atanh(1.0 / (Math.abs(zP) + Math.abs(1-zP)));
	}

	/**
	 * special case for rho'=0; near-field (excluding wire)
	 * @param rhoP
	 * @param zP
	 * @return
	 */
	static double A_z_2b(double rhoP, double zP) {
		return Math.signum(zP) * Math.log(Math.abs(zP/(1-zP)))/2;
	}





	/**
	 * combined solution for zP=0 or zP = 1
	 * @param rhoP
	 * @param zP
	 * @return
	 */
	static double A_z_along_zP_0_or_1(double rhoP, double zP) {

		if (rhoP > 1.0) {
			return A_z_3(rhoP, 0.0);
		} else {
			// rhoP <= 1
			return A_z_3b(rhoP, 0.0);
		}
	}

	/**
	 * special case for z'=0
	 * @param rhoP
	 * @param zP
	 * @return
	 */
	static double A_z_3(double rhoP, double zP) {
		return FastMath.atanh(1.0 / (rhoP + Math.sqrt(rhoP*rhoP + 1)));
	}

	/**
	 * special case for z'=0
	 * @param rhoP
	 * @param zP
	 * @return
	 */
	static double A_z_3b(double rhoP, double zP) {
		// a little bit more robust --> around rho'=1 +/- one test point we have 15 digits
		double cat = 1/Math.sqrt(rhoP*rhoP + 1);
		double sat = Math.sin(Math.atan(rhoP)/2);
		double num = rhoP*cat + 1 + cat;
		double den = rhoP*cat + 2*sat*sat;
		return Math.log(num/den)/2;
	}








	// (1): rho' < 1e-15, |z'|>=1
	static double A_z_6a(double rhoP, double zP) {

		double ang = Math.atan2(rhoP, zP);
		double s = Math.sin(ang/2);
		double c = Math.cos(ang);

		// nutritious zero: R_i - 1 == (R_i - z') + (z' - 1)
		double Ri_zP = zP * 2 * s*s / c; // R_i - z'

		double zpM1 = zP - 1.0;
		double Rf = Math.sqrt(rhoP * rhoP + zpM1*zpM1);

		double n = Ri_zP + Rf + zpM1;

		return (Math.log(2 + n) - Math.log(n)) / 2;
	}


	static double A_z_6b(double rhoP, double zP) {

		double alpha = Math.atan2(rhoP, zP);
		double sinAlphaHalf = Math.sin(alpha/2);
		double cosAlpha = Math.cos(alpha);

		// nutritious zero: R_i - 1 == (R_i - z') + (z' - 1)
		double Ri_zP = 2 * zP * sinAlphaHalf*sinAlphaHalf / cosAlpha; // R_i - z'

		double omz = 1.0 - zP;
		double beta = Math.atan2(rhoP, omz);
		double sinBetaHalf = Math.sin(beta/2);
		double cosBeta = Math.cos(beta);

		double Rf_p_zM1 = 2 * omz * sinBetaHalf*sinBetaHalf / cosBeta; // R_f - 1 + z'

		double n = Ri_zP + Rf_p_zM1;

		return (Math.log(2 + n) - Math.log(n)) / 2;
	}

	static double A_z_6c(double rhoP, double zP) {

		double alpha = Math.atan2(rhoP, 1-zP);
		double sinAlphaHalf = Math.sin(alpha/2);

		double R_i = Math.sqrt(rhoP*rhoP + zP*zP);
		double R_f = Math.sqrt(rhoP*rhoP + (1.0-zP)*(1.0-zP));

		double Rf_m_1 = 2.0 * R_f * sinAlphaHalf*sinAlphaHalf - zP;

		double n = R_i + Rf_m_1;
		double omEps = n/(n + 1);
		double opEps = 2 - omEps;
		return (Math.log(opEps) - Math.log(omEps))/2;
	}

	/**
	 * special case for rho'=0
	 * @param rhoP
	 * @param zP
	 * @return
	 */
	static double B_phi_2(double rhoP, double zP) {
		// works everywhere, although only derived for zP < 0 ???
		double zPM1 = 1 - zP;
		return (1 / (zPM1 * zPM1) - 1 / (zP * zPM1)) / 4;

	}

	/**
	 * special case for zP=0 or zP=1
	 * @param rhoP
	 * @param zP
	 * @return
	 */
	static double B_phi_3(double rhoP, double zP) {
		double zPM1 = 1 - zP;
		return 1 / (2 * rhoP * Math.hypot(rhoP, zPM1));
	}


	// near-field:
	// zp>=1 or zp<0, all rhoP
	// zP from 0 to 1/2, rhoP from 1e-30 at zp=0 to rhoP=1 at zP=1/2
	// zP from 1/2 to 1, rhoP from 1e-15 at zp=1 to rhoP=1 at zP=1/2
	static double B_phi_4(double rhoP, double zP) {

		double Ri = Math.hypot(rhoP, zP);
		double Rf = Math.hypot(rhoP, 1 - zP);

		double t1 = Ri/Rf + 1;

		double den = rhoP*rhoP + zP*(zP - 1) + Ri*Rf;

		return t1 / (2 * den);
	}

	// near-field: zP approx 1
	static double B_phi_5(double rhoP, double zP) {

		double Ri = Math.hypot(rhoP, zP);
		double Rf = Math.hypot(rhoP, 1 - zP);

		double t1 = Ri/Rf + 1;

		double beta = Math.atan2(rhoP, 1-zP);
		double cosBeta = Math.cos(beta);
		double sinBetaHalf = Math.sin(beta/2.0);

		double gamma = Math.atan2(rhoP, zP);
		double cosGamma = Math.cos(gamma);
		double sinGammaHalf = Math.sin(gamma/2.0);

		// (a*b - 1)
		double abm1 = 2.0/cosGamma * (sinBetaHalf*sinBetaHalf/cosBeta + sinGammaHalf*sinGammaHalf);

		// R_i*R_f - zP*(1-zP) == zP*(1-zP) * (a*b - 1)
		double den = 2.0 * (rhoP*rhoP + zP*(1.0 - zP) * abm1);

		return t1 / den;
	}

	static double A_phi_1(double rhoP, double zP) {
		// complementary modulus k_c
		double kCSqNum = zP*zP + 1 + rhoP * rhoP - 2.0 * rhoP;
		double kCSqDen = zP*zP + 1 + rhoP * rhoP + 2.0 * rhoP;

		double sqrt_kCSqNum = Math.sqrt(kCSqNum);
		double sqrt_kCSqDen = Math.sqrt(kCSqDen);
		double kC = sqrt_kCSqNum/sqrt_kCSqDen;

		double kSq = 4.0*rhoP / kCSqDen;

		double celPrefactor = 1.0/sqrt_kCSqDen;

		// Walstrom
		double arg1 = 2*Math.sqrt(kC)/(1+kC);
		double a2d =  1+kC;
		double arg2 = 2/(a2d*a2d*a2d);
		double C = CompleteEllipticIntegral.cel(arg1, 1, 0, arg2);

		return celPrefactor * kSq * C;
	}


	/**
	 * special case for rho'=1 and z' close to 0
	 * @param rhoP
	 * @param zP
	 * @return
	 */
	static double A_phi_5(double rhoP, double zP) {
		double kc = Math.sqrt(4.0 + zP * zP) / zP;
		return 1.0 / zP * CompleteEllipticIntegral.cel(kc, 1, 1, -1);
	}


	/**
	 *
	 * @param rhoP
	 * @param zP
	 * @return
	 */
	static double A_phi_6(double rhoP, double zP) {
		double n = zP/(rhoP-1);
		double m = 1+2/(rhoP-1);
		double den = n*n + m*m;
		double num = n*n + 1;
		double kcsq = num/den;

		double prefac = 1/(Math.abs(rhoP-1) * Math.sqrt(den));
		double celPart = CompleteEllipticIntegral.cel(Math.sqrt(kcsq), 1, -1, 1);
		return prefac * celPart;
	}


	static double B_rho_1(double rhoP, double zP) {
		double n = zP/(rhoP-1);
		double m = 1 + 2/(rhoP-1);
		double den = n*n + m*m;
		double num = n*n + 1;
		double kCSq = num/den;

		double kC = Math.sqrt(kCSq);

		double D = CompleteEllipticIntegral.cel(kC, 1, 0, 1);

		double arg1 = 2*Math.sqrt(kC)/(1+kC);
		double a2d =  1+kC;
		double arg2 = 2/(a2d*a2d*a2d);
		double C = CompleteEllipticIntegral.cel(arg1, 1, 0, arg2);

		double rd = rhoP - 1;
		double rd2 = rd*rd;

		double fac1 = 4*rhoP/(rd2*rd2) * Math.abs(n);

		return fac1 * (D-C) / (den*Math.sqrt(den) * num);
	}


	static double B_rho_3(double rhoP, double zP) {


		double sqrt_kCSqNum = Math.hypot(zP, 1 - rhoP);
		double sqrt_kCSqDen = Math.hypot(zP, 1 + rhoP);
		double k = 2 * Math.sqrt(rhoP) / sqrt_kCSqDen;
		double kSq = k * k;
		double kCSq = 1.0 - kSq;
		double kC = Math.sqrt(kCSq);

		double D = CompleteEllipticIntegral.cel(kC, 1, 0, 1);

		double arg1 = 2*Math.sqrt(kC)/(1+kC);
		double a2d =  1+kC;
		double arg2 = 2/(a2d*a2d*a2d);
		double C = CompleteEllipticIntegral.cel(arg1, 1, 0, arg2);

		double prefac = 4*rhoP*zP / (sqrt_kCSqDen*sqrt_kCSqDen*sqrt_kCSqDen * sqrt_kCSqNum*sqrt_kCSqNum);

		return prefac * (D - C);
	}

	/** special case for rhoP=1, zP --> 0 */
	static double B_rho_4(double rhoP, double zP) {

		double zPSq = zP*zP;
		double pfd = 1 + 4/zPSq;
		double kCSq = 1/pfd;
		double kC = Math.sqrt(kCSq);

		double E = CompleteEllipticIntegral.cel(kC, 1, 1, kCSq);
		double K = CompleteEllipticIntegral.cel(kC, 1, 1, 1);

		return 1/(2 * Math.sqrt(pfd)) * (E/pfd * (1 + (6 + 8/zPSq)/zPSq) - K );
	}

	static double B_z_1(double rhoP, double zP) {

		double sqrt_kCSqNum = Math.hypot(zP, 1 - rhoP);
		double kCSqNum = sqrt_kCSqNum * sqrt_kCSqNum;

		double sqrt_kCSqDen = Math.hypot(zP, 1 + rhoP);

		double k = 2 * Math.sqrt(rhoP) / sqrt_kCSqDen;
		double kSq = k * k;
		double kCSq = 1.0 - kSq;
		double kC = Math.sqrt(kCSq);

		double E = CompleteEllipticIntegral.cel(kC, 1, 1, kCSq);
		double K = CompleteEllipticIntegral.cel(kC, 1, 1, 1);
		double D = CompleteEllipticIntegral.cel(kC, 1, 0, 1);

		double comb = (E - 2*K + 2*D);

		double prefac = 1.0 / (sqrt_kCSqDen * kCSqNum);

		return prefac * ( E + rhoP * comb );
	}

	static double B_z_2(double rhoP, double zP) {
		// large rho'

		double sqrt_kCSqDen = Math.hypot(zP, 1 + rhoP);
		double k = 2 * Math.sqrt(rhoP) / sqrt_kCSqDen;
		double kSq = k * k;
		double kCSq = 1.0 - kSq;
		double kC = Math.sqrt(kCSq);

		double E = CompleteEllipticIntegral.cel(kC, 1, 1, kCSq);
		double D = CompleteEllipticIntegral.cel(kC, 1, 0, 1);

		double zp2_1 = zP*zP + 1;
		double r6 = zp2_1 * zp2_1 * zp2_1;
		double r5 = r6/rhoP - 2*zp2_1*zp2_1;
		double r4 = r5/rhoP + 3*zp2_1*zp2_1 - 4*zp2_1;
		double r3 = r4/rhoP - 4*zP*zP + 4;
		double r2 = r3/rhoP + 3*zP*zP - 1;
		double r1 = r2/rhoP - 2;
		double r0 = r1/rhoP + 1;

		// use C-D for (2D-E)/kSq: this finally works without expansion !!!
		double arg1 = 2*Math.sqrt(kC)/(1+kC);
		double a2d =  1+kC;
		double arg2 = 2/(a2d*a2d*a2d);
		double C = CompleteEllipticIntegral.cel(arg1, 1, 0, arg2);

		double cdScale = 1 + (2 + (zP*zP+1)/rhoP)/rhoP;

		return 1.0 /(Math.sqrt(r0) * rhoP*rhoP*rhoP) * (E + 4*(C-D)/cdScale);
	}

	static double B_z_4(double rhoP, double zP) {
		// special case for rhoP=1, zp->0

		double kCSq = zP*zP/(4 + zP*zP);
		double kC = Math.sqrt(kCSq);

		double f = zP*zP + 4;
		double prefac = 1/(f*Math.sqrt(f));

		return prefac * CompleteEllipticIntegral.cel(kC, kCSq, 2, 0);
	}

	static double B_z_5(double rhoP, double zP) {
		// special case for near-field: rhoP->1, zP->0; but not rhoP=1 or zP=0

		double rp1 = rhoP - 1;

		double n = zP/rp1;
		double m = 1 + 2/rp1;
		double den = n*n + m*m;
		double num = n*n + 1;
		double kCSq = num/den;

		double prefac = Math.abs(n)/(rp1*rp1 * den*Math.sqrt(den));

		double ca1 = (1+rhoP)/zP;
		double ca2 = (1-rhoP)/zP;
		double cp = CompleteEllipticIntegral.cel(Math.sqrt(kCSq), kCSq, ca1, ca2);

		return prefac*cp;
	}

	static double B_z_6(double rhoP, double zP) {

		double rp = rhoP + 1;

		double kC = (1-rhoP)/rp;
		double kCSq = kC*kC;

		double ca1 = 1+rhoP;
		double ca2 = 1-rhoP;
		double cp = CompleteEllipticIntegral.cel(Math.sqrt(kCSq), kCSq, ca1, ca2);

		return 1/Math.abs(rp*rp*rp) * cp;
	}
}
