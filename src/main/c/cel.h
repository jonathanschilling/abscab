#ifndef CEL_H
#define CEL_H

#include <float.h>
#include <limits.h>
#include <math.h>

// C99 does not define M_PI
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

/** half of pi */
const double PI_2 = M_PI / 2.0;

/** sqrt of machine epsilon */
const double SQRT_EPS = sqrt(DBL_EPSILON);

/**
 * Compute the complete elliptic integral introduced in
 * "Numerical Calculation of Elliptic Integrals and Elliptic Functions. III"
 * by R. Bulirsch in "Numerische Mathematik" 13, 305-315 (1969):
 * cel(k_c, p, a, b) =
 * \int_0^{\pi/2} \frac{a \cos^2{\varphi} + b \sin^2{\varphi}}
 *                     {  \cos^2{\varphi} + p \sin^2{\varphi}}
 *                \frac{\mathrm{d}\varphi}
 *                     {\sqrt{\cos^2{\varphi} + k_c^2 \sin^2{\varphi}}}
 * @param k_c parameter k_c of cel(); absolute value must not be 0
 * @param p   parameter p of cel()
 * @param a   parameter a of cel()
 * @param b   parameter b of cel()
 * @return the value of cel(k_c, p, a, b)
 */
double cel(double k_c, double p, double a, double b) {
	if (k_c == 0.0) {
		if (b != 0.0) {
			// when k_c is zero and b != 0, cel diverges (?)
			return HUGE_VAL;
		} else {
			k_c = SQRT_EPS*SQRT_EPS;
		}
	} else {
		k_c = fabs(k_c);
	}

	double m = 1.0; // \mu
	double e = k_c; // \nu * \mu
	// In the iterations, \nu_i is stored in k_c.

	double f, g;

	// initialization
	if (p > 0.0) {
		p = sqrt(p);
		b = b / p;
	} else { // p <= 0
		double q;

		f = k_c*k_c;        // f = kc^2 (re-used here; later f = a_i)
		q = 1.0 - f;        // 1 - kc^2
		g = 1.0 - p;
		f -= p;             // kc^2 - p
		q *= b-a*p;         // (1 - kc^2)*(b-a*p)
		p = sqrt(f/g); // sqrt((kc^2 - p)/(1-p))              --> p0 done here
		a = (a-b)/g;        // (a-b)/(1-p)                         --> a0 done here
		b = -q/(g*g*p)+a*p; // -(1 - kc^2)*(b-a*p)/( (1-p)^2 * p ) --> b0 done here
	}

	// iteration until convergence
	while (1) {
		f = a;    // f = a_i
		a += b/p; // a_{i+1} <-- a_i + b_i/p_i
		g = e/p;  // g = (\nu_i*\mu_i) / p_i
		b += f*g; // b_i + a_i * (\nu_i*\mu_i) /p_i
		b += b;   // b_{i+1} <-- 2*( b_i + a_i * (\nu_i*\mu_i) /p_i )
		p += g;   // p_{i+1} <-- p_i + (\nu_i*\mu_i) / p_i
		g = m;    // g = mu_i
		m += k_c; // mu_{i+1} <-- mu_i + \nu_i
		if (fabs(g - k_c) > g * SQRT_EPS) { // |\mu_i - \nu_i| > \mu_i*EPS == |1 - \nu_i / \mu_i| > EPS, but more robust!
			// not converged yet...
			k_c = sqrt(e); // k_c = sqrt(\nu_i * \mu_i)
			k_c += k_c;         // \nu_{i+1} <-- 2*sqrt(\nu_i * \mu_i)
			e = k_c*m; // (\nu * \mu)_{i+1} <-- \nu_{i+1} * mu_{i+1} (also update product explicitly)
		} else {
			break;
		}
	}

	// final approximation:
	// \pi/2 * (a * \mu + b)/(\mu * (\mu + p))
	return PI_2 * (a*m+b) / (m*(m+p));
}

#endif // CEL_H
