#ifndef COMPSUM_H
#define COMPSUM_H

#include <math.h>

/**
 * Add a single contribution to the sum.
 * The compensated sum is obtained by summing the final values of s, cs and ccs
 * after this method has been called for all contributions.
 *
 * @param contribution contribution to add to the sum
 * @param compSum[3]: {s, cs, ccs}: target for output
 */
void compAdd(double contribution, double *compSum) {
	double s   = compSum[0];
	double cs  = compSum[1];

	double t = s + contribution;
	double c;
	if (fabs(s) >= fabs(contribution)) {
		c = (s - t) + contribution;
	} else {
		c = (contribution - t) + s;
	}
	compSum[0] = t;

	double t2 = cs + c;
	double cc;
	if (fabs(cs) >= fabs(c)) {
		cc = (cs - t2) + c;
	} else {
		cc = (c - t2) + cs;
	}
	compSum[1] = t2;
	compSum[2] += cc;
}

#endif // COMPSUM_H
