#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <math.h>

/**
 * Check if two values are approximately equal within a prescribed tolerance.
 * For values much smaller than 1, this is similar to a comparison of the
 * absolute values. For values much greater than 1, this is similar to a
 * comparison of the relative values.
 *
 * This method is described in Gill, Murray & Wright, "Practical Optimization" (1984).
 *
 * @param expected  expected result
 * @param actual    actual result
 * @param tolerance relative or absolute tolerance on the mismatch between the
 *                  expected and the actual values
 * @return 0 if the values match within the prescribed tolerance; 1 otherwise
 */
int assertRelAbsEquals(double expected, double actual, double tolerance) {
	double relAbsError = fabs(actual - expected) / (1.0 + fabs(expected));
	if (relAbsError > tolerance) {
		printf("expected %g, actual %g (rel/abs error %g, tolerance %g)\n",
				expected, actual, relAbsError, tolerance);
		return 1;
	}
	return 0;
}

#endif // UTIL_H
