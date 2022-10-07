package de.labathome.abscab;

import java.util.Locale;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

public class TestCEL {

	/**
	 * Check if two values are approximately equal within a prescribed tolerance.
	 * For values much smaller than 1, this is similar to a comparison of the
	 * absolute values. For values much greater than 1, this is similar to a
	 * comparison of the relative values.
	 *
	 * This method is described in Gill, Murray & Wright, "Practical Optimization"
	 * (1984).
	 *
	 * @param expected  expected result
	 * @param actual    actual result
	 * @param tolerance relative or absolute tolerance on the mismatch between the
	 *                  expected and the actual values
	 * @return true if the values match within the prescribed tolerance; false
	 *         otherwise
	 */
	public static boolean assertRelAbsEquals(double expected, double actual, double tolerance) {
		final double relAbsError = Math.abs(actual - expected) / (1.0 + Math.abs(expected));
		if (relAbsError > tolerance) {
			Assertions.fail(String.format(Locale.ENGLISH, "expected %g, actual %g (rel/abs error %g, tolerance %g)",
					expected, actual, relAbsError, tolerance));
			return false;
		}
		return true;
	}

	/**
	 * test case for cel() implementation as described in section 4.2 of the 1969 Bulirsch article
	 */
	@Test
	public void testCel() {
		final double tolerance = 1.0e-15;

		final double k_c = 0.1;
		final double p1 =  4.1;
		final double p2 = -4.1;
		final double a = 1.2;
		final double b = 1.1;

		final double cel1 =  1.5464442694017956;
		final double cel2 = -6.7687378198360556e-1;

		final double c1 = CompleteEllipticIntegral.cel(k_c, p1, a, b);
		final double c2 = CompleteEllipticIntegral.cel(k_c, p2, a, b);

//		final double ra1 = Math.abs(cel1 - c1)/(1.0 + Math.abs(cel1));
//		final double ra2 = Math.abs(cel2 - c2)/(1.0 + Math.abs(cel2));
//		System.out.println("case 1: rel/abs deviation = " + ra1);
//		System.out.println("case 2: rel/abs deviation = " + ra2);

		assertRelAbsEquals(cel1, c1, tolerance);
		assertRelAbsEquals(cel2, c2, tolerance);
	}

	@Test
	public void testCompleteEllipticIntegrals() {
		final double tolerance = 1.0e-15;

		// analytical test cases
		double kAt0 = CompleteEllipticIntegral.ellipticK(0.0);
		double eAt0 = CompleteEllipticIntegral.ellipticE(0.0);
//		System.out.printf("  k=0.0  => K(k)=%g  E(k)=%g\n", kAt0, eAt0);
		assertRelAbsEquals(Math.PI/2.0, kAt0, tolerance);
		assertRelAbsEquals(Math.PI/2.0, eAt0, tolerance);

		double kAt1 = CompleteEllipticIntegral.ellipticK(1.0);
		double eAt1 = CompleteEllipticIntegral.ellipticE(1.0);
//		System.out.printf("  k=1.0  => K(k)=%g E(k)=%g\n", kAt1, eAt1);
		assertRelAbsEquals(Double.POSITIVE_INFINITY, kAt1, tolerance);
		assertRelAbsEquals(1.0,                      eAt1, tolerance);

		// from J. M. Hammersley, "Tables of Complete Elliptic Integrals",
		// J. Res. Nat. Bureau Std., Vol. 50, No. 1, Jan. 1953 (Research Paper 2386)
		double kAtHalf = CompleteEllipticIntegral.ellipticK(0.5);
		double eAtHalf = CompleteEllipticIntegral.ellipticE(0.5);
//		System.out.printf("  k=0.5  => K(k)=%g  E(k)=%g\n", kAtHalf, eAtHalf);
		assertRelAbsEquals(1.685750355, kAtHalf, 1.0e-9);
		assertRelAbsEquals(1.467462209, eAtHalf, 1.0e-9);

		double kClose1 = CompleteEllipticIntegral.ellipticK(1.0/1.01);
		double eClose1 = CompleteEllipticIntegral.ellipticE(1.0/1.01);
//		System.out.printf("1/k=1.01 => K(k)=%g  E(k)=%g\n", kClose1, eClose1);
		assertRelAbsEquals(3.361458120, kClose1, 1.0e-9);
		assertRelAbsEquals(1.028242731, eClose1, 1.0e-9);
	}
}
