package de.labathome;

import java.util.LinkedList;
import java.util.List;
import java.util.Locale;

public class GenerateTestKnots {

	/** machine precision: approximately 2.22e-16 for 64-bit double precision */
	public static final double EPS = Math.ulp(1.0);

	public static void main(String[] args) {
		Locale.setDefault(Locale.ENGLISH);

		System.out.printf("machine precision: %.17e\n", EPS);

		generateTestKnotsStraightWireSegment();
		generateTestKnotsCircularWireLoop();
	}

	public static void generateTestKnotsStraightWireSegment() {

		/** assemble list of test knots in radial direction */
		List<Double> testKnotsRpLst = new LinkedList<>();

		// 0
		testKnotsRpLst.add(0.0);

		// 1e-30 ... 1e30
		for (int exponent = -30; exponent <= 30; ++exponent) {
			testKnotsRpLst.add(Math.pow(10.0, exponent));
		}

		final double[] testKnotsRp = testKnotsRpLst.stream().mapToDouble(d -> d).toArray();

		System.out.printf("number of test knots in  radial  direction: %d\n", testKnotsRp.length);

		/** assemble list of test knots in vertical direction */
		List<Double> testKnotsZpLst = new LinkedList<>();

		// -1e30 ... -1e-30
		for (int exponent = 30; exponent >= -30; --exponent) {
			testKnotsZpLst.add(-Math.pow(10.0, exponent));
		}

		// 0
		testKnotsZpLst.add(0.0);

		// 1e-30 ... 1e-1
		for (int exponent = -30; exponent <= -1; ++exponent) {
			testKnotsZpLst.add(Math.pow(10.0, exponent));
		}

		// 1/2
		testKnotsZpLst.add(0.5);

		// 1 - 1e-1 ... 1 - 1e-15
		for (int exponent = -1; exponent >= -15; --exponent) {
			testKnotsZpLst.add(1.0 - Math.pow(10.0, exponent));
		}

		// 1 - eps/2 (next lower double precision number from 1)
		testKnotsZpLst.add(Math.nextDown(1.0));

		// 1
		testKnotsZpLst.add(1.0);

		// 1 + eps (next higher double precision number from 1)
		testKnotsZpLst.add(Math.nextUp(1.0));

		// 1 + 1e-15 ... 1 + 1e-1
		for (int exponent = -15; exponent <= -1; ++exponent) {
			testKnotsZpLst.add(1.0 + Math.pow(10.0, exponent));
		}

		// 2
		testKnotsZpLst.add(2.0);

		// 1e1 ... 1e30
		for (int exponent = 1; exponent <= 30; ++exponent) {
			testKnotsZpLst.add(Math.pow(10.0, exponent));
		}

		final double[] testKnotsZp = testKnotsZpLst.stream().mapToDouble(d -> d).toArray();

		System.out.printf("number of test knots in vertical direction: %d\n", testKnotsZp.length);

		/** assemble test points from knots; exclude locations on wire segment */

	}

	public static void generateTestKnotsCircularWireLoop() {

	}
}
