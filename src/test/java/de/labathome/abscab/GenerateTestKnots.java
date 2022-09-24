package de.labathome.abscab;

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

	/**
	 * Generate the set of test points used to test the straight wire segment
	 * methods.
	 */
	public static void generateTestKnotsStraightWireSegment() {

		/** assemble list of test knots in radial direction */
		List<Double> testKnotsRpLst = new LinkedList<>();

		// 0
		testKnotsRpLst.add(0.0);

		// 1e-30 ... 1e30
		for (int exponent = -30; exponent <= 30; ++exponent) {
			testKnotsRpLst.add(Math.pow(10.0, exponent));
		}

		double[] testKnotsRp = testKnotsRpLst.stream().mapToDouble(d -> d).toArray();

		System.out.printf("number of test knots in  radial  direction: %d\n", testKnotsRp.length);

//		System.out.println("test knots in radial direction:");
//		for (int i = 0; i < testKnotsRp.length; ++i) {
//			System.out.printf("% .17e\n", testKnotsRp[i]);
//		}

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

		double[] testKnotsZp = testKnotsZpLst.stream().mapToDouble(d -> d).toArray();

		System.out.printf("number of test knots in vertical direction: %d\n", testKnotsZp.length);

//		System.out.println("test knots in vertical direction:");
//		for (int i = 0; i < testKnotsZp.length; ++i) {
//			System.out.printf("% .17e\n", testKnotsZp[i]);
//		}

		UtilsTestABSCAB.dumpToFile(testKnotsRp, "src/test/resources/testKnotsRpStraightWireSegment.dat");
		UtilsTestABSCAB.dumpToFile(testKnotsZp, "src/test/resources/testKnotsZpStraightWireSegment.dat");

		/** assemble test points from knots; exclude locations on wire segment */
		List<Double> testPointsRpLst = new LinkedList<>();
		List<Double> testPointsZpLst = new LinkedList<>();

		List<Integer> idxRpLst = new LinkedList<>();
		List<Integer> idxZpLst = new LinkedList<>();

		int numCases = 0;
		for (int idxZ = 0; idxZ < testKnotsZp.length; ++idxZ) {
			double zp = testKnotsZp[idxZ];

			for (int idxR = 0; idxR < testKnotsRp.length; ++idxR) {
				double rp = testKnotsRp[idxR];

				// exclude locations on the wire segment
				if (rp != 0.0 || zp < 0.0 || zp > 1.0) {
					testPointsRpLst.add(rp);
					testPointsZpLst.add(zp);

					idxRpLst.add(idxR);
					idxZpLst.add(idxZ);

					numCases++;
				}
			}
		}

		System.out.printf("number of test cases: %d\n", numCases);

		double[] testPointsRp = testPointsRpLst.stream().mapToDouble(d -> d).toArray();
		double[] testPointsZp = testPointsZpLst.stream().mapToDouble(d -> d).toArray();

		int[] idxRp = idxRpLst.stream().mapToInt(i -> i).toArray();
		int[] idxZp = idxZpLst.stream().mapToInt(i -> i).toArray();

		UtilsTestABSCAB.dumpToFile(testPointsRp, "src/test/resources/testPointsRpStraightWireSegment.dat");
		UtilsTestABSCAB.dumpToFile(testPointsZp, "src/test/resources/testPointsZpStraightWireSegment.dat");

		UtilsTestABSCAB.dumpToFile(idxRp, "src/test/resources/idxRpStraightWireSegment.dat");
		UtilsTestABSCAB.dumpToFile(idxZp, "src/test/resources/idxZpStraightWireSegment.dat");

		UtilsTestABSCAB.dumpTestPoints(testPointsRp, testPointsZp, "src/test/resources/testPointsStraightWireSegment.dat");
	}

	/**
	 * Generate the set of test points used to test the circular wire loop methods.
	 */
	public static void generateTestKnotsCircularWireLoop() {

		/** assemble list of test knots in radial direction */
		List<Double> testKnotsRpLst = new LinkedList<>();

		// 0
		testKnotsRpLst.add(0.0);

		// 1e-30 ... 1e-1
		for (int exponent = -30; exponent <= -1; ++exponent) {
			testKnotsRpLst.add(Math.pow(10.0, exponent));
		}

		// 1/2
		testKnotsRpLst.add(0.5);

		// 1 - 1e-1 ... 1 - 1e-15
		for (int exponent = -1; exponent >= -15; --exponent) {
			testKnotsRpLst.add(1.0 - Math.pow(10.0, exponent));
		}

		// 1 - eps/2
		testKnotsRpLst.add(Math.nextDown(1.0));

		// 1
		testKnotsRpLst.add(1.0);

		// 1 + eps
		testKnotsRpLst.add(Math.nextUp(1.0));

		// 1 + 1e-15 ... 1 + 1e-1
		for (int exponent = -15; exponent <= -1; ++exponent) {
			testKnotsRpLst.add(1.0 + Math.pow(10.0, exponent));
		}

		// 2
		testKnotsRpLst.add(2.0);

		// 1e1 ... 1e30
		for (int exponent = 1; exponent <= 30; ++exponent) {
			testKnotsRpLst.add(Math.pow(10.0, exponent));
		}

		double[] testKnotsRp = testKnotsRpLst.stream().mapToDouble(d -> d).toArray();

		System.out.printf("number of test knots in  radial  direction: %d\n", testKnotsRp.length);

//		System.out.println("test knots in radial direction:");
//		for (int i = 0; i < testKnotsRp.length; ++i) {
//			System.out.printf("% .17e\n", testKnotsRp[i]);
//		}

		/** assemble list of test knots in vertical direction */
		List<Double> testKnotsZpLst = new LinkedList<>();

		// 0
		testKnotsZpLst.add(0.0);

		// 1e-30 ... 1e30
		for (int exponent = -30; exponent <= 30; ++exponent) {
			testKnotsZpLst.add(Math.pow(10.0, exponent));
		}

		double[] testKnotsZp = testKnotsZpLst.stream().mapToDouble(d -> d).toArray();

		System.out.printf("number of test knots in vertical direction: %d\n", testKnotsZp.length);

//		System.out.println("test knots in vertical direction:");
//		for (int i = 0; i < testKnotsZp.length; ++i) {
//			System.out.printf("% .17e\n", testKnotsZp[i]);
//		}

		UtilsTestABSCAB.dumpToFile(testKnotsRp, "src/test/resources/testKnotsRpCircularWireLoop.dat");
		UtilsTestABSCAB.dumpToFile(testKnotsZp, "src/test/resources/testKnotsZpCircularWireLoop.dat");

		/** assemble test points from knots; exclude locations on wire segment */
		List<Double> testPointsRpLst = new LinkedList<>();
		List<Double> testPointsZpLst = new LinkedList<>();

		List<Integer> idxRpLst = new LinkedList<>();
		List<Integer> idxZpLst = new LinkedList<>();

		int numCases = 0;
		for (int idxZ = 0; idxZ < testKnotsZp.length; ++idxZ) {
			double zp = testKnotsZp[idxZ];

			for (int idxR = 0; idxR < testKnotsRp.length; ++idxR) {
				double rp = testKnotsRp[idxR];

				// exclude location on the wire loop
				if (rp != 1.0 || zp != 0.0) {
					testPointsRpLst.add(rp);
					testPointsZpLst.add(zp);

					idxRpLst.add(idxR);
					idxZpLst.add(idxZ);

					numCases++;
				}
			}
		}

		System.out.printf("number of test cases: %d\n", numCases);

		double[] testPointsRp = testPointsRpLst.stream().mapToDouble(d -> d).toArray();
		double[] testPointsZp = testPointsZpLst.stream().mapToDouble(d -> d).toArray();

		int[] idxRp = idxRpLst.stream().mapToInt(i -> i).toArray();
		int[] idxZp = idxZpLst.stream().mapToInt(i -> i).toArray();

		UtilsTestABSCAB.dumpToFile(testPointsRp, "src/test/resources/testPointsRpCircularWireLoop.dat");
		UtilsTestABSCAB.dumpToFile(testPointsZp, "src/test/resources/testPointsZpCircularWireLoop.dat");

		UtilsTestABSCAB.dumpToFile(idxRp, "src/test/resources/idxRpCircularWireLoop.dat");
		UtilsTestABSCAB.dumpToFile(idxZp, "src/test/resources/idxZpCircularWireLoop.dat");

		UtilsTestABSCAB.dumpTestPoints(testPointsRp, testPointsZp, "src/test/resources/testPointsCircularWireLoop.dat");
	}
}
