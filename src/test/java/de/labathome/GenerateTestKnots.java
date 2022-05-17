package de.labathome;

import java.io.File;
import java.io.PrintWriter;
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

		/** assemble test points from knots; exclude locations on wire segment */
		List<Double> testPointsRpLst = new LinkedList<>();
		List<Double> testPointsZpLst = new LinkedList<>();

		int numCases = 0;
		for (int idxZ = 0; idxZ < testKnotsZp.length; ++idxZ) {
			double zp = testKnotsZp[idxZ];

			for (int idxR = 0; idxR < testKnotsRp.length; ++idxR) {
				double rp = testKnotsRp[idxR];

				// exclude locations on the wire segment
				if (rp != 0.0 || zp < 0.0 || zp > 1.0) {
					testPointsRpLst.add(rp);
					testPointsZpLst.add(zp);

					numCases++;
				}
			}
		}

		System.out.printf("number of test cases: %d\n", numCases);

		double[] testPointsRp = testPointsRpLst.stream().mapToDouble(d -> d).toArray();
		double[] testPointsZp = testPointsZpLst.stream().mapToDouble(d -> d).toArray();

		dumpToFile(testPointsRp, "src/test/resources/testPointsRpStraightWireSegment.dat");
		dumpToFile(testPointsZp, "src/test/resources/testPointsZpStraightWireSegment.dat");

		dumpTestPoints(testPointsRp, testPointsZp, "src/test/resources/testPointsStraightWireSegment.dat");
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

		/** assemble test points from knots; exclude locations on wire segment */
		List<Double> testPointsRpLst = new LinkedList<>();
		List<Double> testPointsZpLst = new LinkedList<>();

		int numCases = 0;
		for (int idxZ = 0; idxZ < testKnotsZp.length; ++idxZ) {
			double zp = testKnotsZp[idxZ];

			for (int idxR = 0; idxR < testKnotsRp.length; ++idxR) {
				double rp = testKnotsRp[idxR];

				// exclude location on the wire loop
				if (rp != 1.0 || zp != 0.0) {
					testPointsRpLst.add(rp);
					testPointsZpLst.add(zp);

					numCases++;
				}
			}
		}

		System.out.printf("number of test cases: %d\n", numCases);

		double[] testPointsRp = testPointsRpLst.stream().mapToDouble(d -> d).toArray();
		double[] testPointsZp = testPointsZpLst.stream().mapToDouble(d -> d).toArray();

		dumpToFile(testPointsRp, "src/test/resources/testPointsRpCircularWireLoop.dat");
		dumpToFile(testPointsZp, "src/test/resources/testPointsZpCircularWireLoop.dat");

		dumpTestPoints(testPointsRp, testPointsZp, "src/test/resources/testPointsCircularWireLoop.dat");
	}

	/**
	 * Write a given vector to a text file; one line per element.
	 *
	 * @param arr      array to write to a file
	 * @param filename file into which to write the given array
	 */
	private static void dumpToFile(double[] arr, String filename) {
		File outFile = new File(filename);
		try (PrintWriter pw = new PrintWriter(outFile)) {
			for (int i = 0; i < arr.length; ++i) {
				pw.printf(Locale.ENGLISH, "%+.20e\n", arr[i]);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * Write the given set of test points, such that the implied values can be
	 * represented exactly in arbitrary-precision software.
	 *
	 * @param testPointsRp [numCases] set of rp test point coordinates
	 * @param testPointsZp [numCases] set of zp test point coordinates
	 * @param filename     file into which to write the given test points
	 */
	private static void dumpTestPoints(double[] testPointsRp, double[] testPointsZp, String filename) {
		int numCases = testPointsRp.length;

		//                          64   60   56   52   48   44   40   36   32   28   24   20   16   12   8    4
		final long signMask     = 0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000L;
		final long exponentMask = 0b0111_1111_1111_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000L;
		final long mantissaMask = 0b0000_0000_0000_1111_1111_1111_1111_1111_1111_1111_1111_1111_1111_1111_1111_1111L;

		File outFile = new File(filename);
		try (PrintWriter pw = new PrintWriter(outFile)) {
			pw.println("# rp: sign bit, exponent E, mantiassa M; zp: sign bit, exponent E, mantiassa M");

			for (int i = 0; i < numCases; ++i) {

				// format of an IEEE754 double precision number F:
				// F = (-1)^s * 2^{E - 1023} * (1 + M/2^{52})
				// where:
				// s = 0, 1            ( 1 bit )
				// E = 0, 1, ..., 2047 (11 bits)
				// M = 0, 1, ...       (52 bits)

				long rpBits = Double.doubleToRawLongBits(testPointsRp[i]);
				long zpBits = Double.doubleToRawLongBits(testPointsZp[i]);

				long signRp     = (rpBits &     signMask) >>> 63;
				long exponentRp = (rpBits & exponentMask) >>> 52;
				long mantissaRp = (rpBits & mantissaMask);

				long signZp     = (zpBits &     signMask) >>> 63;
				long exponentZp = (zpBits & exponentMask) >>> 52;
				long mantissaZp = (zpBits & mantissaMask);

				pw.printf(Locale.ENGLISH, "%1d %4d %16d %1d %4d %16d\n",
						signRp, exponentRp, mantissaRp,
						signZp, exponentZp, mantissaZp);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
