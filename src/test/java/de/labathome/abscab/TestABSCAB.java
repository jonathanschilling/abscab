package de.labathome.abscab;

import java.util.Locale;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

public class TestABSCAB {

	@Test
	public void testMagneticFieldInfiniteLineFilament() {
		double tolerance = 1.0e-15;

		// Demtroeder 2, Sec. 3.2.2 ("Magnetic field of a straight wire")
		// B(r) = mu_0 * I / (2 pi r)
		// Test this here with:
		// I = 123.0 A
		// r = 0.132 m
		// => B = 0.186 mT
		double current = 123.0;
		double r = 0.132;
		double bPhiRef = ABSCAB.MU_0 * current / (2.0 * Math.PI * r);
//		System.out.printf("ref bPhi = %.5e\n", bPhiRef);

		double[][] vertices = {
				{0.0, 0.0},
				{0.0, 0.0},
				{-1e6, 1e6}
		};

		double[][] evalPos = {
				{r},
				{0.0},
				{0.0}
		};

		// y component is B_phi
		double bPhi = ABSCAB.magneticFieldPolygonFilament(vertices, current, evalPos)[1][0];
//		System.out.printf("act bPhi = %.5e\n", bPhi);

		double relAbsErr = Math.abs(bPhi - bPhiRef) / (1.0 + Math.abs(bPhiRef));
//		System.out.printf("raErr = %.5e\n", relAbsErr);

		Assertions.assertTrue(relAbsErr < tolerance);
	}

	@Test
	public void testBPhiInfiniteLineFilament() {
		double tolerance = 1.0e-15;

		// Demtroeder 2, Sec. 3.2.2 ("Magnetic field of a straight wire")
		// B(r) = mu_0 * I / (2 pi r)
		// Test this here with:
		// I = 123.0 A
		// r = 0.132 m
		// => B = 0.186 mT
		double current = 123.0;
		double r = 0.132;
		double bPhiRef = ABSCAB.MU_0 * current / (2.0 * Math.PI * r);
//		System.out.printf("ref bPhi = %.5e\n", bPhiRef);

		// half the length of the wire segment
		double halfL = 1e6;
		double L = 2*halfL;
		double rhoP = r / L;
		double zP = halfL / L;
		double bPhi = ABSCAB.MU_0 * current / (4.0 * Math.PI * L) * ABSCAB.straightWireSegment_B_phi(rhoP, zP);
//		System.out.printf("act bPhi = %.5e\n", bPhi);

		double relAbsErr = Math.abs(bPhi - bPhiRef) / (1.0 + Math.abs(bPhiRef));
//		System.out.printf("raErr = %.5e\n", relAbsErr);

		Assertions.assertTrue(relAbsErr < tolerance);
	}

	@Test
	public void testMagneticFieldInsideLongCoil() {
		double tolerance = 1.0e-4;

		// Demtroeder 2, Sec. 3.2.3 ("Magnetic field of a long coil")
		// B_z = mu_0 * n * I
		// where n is the winding density: n = N / L
		// of a coil of N windings over a length L
		// Example (which is tested here):
		// n = 1e3 m^{-1}
		// I = 10 A
		// => B = 0.0126T
		double bZRef = 0.0126;

		double bZ = 0.0;

		int N = 50000; // windings
		double L = 50.0; // total length of coil in m
		double n = N/L;

		double current = 10.0; // A
		double radius = 1.0; // m

		for (int i=0; i<N; ++i) {

			// axial position of coil
			double z0 = -L/2.0 + (i + 0.5) / n;

			// compute magnetic field
			//double prefac = ABSCAB.MU_0 * current / (2.0*Math.PI * radius);
			double prefac = ABSCAB.MU_0 * current / (Math.PI * radius); // TODO: factor of 2 ???
			double bZContrib = prefac * ABSCAB.circularWireLoop_B_z(0.0, z0);

//			System.out.printf("coil %d at z0 = % .3e => contrib = %.3e\n", i, z0, bZContrib);

			bZ += bZContrib;
		}
//		System.out.printf("B_z = %.5e\n", bZ);

		double relAbsErr = Math.abs(bZ - bZRef) / (1.0 + Math.abs(bZRef));
//		System.out.printf("raErr = %.5e\n", relAbsErr);

		Assertions.assertTrue(relAbsErr < tolerance);
	}

	@Test
	public void testMagneticFieldAtCenterOfWireLoop() {
		double tolerance = 1.0e-13;

		// Demtroeder 2, Sec. 3.2.6b ("Magnetic field of a circular wire loop")
		// B_z = mu_0 * I / (2 a)
		// Test this here with:
		// I = 123.0 A (loop current)
		// a = 13.2 m (loop radius)
		// => B = 5.85 uT
		double current = 123.0;
		double a = 13.2;
		double bZRef = ABSCAB.MU_0 * current / (2.0 * a);
//		System.out.printf("ref bZ = %.5e\n", bZRef);

		double[] center = { 0.0, 0.0, 0.0 };
		double[] normal = { 0.0, 0.0, 1.0 };

		double[][] evalPos = {
				{0.0},
				{0.0},
				{0.0}
		};

		double bZ = ABSCAB.magneticFieldCircularFilament(center, normal, a, current, evalPos)[2][0];
//		System.out.printf("act bZ = %.5e\n", bZ);

		double relAbsErr = Math.abs(bZ - bZRef) / (1.0 + Math.abs(bZRef));
//		System.out.printf("raErr = %.5e\n", relAbsErr);

		Assertions.assertTrue(relAbsErr < tolerance);
	}

	@Test
	public void testStraightWireSegment_A_z() {
		final double tolerance = 1.0e-15;

		// load set of test points
		double[] testPointsRp = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/testPointsRpStraightWireSegment.dat")[0];
		double[] testPointsZp = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/testPointsZpStraightWireSegment.dat")[0];

		int numCases = testPointsRp.length;

		// load reference data
		double[] ref_A_z = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/StraightWireSegment_A_z_ref.dat")[0];

		// compute A_z at test points and compare against reference
		for (int i=0; i<numCases; ++i) {
			double rhoP = testPointsRp[i];
			double zP   = testPointsZp[i];

			double aZ = ABSCAB.straightWireSegment_A_z(rhoP, zP);

			if (ref_A_z[i] == 0.0) {
				// exact zero has to be reproduced exactly
				Assertions.assertEquals(0.0, aZ);
			} else {
				double relErr = Math.abs((aZ - ref_A_z[i])/ref_A_z[i]);
				if (relErr >= tolerance) {
					System.out.printf(Locale.ENGLISH, "case %4d (rhoP=%.20e zP=%.20e) => relErr = %.3e\n",
							i, rhoP, zP, relErr);
				}
				Assertions.assertTrue(relErr < tolerance);
			}
		}
	}

	@Test
	public void testStraightWireSegment_B_phi() {
		final double tolerance = 1.0e-15;

		// load set of test points
		double[] testPointsRp = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/testPointsRpStraightWireSegment.dat")[0];
		double[] testPointsZp = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/testPointsZpStraightWireSegment.dat")[0];

		int numCases = testPointsRp.length;

		// load reference data
		double[] ref_B_phi = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/StraightWireSegment_B_phi_ref.dat")[0];

		// compute B_phi at test points and compare against reference
		for (int i=0; i<numCases; ++i) {
			double rhoP = testPointsRp[i];
			double zP   = testPointsZp[i];

			double bPhi = ABSCAB.straightWireSegment_B_phi(rhoP, zP);

			if (ref_B_phi[i] == 0.0) {
				// exact zero has to be reproduced exactly
				Assertions.assertEquals(0.0, bPhi);
			} else {
				double relErr = Math.abs((bPhi - ref_B_phi[i])/ref_B_phi[i]);
				if (relErr >= tolerance) {
					System.out.printf(Locale.ENGLISH, "case %4d (rhoP=%.20e zP=%.20e) => relErr = %.3e\n",
							i, rhoP, zP, relErr);
				}
				Assertions.assertTrue(relErr < tolerance);
			}
		}
	}

	@Test
	public void testCircularWireLoop_A_phi() {
		final double tolerance = 1.0e-15;

		// load set of test points
		double[] testPointsRp = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/testPointsRpCircularWireLoop.dat")[0];
		double[] testPointsZp = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/testPointsZpCircularWireLoop.dat")[0];

		int numCases = testPointsRp.length;

		// load reference data
		double[] ref_A_phi = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/CircularWireLoop_A_phi_ref.dat")[0];

		// compute A_phi at test points and compare against reference
		for (int i=0; i<numCases; ++i) {
			double rhoP = testPointsRp[i];
			double zP   = testPointsZp[i];

			// A_phi is equal for sign change of z'

			double aPhi     = ABSCAB.circularWireLoop_A_phi(rhoP,  zP);
			double aPhiNegZ = ABSCAB.circularWireLoop_A_phi(rhoP, -zP);

			if (ref_A_phi[i] == 0.0) {
				// exact zero has to be reproduced exactly
				Assertions.assertEquals(0.0, aPhi);
				Assertions.assertEquals(0.0, aPhiNegZ);
			} else {
				double relErr = Math.abs((aPhi - ref_A_phi[i])/ref_A_phi[i]);
				if (relErr >= tolerance) {
					System.out.printf(Locale.ENGLISH, "case %4d (rhoP=%.20e zP=%.20e) => relErr = %.3e\n",
							i, rhoP, zP, relErr);
				}
				Assertions.assertTrue(relErr < tolerance);

				double relErrNegZ = Math.abs((aPhiNegZ - ref_A_phi[i])/ref_A_phi[i]);
				if (relErrNegZ >= tolerance) {
					System.out.printf(Locale.ENGLISH, "case %4d (rhoP=%.20e zP=%.20e) => relErr = %.3e\n",
							i, rhoP, -zP, relErrNegZ);
				}
				Assertions.assertTrue(relErrNegZ < tolerance);
			}
		}
	}

	@Test
	public void testCircularWireLoop_B_rho() {
		final double tolerance = 1.0e-13;

		// load set of test points
		double[] testPointsRp = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/testPointsRpCircularWireLoop.dat")[0];
		double[] testPointsZp = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/testPointsZpCircularWireLoop.dat")[0];

		int numCases = testPointsRp.length;

		// load reference data
		double[] ref_B_rho = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/CircularWireLoop_B_rho_ref.dat")[0];

		// compute B_rho at test points and compare against reference
		for (int i=0; i<numCases; ++i) {
			double rhoP = testPointsRp[i];
			double zP   = testPointsZp[i];

			// B_rho switches sign for sign change of z'

			double bRho     = ABSCAB.circularWireLoop_B_rho(rhoP,  zP);
			double bRhoNegZ = ABSCAB.circularWireLoop_B_rho(rhoP, -zP);

			if (ref_B_rho[i] == 0.0) {
				// exact zero has to be reproduced exactly
				Assertions.assertEquals(0.0, bRho);
				Assertions.assertEquals(0.0, bRhoNegZ);
			} else {
				double relErr = Math.abs((bRho - ref_B_rho[i])/ref_B_rho[i]);
				if (relErr >= tolerance) {
					System.out.printf(Locale.ENGLISH, "case %4d (rhoP=%.20e zP=%.20e) => relErr = %.3e\n",
							i, rhoP, zP, relErr);
				}
				Assertions.assertTrue(relErr < tolerance);

				double relErrNegZ = Math.abs((-bRhoNegZ - ref_B_rho[i])/ref_B_rho[i]);
				if (relErrNegZ >= tolerance) {
					System.out.printf(Locale.ENGLISH, "case %4d (rhoP=%.20e zP=%.20e) => relErr = %.3e\n",
							i, rhoP, -zP, relErrNegZ);
				}
				Assertions.assertTrue(relErrNegZ < tolerance);
			}
		}
	}

	@Test
	public void testCircularWireLoop_B_z() {
		final double tolerance = 1.0e-14;

		// load set of test points
		double[] testPointsRp = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/testPointsRpCircularWireLoop.dat")[0];
		double[] testPointsZp = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/testPointsZpCircularWireLoop.dat")[0];

		int numCases = testPointsRp.length;

		// load reference data
		double[] ref_B_z = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/CircularWireLoop_B_z_ref.dat")[0];

		// compute B_z at test points and compare against reference
		for (int i=0; i<numCases; ++i) {
			double rhoP = testPointsRp[i];
			double zP   = testPointsZp[i];

			// B_z is equal for sign change of z'

			double bZ     = ABSCAB.circularWireLoop_B_z(rhoP,  zP);
			double bZNegZ = ABSCAB.circularWireLoop_B_z(rhoP, -zP);

			if (ref_B_z[i] == 0.0) {
				// exact zero has to be reproduced exactly
				Assertions.assertEquals(0.0, bZ);
				Assertions.assertEquals(0.0, bZNegZ);
			} else {
				double relErr = Math.abs((bZ - ref_B_z[i])/ref_B_z[i]);
				if (relErr >= tolerance) {
					System.out.printf(Locale.ENGLISH, "case %4d (rhoP=%.20e zP=%.20e) => relErr = %.3e\n",
							i, rhoP, zP, relErr);
				}
				Assertions.assertTrue(relErr < tolerance);

				double relErrNegZ = Math.abs((bZNegZ - ref_B_z[i])/ref_B_z[i]);
				if (relErrNegZ >= tolerance) {
					System.out.printf(Locale.ENGLISH, "case %4d (rhoP=%.20e zP=%.20e) => relErr = %.3e\n",
							i, rhoP, -zP, relErrNegZ);
				}
				Assertions.assertTrue(relErrNegZ < tolerance);
			}
		}
	}
}
