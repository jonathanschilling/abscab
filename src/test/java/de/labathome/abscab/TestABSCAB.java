package de.labathome.abscab;

import java.util.Locale;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

public class TestABSCAB {

	@Test
	public void testStraightWireSegment_A_z() {
		final double tolerance = 1.0e-15;

		// load set of test points
		double[] testPointsRp = Util.loadColumnsFromResource(DemoABSCAB.class, "/testPointsRpStraightWireSegment.dat")[0];
		double[] testPointsZp = Util.loadColumnsFromResource(DemoABSCAB.class, "/testPointsZpStraightWireSegment.dat")[0];

		int numCases = testPointsRp.length;

		// load reference data
		double[] ref_A_z = Util.loadColumnsFromResource(DemoABSCAB.class, "/StraightWireSegment_A_z_ref.dat")[0];

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
		double[] testPointsRp = Util.loadColumnsFromResource(DemoABSCAB.class, "/testPointsRpStraightWireSegment.dat")[0];
		double[] testPointsZp = Util.loadColumnsFromResource(DemoABSCAB.class, "/testPointsZpStraightWireSegment.dat")[0];

		int numCases = testPointsRp.length;

		// load reference data
		double[] ref_B_phi = Util.loadColumnsFromResource(DemoABSCAB.class, "/StraightWireSegment_B_phi_ref.dat")[0];

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
		double[] testPointsRp = Util.loadColumnsFromResource(DemoABSCAB.class, "/testPointsRpCircularWireLoop.dat")[0];
		double[] testPointsZp = Util.loadColumnsFromResource(DemoABSCAB.class, "/testPointsZpCircularWireLoop.dat")[0];

		int numCases = testPointsRp.length;

		// load reference data
		double[] ref_A_phi = Util.loadColumnsFromResource(DemoABSCAB.class, "/CircularWireLoop_A_phi_ref.dat")[0];

		// compute A_phi at test points and compare against reference
		for (int i=0; i<numCases; ++i) {
			double rhoP = testPointsRp[i];
			double zP   = testPointsZp[i];

			double aPhi = ABSCAB.circularWireLoop_A_phi(rhoP, zP);

			if (ref_A_phi[i] == 0.0) {
				// exact zero has to be reproduced exactly
				Assertions.assertEquals(0.0, aPhi);
			} else {
				double relErr = Math.abs((aPhi - ref_A_phi[i])/ref_A_phi[i]);
				if (relErr >= tolerance) {
					System.out.printf(Locale.ENGLISH, "case %4d (rhoP=%.20e zP=%.20e) => relErr = %.3e\n",
							i, rhoP, zP, relErr);
				}
				Assertions.assertTrue(relErr < tolerance);
			}
		}
	}

	@Test
	public void testCircularWireLoop_B_rho() {
		final double tolerance = 1.0e-13;

		// load set of test points
		double[] testPointsRp = Util.loadColumnsFromResource(DemoABSCAB.class, "/testPointsRpCircularWireLoop.dat")[0];
		double[] testPointsZp = Util.loadColumnsFromResource(DemoABSCAB.class, "/testPointsZpCircularWireLoop.dat")[0];

		int numCases = testPointsRp.length;

		// load reference data
		double[] ref_B_rho = Util.loadColumnsFromResource(DemoABSCAB.class, "/CircularWireLoop_B_rho_ref.dat")[0];

		// compute B_rho at test points and compare against reference
		for (int i=0; i<numCases; ++i) {
			double rhoP = testPointsRp[i];
			double zP   = testPointsZp[i];

			double bRho = ABSCAB.circularWireLoop_B_rho(rhoP, zP);

			if (ref_B_rho[i] == 0.0) {
				// exact zero has to be reproduced exactly
				Assertions.assertEquals(0.0, bRho);
			} else {
				double relErr = Math.abs((bRho - ref_B_rho[i])/ref_B_rho[i]);
				if (relErr >= tolerance) {
					System.out.printf(Locale.ENGLISH, "case %4d (rhoP=%.20e zP=%.20e) => relErr = %.3e\n",
							i, rhoP, zP, relErr);
				}
				Assertions.assertTrue(relErr < tolerance);
			}
		}
	}

	@Test
	public void testCircularWireLoop_B_z() {
		final double tolerance = 1.0e-14;

		// load set of test points
		double[] testPointsRp = Util.loadColumnsFromResource(DemoABSCAB.class, "/testPointsRpCircularWireLoop.dat")[0];
		double[] testPointsZp = Util.loadColumnsFromResource(DemoABSCAB.class, "/testPointsZpCircularWireLoop.dat")[0];

		int numCases = testPointsRp.length;

		// load reference data
		double[] ref_B_z = Util.loadColumnsFromResource(DemoABSCAB.class, "/CircularWireLoop_B_z_ref.dat")[0];

		// compute B_z at test points and compare against reference
		for (int i=0; i<numCases; ++i) {
			double rhoP = testPointsRp[i];
			double zP   = testPointsZp[i];

			double bZ = ABSCAB.circularWireLoop_B_z(rhoP, zP);

			if (ref_B_z[i] == 0.0) {
				// exact zero has to be reproduced exactly
				Assertions.assertEquals(0.0, bZ);
			} else {
				double relErr = Math.abs((bZ - ref_B_z[i])/ref_B_z[i]);
				if (relErr >= tolerance) {
					System.out.printf(Locale.ENGLISH, "case %4d (rhoP=%.20e zP=%.20e) => relErr = %.3e\n",
							i, rhoP, zP, relErr);
				}
				Assertions.assertTrue(relErr < tolerance);
			}
		}
	}
}
