package de.labathome.abscab;

public class DemoABSCAB {

	public static void main(String[] args) {
		demoStraightWireSegment();
		demoCircularWireLoop();

		dumpInternalResultsStraightWireSegment();
		dumpInternalResultsCircularWireLoop();
	}

	public static void demoStraightWireSegment() {
		// load set of test points
		double[] testPointsRp = Util.loadColumnsFromResource(DemoABSCAB.class, "/testPointsRpStraightWireSegment.dat")[0];
		double[] testPointsZp = Util.loadColumnsFromResource(DemoABSCAB.class, "/testPointsZpStraightWireSegment.dat")[0];

		int numCases = testPointsRp.length;

		// compute A_z and B_phi at test points
		double[] A_z   = new double[numCases];
		double[] B_phi = new double[numCases];

		for (int i=0; i<numCases; ++i) {
			double rhoP = testPointsRp[i];
			double zP   = testPointsZp[i];

			A_z[i]   = ABSCAB.straightWireSegment_A_z(rhoP, zP);
			B_phi[i] = ABSCAB.straightWireSegment_B_phi(rhoP, zP);
		}

		// write to output file
		Util.dumpToFile(A_z,   "data/StraightWireSegment_A_z_Java.dat");
		Util.dumpToFile(B_phi, "data/StraightWireSegment_B_phi_Java.dat");
	}

	public static void demoCircularWireLoop() {
		// load set of test points
		double[] testPointsRp = Util.loadColumnsFromResource(DemoABSCAB.class, "/testPointsRpCircularWireLoop.dat")[0];
		double[] testPointsZp = Util.loadColumnsFromResource(DemoABSCAB.class, "/testPointsZpCircularWireLoop.dat")[0];

		int numCases = testPointsRp.length;

		// compute A_phi, B_rho and B_z at test points
		double[] A_phi = new double[numCases];
		double[] B_rho = new double[numCases];
		double[] B_z   = new double[numCases];

		for (int i=0; i<numCases; ++i) {
			double rhoP = testPointsRp[i];
			double zP   = testPointsZp[i];

			A_phi[i] = ABSCAB.circularWireLoop_A_phi(rhoP, zP);
			B_rho[i] = ABSCAB.circularWireLoop_B_rho(rhoP, zP);
			B_z[i]   = ABSCAB.circularWireLoop_B_z(rhoP, zP);
		}

		// write to output file
		Util.dumpToFile(A_phi, "data/CircularWireLoop_A_phi_Java.dat");
		Util.dumpToFile(B_rho, "data/CircularWireLoop_B_rho_Java.dat");
		Util.dumpToFile(B_z,   "data/CircularWireLoop_B_z_Java.dat");
	}

	public static void dumpInternalResultsStraightWireSegment() {
		// load set of test points
		double[] testPointsRp = Util.loadColumnsFromResource(DemoABSCAB.class, "/testPointsRpStraightWireSegment.dat")[0];
		double[] testPointsZp = Util.loadColumnsFromResource(DemoABSCAB.class, "/testPointsZpStraightWireSegment.dat")[0];

		int numCases = testPointsRp.length;

		// compute A_z and B_phi at test points
		double[] A_z_along_rhoP_0    = new double[numCases];
		double[] A_z_along_zP_0_or_1 = new double[numCases];
		double[] A_z_6a              = new double[numCases];
		double[] A_z_6b              = new double[numCases];
		double[] A_z_6c              = new double[numCases];
		double[] A_z_1               = new double[numCases];

		double[] B_phi_2 = new double[numCases];
		double[] B_phi_3 = new double[numCases];
		double[] B_phi_4 = new double[numCases];
		double[] B_phi_5 = new double[numCases];

		for (int i=0; i<numCases; ++i) {
			double rhoP = testPointsRp[i];
			double zP   = testPointsZp[i];

			A_z_along_rhoP_0[i]    = ABSCAB.A_z_along_rhoP_0(rhoP, zP);
			A_z_along_zP_0_or_1[i] = ABSCAB.A_z_along_zP_0_or_1(rhoP, zP);
			A_z_6a[i]              = ABSCAB.A_z_6a(rhoP, zP);
			A_z_6b[i]              = ABSCAB.A_z_6b(rhoP, zP);
			A_z_6c[i]              = ABSCAB.A_z_6c(rhoP, zP);
			A_z_1[i]               = ABSCAB.A_z_1(rhoP, zP);

			B_phi_2[i] = ABSCAB.B_phi_2(rhoP, zP);
			B_phi_3[i] = ABSCAB.B_phi_3(rhoP, zP);
			B_phi_4[i] = ABSCAB.B_phi_4(rhoP, zP);
			B_phi_5[i] = ABSCAB.B_phi_5(rhoP, zP);
		}

		// write to output file
		Util.dumpToFile(A_z_along_rhoP_0,    "data/StraightWireSegment_A_z_along_rhoP_0_Java.dat");
		Util.dumpToFile(A_z_along_zP_0_or_1, "data/StraightWireSegment_A_z_along_zP_0_or_1_Java.dat");
		Util.dumpToFile(A_z_6a,              "data/StraightWireSegment_A_z_6a_Java.dat");
		Util.dumpToFile(A_z_6b,              "data/StraightWireSegment_A_z_6b_Java.dat");
		Util.dumpToFile(A_z_6c,              "data/StraightWireSegment_A_z_6c_Java.dat");
		Util.dumpToFile(A_z_1,               "data/StraightWireSegment_A_z_1_Java.dat");

		Util.dumpToFile(B_phi_2, "data/StraightWireSegment_B_phi_2_Java.dat");
		Util.dumpToFile(B_phi_3, "data/StraightWireSegment_B_phi_3_Java.dat");
		Util.dumpToFile(B_phi_4, "data/StraightWireSegment_B_phi_4_Java.dat");
		Util.dumpToFile(B_phi_5, "data/StraightWireSegment_B_phi_5_Java.dat");
	}


	public static void dumpInternalResultsCircularWireLoop() {
		// load set of test points
		double[] testPointsRp = Util.loadColumnsFromResource(DemoABSCAB.class, "/testPointsRpCircularWireLoop.dat")[0];
		double[] testPointsZp = Util.loadColumnsFromResource(DemoABSCAB.class, "/testPointsZpCircularWireLoop.dat")[0];

		int numCases = testPointsRp.length;

		// compute A_phi, B_rho and B_z at test points
		double[] A_phi_1 = new double[numCases];
		double[] A_phi_6 = new double[numCases];
		double[] A_phi_5 = new double[numCases];

		double[] B_rho_3 = new double[numCases];
		double[] B_rho_1 = new double[numCases];
		double[] B_rho_4 = new double[numCases];

		double[] B_z_1 = new double[numCases];
		double[] B_z_2 = new double[numCases];
		double[] B_z_4 = new double[numCases];
		double[] B_z_5 = new double[numCases];
		double[] B_z_6 = new double[numCases];

		for (int i=0; i<numCases; ++i) {
			double rhoP = testPointsRp[i];
			double zP   = testPointsZp[i];

			A_phi_1[i] = ABSCAB.A_phi_1(rhoP, zP);
			A_phi_6[i] = ABSCAB.A_phi_6(rhoP, zP);
			A_phi_5[i] = ABSCAB.A_phi_5(rhoP, zP);

			B_rho_3[i] = ABSCAB.B_rho_3(rhoP, zP);
			B_rho_1[i] = ABSCAB.B_rho_1(rhoP, zP);
			B_rho_4[i] = ABSCAB.B_rho_4(rhoP, zP);

			B_z_1[i] = ABSCAB.B_z_1(rhoP, zP);
			B_z_2[i] = ABSCAB.B_z_2(rhoP, zP);
			B_z_4[i] = ABSCAB.B_z_4(rhoP, zP);
			B_z_5[i] = ABSCAB.B_z_5(rhoP, zP);
			B_z_6[i] = ABSCAB.B_z_6(rhoP, zP);
		}

		// write to output file
		Util.dumpToFile(A_phi_1, "data/CircularWireLoop_A_phi_1_Java.dat");
		Util.dumpToFile(A_phi_6, "data/CircularWireLoop_A_phi_6_Java.dat");
		Util.dumpToFile(A_phi_5, "data/CircularWireLoop_A_phi_5_Java.dat");

		Util.dumpToFile(B_rho_3, "data/CircularWireLoop_B_rho_3_Java.dat");
		Util.dumpToFile(B_rho_1, "data/CircularWireLoop_B_rho_1_Java.dat");
		Util.dumpToFile(B_rho_4, "data/CircularWireLoop_B_rho_4_Java.dat");

		Util.dumpToFile(B_z_1,   "data/CircularWireLoop_B_z_1_Java.dat");
		Util.dumpToFile(B_z_2,   "data/CircularWireLoop_B_z_2_Java.dat");
		Util.dumpToFile(B_z_4,   "data/CircularWireLoop_B_z_4_Java.dat");
		Util.dumpToFile(B_z_5,   "data/CircularWireLoop_B_z_5_Java.dat");
		Util.dumpToFile(B_z_6,   "data/CircularWireLoop_B_z_6_Java.dat");
	}
}
