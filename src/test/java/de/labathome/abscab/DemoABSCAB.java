package de.labathome.abscab;

public class DemoABSCAB {

	public static void main(String[] args) {
		demoStraightWireSegment();
		demoCircularWireLoop();
	}

	public static void demoStraightWireSegment() {
		// load set of test points
		double[] testPointsRp = Util.loadColumnsFromFile("src/test/resources/testPointsRpStraightWireSegment.dat")[0];
		double[] testPointsZp = Util.loadColumnsFromFile("src/test/resources/testPointsZpStraightWireSegment.dat")[0];

		int numCases = testPointsRp.length;

		// compute A_z and B_phi at test points
		double[] A_z = new double[numCases];
		double[] B_phi = new double[numCases];

		for (int i=0; i<numCases; ++i) {
			double rhoP = testPointsRp[i];
			double zP = testPointsZp[i];

			A_z[i] = ABSCAB.straightWireSegment_A_z(rhoP, zP);
			B_phi[i] = ABSCAB.straightWireSegment_B_phi(rhoP, zP);
		}

		// write to output file
		Util.dumpToFile(A_z,   "data/StraightWireSegment_A_z_Java.dat");
		Util.dumpToFile(B_phi, "data/StraightWireSegment_B_phi_Java.dat");
	}

	public static void demoCircularWireLoop() {
		// load set of test points
		double[] testPointsRp = Util.loadColumnsFromFile("src/test/resources/testPointsRpCircularWireLoop.dat")[0];
		double[] testPointsZp = Util.loadColumnsFromFile("src/test/resources/testPointsZpCircularWireLoop.dat")[0];

		int numCases = testPointsRp.length;

		// compute A_phi, B_rho and B_z at test points
		double[] A_phi = new double[numCases];
		double[] B_rho = new double[numCases];
		double[] B_z = new double[numCases];

		for (int i=0; i<numCases; ++i) {
			double rhoP = testPointsRp[i];
			double zP = testPointsZp[i];

			A_phi[i] = ABSCAB.circularWireLoop_A_phi(rhoP, zP);
			B_rho[i] = ABSCAB.circularWireLoop_B_rho(rhoP, zP);
			B_z[i] = ABSCAB.circularWireLoop_B_z(rhoP, zP);
		}

		// write to output file
		Util.dumpToFile(A_phi, "data/CircularWireLoop_A_phi_Java.dat");
		Util.dumpToFile(B_rho, "data/CircularWireLoop_B_rho_Java.dat");
		Util.dumpToFile(B_z,   "data/CircularWireLoop_B_z_Java.dat");
	}
}
