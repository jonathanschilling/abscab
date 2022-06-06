package de.labathome.abscab;

import aliceinnets.python.jyplot.JyPlot;

public class DemoABSCAB {

	public static void main(String[] args) {
		//		run();

		demoMcGreivy();
//		demoFiniteCoil();
//		demoAntiHelmholtzCoilField();
//		demoHelmholtzCoilField();
//		demoMagneticFieldOnAxisOfCircularWireLoop();

//		demoStraightWireSegment();
//		demoCircularWireLoop();
//
//		dumpInternalResultsStraightWireSegment();
//		dumpInternalResultsCircularWireLoop();
	}

	public static void run() {
		int s = 0;
		int E = 127;
		int M = 1432344328;

		double f = Math.pow(-1, s) * Math.pow(2.0, E - 1023) * (1 + M / Math.pow(2, 52));

		int exponent = Math.getExponent(f) + 1023;

		long[] fParts = Util.doubleParts(f);

		System.out.printf("s: %d =?= %d\n", s, fParts[0]);
		System.out.printf("E: %d =?= %d =?= %d\n", E, fParts[1], exponent);
		System.out.printf("M: %d =?= %d\n", M, fParts[2]);
	}

	public static void demoMcGreivy() {

		double radius = 1.23; // m
		double current = 17.0; // A

		double[] center = { 0.0, 0.0, 0.0 };
		double[] normal = { 0.0, 0.0, 1.0 };

		double[][] evalPos = {
				{10.0},
				{5.0},
				{0.0}
		};

		double bZRef = ABSCAB.magneticFieldCircularFilament(center, normal, radius, current, evalPos)[2][0];
		System.out.printf("ref B_z = %.3e\n", bZRef);

		// mimic circular wire loop as:
		// a) Polygon with points on the circule to be mimiced
		// b) Polygon with points slightly offset radially outward (McGreivy correction)
		// --> a) should have 2nd-order convergence;
		//     b) should have 4th-order convergence wrt. number of Polygon points

		int[] allNumPhi = {
				10, 30, 100, 300, 1000, 3000,
				10_000, 30_000, 100_000, 300_000,
				1_000_000, 3_000_000,
				10_000_000, 30_000_000,
				100_000_000, 300_000_000, 1_000_000_000
		};
//		int[] allNumPhi = {10, 30, 100, 300, 1000, 3000};
		int numCases = allNumPhi.length;

		double[] allBzStdErr = new double[numCases];
		double[] allBzMcGErr = new double[numCases];

		double[][] resultTable = new double[3][numCases];

		int numProcessors = Runtime.getRuntime().availableProcessors();
		boolean useCompensatedSummation = false;

		for (int i=0; i<numCases; ++i) {

			int numPhi = allNumPhi[i];
			System.out.printf("numPhi = %d\n", numPhi);

			double[][] verticesStd = polygonCircleAround0(radius, numPhi);
			double bZStd = ABSCAB.magneticFieldPolygonFilament(verticesStd, current, evalPos, numProcessors, useCompensatedSummation)[2][0];
			System.out.printf("ABSCAB B_z = %.3e\n", bZStd);
			allBzStdErr[i] = Math.abs((bZStd - bZRef)/bZRef);

			// McGreivy radius correction
			double dPhi = 2.0 * Math.PI / (numPhi - 1); // spacing between points

			// TODO: understand derivation of alpha for special case of closed circle
			// |dr/ds| = 2*pi
			// --> alpha = 1/R * (dr)^2 / 12
			// == 4 pi^2 / (12 R)
			//double rCorr = radius * (1.0 + 4 * Math.PI * Math.PI * dPhi*dPhi/ 12);
			double rCorr = radius * (1.0 + dPhi*dPhi/ 12);

			double[][] verticesMcG = polygonCircleAround0(rCorr, numPhi);
			double bZMcG = ABSCAB.magneticFieldPolygonFilament(verticesMcG, current, evalPos, numProcessors, useCompensatedSummation)[2][0];
			System.out.printf("McGrvy B_z = %.3e\n", bZMcG);

			allBzMcGErr[i] = Math.abs((bZMcG - bZRef)/bZRef);

			resultTable[0][i] = numPhi;
			resultTable[1][i] = allBzStdErr[i];
			resultTable[2][i] = allBzMcGErr[i];
		}

		if (useCompensatedSummation) {
			Util.dumpToFile(resultTable, "data/convergenceMcGreivy_StandardSummation.dat");
		} else {
			Util.dumpToFile(resultTable, "data/convergenceMcGreivy_CompensatedSummation.dat");
		}

		JyPlot plt = new JyPlot();

		plt.figure();
		plt.loglog(allNumPhi, allBzStdErr, ".-", "label='standard'");
		plt.loglog(allNumPhi, allBzMcGErr, ".-", "label='McGreivy'");
		plt.grid(true);
		plt.xlabel("numPhi");
		plt.ylabel("rel. err");
		plt.legend("loc='upper right'");
		plt.title("McGreivy method for circular loop");
		plt.ylim(new double[] {1e-16, 1});
		plt.tight_layout();

		plt.show();
		plt.exec();
	}

	private static double[][] polygonCircleAround0(double radius, int numPhi) {
		double[][] ret = new double[3][numPhi];
		double omega = 2.0*Math.PI / (numPhi-1);
		for (int i=0; i<numPhi; ++i) {
			double phi = omega * i;
			ret[0][i] = radius * Math.cos(phi);
			ret[1][i] = radius * Math.sin(phi);
		}

		return ret;
	}

	public static void demoFiniteCoil() {
		// Demtroeder 2, Sec. 3.2.6d ("Magnetic field of a cylindrical coil")

		double radius = 1.23; // m
		double current = 17.0; // A
		int N = 1000; // windings of coil

		int n = 100;
		double[] zP = new double[n];

		// coil with aspect ratio L = 6 R
		double[] bZRef6 = new double[n];
		double[] bZ6 = new double[n];

		// coil with aspect ratio L = 12 R
		double[] bZRef12 = new double[n];
		double[] bZ12 = new double[n];

		for (int i=0; i<n; ++i) {

			// axial evaluation position
			zP[i] = -10.0 + i * 20.0/(n-1);

			// reference from Demtroeder
			double coilLength6 = 6.0 * radius;
			double windingDensity6 = N/coilLength6;
			double prefac6 = ABSCAB.MU_0 * windingDensity6 * current / 2.0;
			double t1_6 = zP[i]*radius + coilLength6/2.0;
			double t2_6 = zP[i]*radius - coilLength6/2.0;
			bZRef6[i] = prefac6 * (t1_6/Math.sqrt(radius*radius + t1_6*t1_6)
					- t2_6/Math.sqrt(radius*radius + t2_6*t2_6));

			double coilLength12 = 12.0 * radius;
			double windingDensity12 = N/coilLength12;
			double prefac12 = ABSCAB.MU_0 * windingDensity12 * current / 2.0;
			double t1_12 = zP[i]*radius + coilLength12/2.0;
			double t2_12 = zP[i]*radius - coilLength12/2.0;
			bZRef12[i] = prefac12 * (t1_12/Math.sqrt(radius*radius + t1_12*t1_12)
					- t2_12/Math.sqrt(radius*radius + t2_12*t2_12));

			// eval using ABSCAB
			bZ6[i]  = B_z_coil(zP[i], 6.0, radius, N, current);
			bZ12[i] = B_z_coil(zP[i], 12.0, radius, N, current);
		}

		JyPlot plt = new JyPlot();

		plt.figure();
		plt.plot(zP, bZRef6, "o-", "label='ref 6'");
		plt.plot(zP, bZ6,    "x--", "label='ABSCAB 6'");
		plt.plot(zP, bZRef12, "o-", "label='ref 12'");
		plt.plot(zP, bZ12,    "x--", "label='ABSCAB 12'");
		plt.grid(true);
		plt.legend("loc='upper right'");
		plt.xlabel("z / R");
		plt.ylabel("B_z / T");
		plt.ticklabel_format("axis='y', style='sci', scilimits=(-2,2)");
		plt.title("B_z along the axis of a cylindrical coil");
		plt.tight_layout();

		plt.show();
		plt.exec();
	}

	/**
	 *
	 * @param zP z/R eval position
	 * @param coilAspectRatio L/R aspect ratio of coil
	 * @param radius R radius of coil
	 * @param N number of windings
	 * @return
	 */
	private static double B_z_coil(double zP, double coilAspectRatio, double radius, int N, double current) {

		double z = zP * radius; // real-space axial eval position in m
		double L = radius * coilAspectRatio; // length of coil in m
		double n = N/L; // winding density in 1/m

		double bZ = 0.0;
		for (int i=0; i<N; ++i) {

			// axial position of i:th winding
			double z0 = -L/2.0 + (i + 0.5) / n;

			// compute magnetic field
			double prefac = ABSCAB.MU_0 * current / (Math.PI * radius);
			double bZContrib = prefac * ABSCAB.circularWireLoop_B_z(0.0, (z - z0)/radius);

			bZ += bZContrib;
		}
		return bZ;
	}

	public static void demoAntiHelmholtzCoilField() {
		double current = 123.0; // A
		double radius = 0.2; // m
		double z0 = -0.1; // m
		double z1 =  0.1; // m

		double[] normal = { 0.0, 0.0, 1.0 };

		int n = 100;
		double deltaZ = (z1-z0)/(n-1);

		double[] z = new double[n];
		double[] B_z_ref = new double[n];
		double[] B_z = new double[n];
		for (int i = 0; i<n; ++i) {
			z[i] = z0 + i * deltaZ;

			// from Demtroeder 2, Sec. 3.2.6c
			double prefac = 48/(25*Math.sqrt(5.0)) * ABSCAB.MU_0 * current / (radius*radius);
			B_z_ref[i] = prefac * z[i];

			double[][] evalPos = {
					{0.0},
					{0.0},
					{z[i]}
			};

			double[] center1 = { 0.0, 0.0, -radius/2.0 };
			double B_z_1 = ABSCAB.magneticFieldCircularFilament(center1, normal, radius, -current, evalPos)[2][0];

			double[] center2 = { 0.0, 0.0,  radius/2.0 };
			double B_z_2 = ABSCAB.magneticFieldCircularFilament(center2, normal, radius,  current, evalPos)[2][0];

			B_z[i] = B_z_1 + B_z_2;
		}

		JyPlot plt = new JyPlot();

		plt.figure();
		plt.plot(z, B_z_ref, "o-", "label='ref'");
		plt.plot(z, B_z,     "x--", "label='ABSCAB'");
		plt.grid(true);
		plt.legend("loc='upper right'");
		plt.xlabel("z / m");
		plt.ylabel("B_z / T");
		plt.ticklabel_format("axis='y', style='sci', scilimits=(-2,2)");
		plt.title("B_z along the axis of an Anti-Helmholtz coil pair");
		plt.tight_layout();

		plt.show();
		plt.exec();
	}

	public static void demoHelmholtzCoilField() {
		double current = 123.0; // A
		double radius = 0.2; // m
		double z0 = -0.1; // m
		double z1 =  0.1; // m

		double[] normal = { 0.0, 0.0, 1.0 };

		int n = 100;
		double deltaZ = (z1-z0)/(n-1);

		double[] z = new double[n];
		double[] B_z_ref = new double[n];
		double[] B_z = new double[n];
		for (int i = 0; i<n; ++i) {
			z[i] = z0 + i * deltaZ;

			// from Demtroeder 2, Sec. 3.2.6c
			double fiveFourth = 5.0 / 4.0;
			double prefac = ABSCAB.MU_0 * current / (Math.sqrt(fiveFourth)*fiveFourth * radius);
			double zP = z[i]/radius;
			B_z_ref[i] = prefac * (1.0 - 144.0/125.0 * zP*zP * zP*zP);

			double[][] evalPos = {
					{0.0},
					{0.0},
					{z[i]}
			};

			double[] center1 = { 0.0, 0.0, -radius/2.0 };
			double B_z_1 = ABSCAB.magneticFieldCircularFilament(center1, normal, radius, current, evalPos)[2][0];

			double[] center2 = { 0.0, 0.0,  radius/2.0 };
			double B_z_2 = ABSCAB.magneticFieldCircularFilament(center2, normal, radius, current, evalPos)[2][0];

			B_z[i] = B_z_1 + B_z_2;
		}

		JyPlot plt = new JyPlot();

		plt.figure();
		plt.plot(z, B_z_ref, "o-", "label='ref'");
		plt.plot(z, B_z,     "x--", "label='ABSCAB'");
		plt.grid(true);
		plt.legend("loc='upper right'");
		plt.xlabel("z / m");
		plt.ylabel("B_z / T");
		plt.ticklabel_format("axis='y', style='sci', scilimits=(-2,2)");
		plt.title("B_z along the axis of a Helmholtz coil pair");
		plt.tight_layout();

		plt.show();
		plt.exec();
	}

	public static void demoMagneticFieldOnAxisOfCircularWireLoop() {

		double current = 123.0; // A
		double radius = 0.2; // m
		double z0 = -1.0; // m
		double z1 =  1.0; // m

		double[] center = { 0.0, 0.0, 0.0 };
		double[] normal = { 0.0, 0.0, 1.0 };

		int n = 100;
		double deltaZ = (z1-z0)/(n-1);

		double[] z = new double[n];
		double[] B_z_ref = new double[n];
		double[] B_z = new double[n];
		for (int i = 0; i<n; ++i) {
			z[i] = z0 + i * deltaZ;

			// from Demtroeder 2, Sec. 3.2.6b
			double d = z[i] * z[i] + radius*radius;
			B_z_ref[i] = ABSCAB.MU_0 * current * radius*radius / (2 * Math.sqrt(d)*d);

			double[][] evalPos = {
					{0.0},
					{0.0},
					{z[i]}
			};

			B_z[i] = ABSCAB.magneticFieldCircularFilament(center, normal, radius, current, evalPos)[2][0];
		}

		JyPlot plt = new JyPlot();

		plt.figure();
		plt.plot(z, B_z_ref, "o-", "label='ref'");
		plt.plot(z, B_z,     "x--", "label='ABSCAB'");
		plt.grid(true);
		plt.legend("loc='upper right'");
		plt.xlabel("z / m");
		plt.ylabel("B_z / T");
		plt.ticklabel_format("axis='y', style='sci', scilimits=(-2,2)");
		plt.title("B_z along the axis of a circular wire loop");
		plt.tight_layout();

		plt.show();
		plt.exec();
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
