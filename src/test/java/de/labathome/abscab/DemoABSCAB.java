package de.labathome.abscab;

import java.util.LinkedList;
import java.util.List;
import java.util.function.IntFunction;

import aliceinnets.python.jyplot.JyPlot;

public class DemoABSCAB {

	public static void main(String[] args) {

		demoCircularWireLoopTEAL();
//		demoStraightWireSegmentAtHalfHeight();

//		demoDoubleParts();

//		demoStraightWireSegmentAlongRhoP0();
//		demoStraightWireSegmentAlongZP01();

//		demoMcGreivy();
//		demoFiniteCoil();
//		demoAntiHelmholtzCoilField();
//		demoHelmholtzCoilField();
//		demoMagneticFieldOnAxisOfCircularWireLoop();

//		demoStraightWireSegment();
//		demoCircularWireLoop();

//		dumpInternalResultsStraightWireSegment();
//		dumpInternalResultsCircularWireLoop();
	}

	public static double cwl_A_phi_TEAL(double rhoP, double zP) {

		double kCSq_num = zP*zP + (1 - rhoP) * (1 - rhoP);
		double kCSq_den = zP*zP + (1 + rhoP) * (1 + rhoP);

		double prefac = 1.0 / Math.sqrt(kCSq_den);
		double kCSq = kCSq_num / kCSq_den;

		return prefac * CompleteEllipticIntegral.cel(Math.sqrt(kCSq), 1, -1, 1);
	}

	public static double cwl_B_rho_TEAL(double rhoP, double zP) {

		double kCSq_num = zP*zP + (1 - rhoP) * (1 - rhoP);
		double kCSq_den = zP*zP + (1 + rhoP) * (1 + rhoP);

		double prefac = zP / (kCSq_den * Math.sqrt(kCSq_den));
		double kCSq = kCSq_num / kCSq_den;

		return prefac * CompleteEllipticIntegral.cel(Math.sqrt(kCSq),kCSq, -1, 1);
	}

	public static double cwl_B_z_TEAL(double rhoP, double zP) {

		double kCSq_num = zP*zP + (1 - rhoP) * (1 - rhoP);
		double kCSq_den = zP*zP + (1 + rhoP) * (1 + rhoP);

		double prefac = 1.0 / (2 * rhoP * Math.sqrt(kCSq_den));
		double kCSq = kCSq_num / kCSq_den;

		double cel1 = CompleteEllipticIntegral.cel(Math.sqrt(kCSq), 1, -1, 1);

		double fac2 = (1 + kCSq - (1 - kCSq) * rhoP) / 2;
		double cel2 = CompleteEllipticIntegral.cel(Math.sqrt(kCSq),kCSq, -1, 1);

		return prefac * (cel1 + fac2 * cel2);
	}

	public static void demoCircularWireLoopTEAL() {

		// load set of test points
		double[] testPointsRp = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/testPointsRpCircularWireLoop.dat")[0];
		double[] testPointsZp = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/testPointsZpCircularWireLoop.dat")[0];

		int numCases = testPointsRp.length;

		// compute A_phi, B_rho and B_z at test points
		double[] A_phi_TEAL = new double[numCases];
		double[] B_rho_TEAL = new double[numCases];
		double[] B_z_TEAL = new double[numCases];

		for (int i=0; i<numCases; ++i) {
			double rhoP = testPointsRp[i];
			double zP   = testPointsZp[i];

			A_phi_TEAL[i] = cwl_A_phi_TEAL(rhoP, zP);
			B_rho_TEAL[i] = cwl_B_rho_TEAL(rhoP, zP);
			B_z_TEAL[i] = cwl_B_z_TEAL(rhoP, zP);
		}

		// write to output file
		UtilsTestABSCAB.dumpToFile(A_phi_TEAL, "data/CircularWireLoop_A_phi_TEAL.dat");
		UtilsTestABSCAB.dumpToFile(B_rho_TEAL, "data/CircularWireLoop_B_rho_TEAL.dat");
		UtilsTestABSCAB.dumpToFile(B_z_TEAL, "data/CircularWireLoop_B_z_TEAL.dat");
	}

	public static void demoStraightWireSegmentAtHalfHeight() {

		int n = 31;

		double rhoValues[] = new double[n];
		double sws_A_z_accurate[] = new double[n];
		double sws_A_z[] = new double[n];
		double error[] = new double[n];

		double z = 0.5;

		for (int i = 0; i < n; ++i) {
			rhoValues[i] = Math.pow(10, i - 15);

			sws_A_z_accurate[i] = ABSCAB.straightWireSegment_A_z(rhoValues[i], z);
			sws_A_z[i] = ABSCAB.sws_A_z_f(rhoValues[i], z);

			error[i] = UtilsTestABSCAB.errorMetric(sws_A_z_accurate[i], sws_A_z[i]);

			System.out.printf("rho' = %5.2e\n" +
					          "  sws_A_z_f = %.18e\n" +
					          "       real = %.18e\n" +
			                  "   => error = %d \n", rhoValues[i], sws_A_z[i], sws_A_z_accurate[i], (int) error[i]);
		}

		JyPlot plt = new JyPlot();

		plt.subplot(2,1,1);
		plt.semilogx(rhoValues, sws_A_z_accurate, "bo-", "label='correct'");
		plt.semilogx(rhoValues, sws_A_z,          "r.--", "label=r'$\\tilde{A}_{z,f}$'");
		plt.legend("loc='upper right'");
		plt.grid(true);
		plt.title("straight wire segment, $L=1m$, $I = 1A$, evaluated at $z=0.5m$");
		plt.ylabel("r'$\\tilde{A}_z$ / a.u.'");
		plt.xticks(new double[] {1.0e-15, 1.0e-12, 1.0e-9, 1.0e-6, 1.0e-3, 1.0, 1.0e3, 1.0e6, 1.0e9, 1.0e12, 1.0e15});
		plt.tick_params("axis='x', labelbottom=False");

		plt.subplot(2,1,2);
		plt.semilogx(rhoValues, error, "k.-", "label='error metric'");
		plt.legend("loc='upper right'");
		plt.grid(true);
		plt.xlabel("$\\rho ~ / ~ m$");
		plt.ylabel("error metric / 1");
		plt.xticks(new double[] {1.0e-15, 1.0e-12, 1.0e-9, 1.0e-6, 1.0e-3, 1.0, 1.0e3, 1.0e6, 1.0e9, 1.0e12, 1.0e15});

		plt.tight_layout();

		plt.show();
		plt.exec();
	}

	public static void demoDoubleParts() {
		int s = 0;
		int E = 127;
		int M = 1432344328;

		double f = Math.pow(-1, s) * Math.pow(2.0, E - 1023) * (1 + M / Math.pow(2, 52));

		int exponent = Math.getExponent(f) + 1023;

		long[] fParts = UtilsTestABSCAB.doubleParts(f);

		System.out.printf("s: %d =?= %d\n", s, fParts[0]);
		System.out.printf("E: %d =?= %d =?= %d\n", E, fParts[1], exponent);
		System.out.printf("M: %d =?= %d\n", M, fParts[2]);
	}

	public static void demoStraightWireSegmentAlongRhoP0() {

		// load set of test points
		double[] testPointsRp = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/testPointsRpStraightWireSegment.dat")[0];
		double[] testPointsZp = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/testPointsZpStraightWireSegment.dat")[0];

		int numCases = testPointsRp.length;

		// load reference data
		double[] ref_A_z   = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/StraightWireSegment_A_z_ref.dat")[0];

		List<Double> sws_A_z_refLst = new LinkedList<>();
		List<Double> sws_A_z_2aLst = new LinkedList<>();
		List<Double> sws_A_z_2bLst = new LinkedList<>();

		// compute B_phi at test points and compare against reference
		for (int i=0; i<numCases; ++i) {
			double rhoP = testPointsRp[i];
			double zP   = testPointsZp[i];

			if (rhoP == 0.0) {
				sws_A_z_refLst.add(ref_A_z[i]);

				sws_A_z_2aLst.add(ABSCAB.sws_A_z_ax_f(zP));
				sws_A_z_2bLst.add(ABSCAB.sws_A_z_ax_n(zP));
			}
		}

		int numCasesAlongRhoP0 = sws_A_z_refLst.size();
		double[] sws_A_z_2a_err = new double[numCasesAlongRhoP0];
		double[] sws_A_z_2b_err = new double[numCasesAlongRhoP0];
		for (int i=0; i<numCasesAlongRhoP0; ++i) {
			sws_A_z_2a_err[i] = UtilsTestABSCAB.errorMetric(sws_A_z_refLst.get(i), sws_A_z_2aLst.get(i));
			sws_A_z_2b_err[i] = UtilsTestABSCAB.errorMetric(sws_A_z_refLst.get(i), sws_A_z_2bLst.get(i));
		}

		JyPlot plt = new JyPlot();

		plt.figure("figsize=(6,2.5)");
		plt.plot(sws_A_z_2a_err, "r.-", "label='A_z_2'");
		plt.plot(sws_A_z_2b_err, "bx--", "label='A_z_2b'");
		plt.xlabel("test cases along $\\rho^{\\prime} = 0$");
		plt.ylabel("$log_{10}$(rel. error)");
		plt.grid(true);
		plt.legend("loc='center right'");
		plt.tight_layout();

		plt.show();
		plt.exec();
	}

	public static void demoStraightWireSegmentAlongZP01() {

		// load set of test points
		double[] testPointsRp = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/testPointsRpStraightWireSegment.dat")[0];
		double[] testPointsZp = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/testPointsZpStraightWireSegment.dat")[0];

		int numCases = testPointsRp.length;

		// load reference data
		double[] ref_A_z   = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/StraightWireSegment_A_z_ref.dat")[0];

		List<Double> sws_A_z_refLst = new LinkedList<>();
		List<Double> sws_A_z_3aLst = new LinkedList<>();
		List<Double> sws_A_z_3bLst = new LinkedList<>();

		// compute B_phi at test points and compare against reference
		for (int i=0; i<numCases; ++i) {
			double rhoP = testPointsRp[i];
			double zP   = testPointsZp[i];

			if (zP == 0.0) {
				sws_A_z_refLst.add(ref_A_z[i]);

				sws_A_z_3aLst.add(ABSCAB.sws_A_z_rad_f(rhoP));
				sws_A_z_3bLst.add(ABSCAB.sws_A_z_rad_n(rhoP));
			}
		}

		int numCasesAlongRhoP0 = sws_A_z_refLst.size();
		double[] sws_A_z_3a_err = new double[numCasesAlongRhoP0];
		double[] sws_A_z_3b_err = new double[numCasesAlongRhoP0];
		for (int i=0; i<numCasesAlongRhoP0; ++i) {
			sws_A_z_3a_err[i] = UtilsTestABSCAB.errorMetric(sws_A_z_refLst.get(i), sws_A_z_3aLst.get(i));
			sws_A_z_3b_err[i] = UtilsTestABSCAB.errorMetric(sws_A_z_refLst.get(i), sws_A_z_3bLst.get(i));
		}

		JyPlot plt = new JyPlot();

		plt.figure("figsize=(6,2.5)");
		plt.plot(sws_A_z_3a_err, "r.-", "label='A_z_3'");
		plt.plot(sws_A_z_3b_err, "bx--", "label='A_z_3b'");
		plt.xlabel("test cases along $z^{\\prime} = 0$");
		plt.ylabel("$log_{10}$(rel. error)");
		plt.grid(true);
		plt.legend("loc='center right'");
		plt.tight_layout();

		plt.show();
		plt.exec();

	}

	// 32 threads: 18 s --> equally fast as optimized C program !
	public static void demoMcGreivy() {

		long startTime = -System.nanoTime();

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

//		int numProcessors = 1;
		int numProcessors = Runtime.getRuntime().availableProcessors();

		boolean useCompensatedSummation = true;

		for (int i=0; i<numCases; ++i) {

			int numPhi = allNumPhi[i];
			System.out.printf("case %2d/%2d: numPhi = %d\n", i+1, numCases, numPhi);

			double omega = 2.0*Math.PI / (numPhi-1);
			IntFunction<double[]> vertexSupplierStd = idxVertex -> {
				double phi = omega * idxVertex;
				double x = radius * Math.cos(phi);
				double y = radius * Math.sin(phi);
				double z = 0.0;
				return new double[] {x, y, z};
			};
			double bZStd = ABSCAB.magneticFieldPolygonFilament(numPhi, vertexSupplierStd, current, evalPos, numProcessors, useCompensatedSummation)[2][0];

//			double[][] verticesStd = polygonCircleAround0(radius, numPhi);
//			double bZStd = ABSCAB.magneticFieldPolygonFilament(verticesStd, current, evalPos, numProcessors, useCompensatedSummation)[2][0];

			allBzStdErr[i] = UtilsTestABSCAB.errorMetric(bZRef, bZStd);
			System.out.printf("ABSCAB B_z = %.3e (err %g)\n", bZStd, allBzStdErr[i]);

			// McGreivy radius correction
			double dPhi = 2.0 * Math.PI / (numPhi - 1); // spacing between points

			// TODO: understand derivation of alpha for special case of closed circle
			// |dr/ds| = 2*pi
			// --> alpha = 1/R * (dr)^2 / 12
			// == 4 pi^2 / (12 R)
			//double rCorr = radius * (1.0 + 4 * Math.PI * Math.PI * dPhi*dPhi/ 12);
			double rCorr = radius * (1.0 + dPhi*dPhi/ 12);

			IntFunction<double[]> vertexSupplierMcGreivy = idxVertex -> {
				double phi = omega * idxVertex;
				double x = rCorr * Math.cos(phi);
				double y = rCorr * Math.sin(phi);
				double z = 0.0;
				return new double[] {x, y, z};
			};
			double bZMcG = ABSCAB.magneticFieldPolygonFilament(numPhi, vertexSupplierMcGreivy, current, evalPos, numProcessors, useCompensatedSummation)[2][0];

//			double[][] verticesMcG = polygonCircleAround0(rCorr, numPhi);
//			double bZMcG = ABSCAB.magneticFieldPolygonFilament(verticesMcG, current, evalPos, numProcessors, useCompensatedSummation)[2][0];

			allBzMcGErr[i] = UtilsTestABSCAB.errorMetric(bZRef, bZMcG);
			System.out.printf("McGrvy B_z = %.3e (err %g)\n", bZMcG, allBzMcGErr[i]);

			resultTable[0][i] = numPhi;
			resultTable[1][i] = allBzStdErr[i];
			resultTable[2][i] = allBzMcGErr[i];
		}

		if (useCompensatedSummation) {
			UtilsTestABSCAB.dumpToFile(resultTable, "data/convergenceMcGreivy_CompensatedSummation.dat");
		} else {
			UtilsTestABSCAB.dumpToFile(resultTable, "data/convergenceMcGreivy_StandardSummation.dat");
		}

		long duration = startTime + System.nanoTime();
		System.out.printf("duration: %.3f s\n", duration/1e9);

		JyPlot plt = new JyPlot();

		plt.figure();
		plt.semilogx(allNumPhi, allBzStdErr, ".-", "label='standard'");
		plt.semilogx(allNumPhi, allBzMcGErr, ".-", "label='McGreivy'");
		plt.grid(true);
		plt.xlabel("numPhi");
		plt.ylabel("rel. err");
		plt.legend("loc='upper right'");
		plt.title("McGreivy method for circular loop");
		plt.tight_layout();

		plt.show();
		plt.exec();
	}

	static double[][] polygonCircleAround0(double radius, int numPhi) {
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
		double[] testPointsRp = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/testPointsRpStraightWireSegment.dat")[0];
		double[] testPointsZp = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/testPointsZpStraightWireSegment.dat")[0];

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
		UtilsTestABSCAB.dumpToFile(A_z,   "data/StraightWireSegment_A_z_Java.dat");
		UtilsTestABSCAB.dumpToFile(B_phi, "data/StraightWireSegment_B_phi_Java.dat");
	}

	public static void demoCircularWireLoop() {
		// load set of test points
		double[] testPointsRp = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/testPointsRpCircularWireLoop.dat")[0];
		double[] testPointsZp = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/testPointsZpCircularWireLoop.dat")[0];

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
		UtilsTestABSCAB.dumpToFile(A_phi, "data/CircularWireLoop_A_phi_Java.dat");
		UtilsTestABSCAB.dumpToFile(B_rho, "data/CircularWireLoop_B_rho_Java.dat");
		UtilsTestABSCAB.dumpToFile(B_z,   "data/CircularWireLoop_B_z_Java.dat");
	}

	public static void dumpInternalResultsStraightWireSegment() {
		// load set of test points
		double[] testPointsRp = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/testPointsRpStraightWireSegment.dat")[0];
		double[] testPointsZp = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/testPointsZpStraightWireSegment.dat")[0];

		int numCases = testPointsRp.length;

		// compute A_z and B_phi at test points
		double[] A_z_ax  = new double[numCases];
		double[] A_z_rad = new double[numCases];
		double[] A_z_n   = new double[numCases];
		double[] A_z_f   = new double[numCases];

		double[] B_phi_rad = new double[numCases];
		double[] B_phi_f   = new double[numCases];
		double[] B_phi_n   = new double[numCases];

		for (int i=0; i<numCases; ++i) {
			double rhoP = testPointsRp[i];
			double zP   = testPointsZp[i];

			A_z_ax[i]  = ABSCAB.sws_A_z_ax(zP);
			A_z_rad[i] = ABSCAB.sws_A_z_rad(rhoP);
			A_z_n[i]   = ABSCAB.sws_A_z_n(rhoP, zP);
			A_z_f[i]   = ABSCAB.sws_A_z_f(rhoP, zP);

			B_phi_rad[i] = ABSCAB.sws_B_phi_rad(rhoP);
			B_phi_f[i]   = ABSCAB.sws_B_phi_f(rhoP, zP);
			B_phi_n[i]   = ABSCAB.sws_B_phi_n(rhoP, zP);
		}

		// write to output file
		UtilsTestABSCAB.dumpToFile(A_z_ax,  "data/StraightWireSegment_A_z_ax_Java.dat");
		UtilsTestABSCAB.dumpToFile(A_z_rad, "data/StraightWireSegment_A_z_rad_Java.dat");
		UtilsTestABSCAB.dumpToFile(A_z_n,   "data/StraightWireSegment_A_z_n_Java.dat");
		UtilsTestABSCAB.dumpToFile(A_z_f,   "data/StraightWireSegment_A_z_f_Java.dat");

		UtilsTestABSCAB.dumpToFile(B_phi_rad, "data/StraightWireSegment_B_phi_rad_Java.dat");
		UtilsTestABSCAB.dumpToFile(B_phi_f,   "data/StraightWireSegment_B_phi_f_Java.dat");
		UtilsTestABSCAB.dumpToFile(B_phi_n,   "data/StraightWireSegment_B_phi_n_Java.dat");
	}


	public static void dumpInternalResultsCircularWireLoop() {
		// load set of test points
		double[] testPointsRp = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/testPointsRpCircularWireLoop.dat")[0];
		double[] testPointsZp = UtilsTestABSCAB.loadColumnsFromResource(DemoABSCAB.class, "/testPointsZpCircularWireLoop.dat")[0];

		int numCases = testPointsRp.length;

		// compute A_phi, B_rho and B_z at test points
		double[] A_phi_f = new double[numCases];
		double[] A_phi_n = new double[numCases];
		double[] A_phi_v = new double[numCases];

		double[] B_rho_f = new double[numCases];
		double[] B_rho_n = new double[numCases];
		double[] B_rho_v = new double[numCases];

		double[] B_z_f1 = new double[numCases];
		double[] B_z_f2 = new double[numCases];
		double[] B_z_n  = new double[numCases];
		double[] B_z_v  = new double[numCases];

		for (int i=0; i<numCases; ++i) {
			double rhoP = testPointsRp[i];
			double zP   = testPointsZp[i];

			A_phi_f[i] = ABSCAB.cwl_A_phi_f(rhoP, zP);
			A_phi_n[i] = ABSCAB.cwl_A_phi_n(rhoP, zP);
			A_phi_v[i] = ABSCAB.cwl_A_phi_v(zP);

			B_rho_f[i] = ABSCAB.cwl_B_rho_f(rhoP, zP);
			B_rho_n[i] = ABSCAB.cwl_B_rho_n(rhoP, zP);
			B_rho_v[i] = ABSCAB.cwl_B_rho_v(zP);

			B_z_f1[i] = ABSCAB.cwl_B_z_f1(rhoP, zP);
			B_z_f2[i] = ABSCAB.cwl_B_z_f2(rhoP, zP);
			B_z_n[i]  = ABSCAB.cwl_B_z_n(rhoP, zP);
			B_z_v[i]  = ABSCAB.cwl_B_z_v(zP);
		}

		// write to output file
		UtilsTestABSCAB.dumpToFile(A_phi_f, "data/CircularWireLoop_A_phi_f_Java.dat");
		UtilsTestABSCAB.dumpToFile(A_phi_n, "data/CircularWireLoop_A_phi_n_Java.dat");
		UtilsTestABSCAB.dumpToFile(A_phi_v, "data/CircularWireLoop_A_phi_v_Java.dat");

		UtilsTestABSCAB.dumpToFile(B_rho_f, "data/CircularWireLoop_B_rho_f_Java.dat");
		UtilsTestABSCAB.dumpToFile(B_rho_n, "data/CircularWireLoop_B_rho_n_Java.dat");
		UtilsTestABSCAB.dumpToFile(B_rho_v, "data/CircularWireLoop_B_rho_v_Java.dat");

		UtilsTestABSCAB.dumpToFile(B_z_f1, "data/CircularWireLoop_B_z_f1_Java.dat");
		UtilsTestABSCAB.dumpToFile(B_z_f2, "data/CircularWireLoop_B_z_f2_Java.dat");
		UtilsTestABSCAB.dumpToFile(B_z_n,  "data/CircularWireLoop_B_z_n_Java.dat");
		UtilsTestABSCAB.dumpToFile(B_z_v,  "data/CircularWireLoop_B_z_v_Java.dat");
	}
}
