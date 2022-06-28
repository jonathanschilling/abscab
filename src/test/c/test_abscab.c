
#include <stdio.h>

#include "util.h"
#include "abscab.h"

/** test method for Straight Wire Segment methods */
int testStraightWireSegment() {

	int rowsRp = 0;
	int colsRp = 0;
	double **test_points_rp = loadColumnsFromFile("../resources/testPointsRpStraightWireSegment.dat", &rowsRp, &colsRp);
	if (rowsRp < 1) {
		printf("error: need at least one row of test point coordinates for rho'\n");
		return 1;
	}
	if (colsRp < 1) {
		printf("error: need at least one column of test point coordinates for rho'\n");
		return 1;
	}

	int rowsZp = 0;
	int colsZp = 0;
	double **test_points_zp = loadColumnsFromFile("../resources/testPointsZpStraightWireSegment.dat", &rowsZp, &colsZp);
	if (rowsZp < 1) {
		printf("error: need at least one row of test point coordinates for z'\n");
		return 1;
	}
	if (colsZp < 1) {
		printf("error: need at least one column of test point coordinates for z'\n");
		return 1;
	}

	if (rowsRp != rowsZp) {
		printf("error: number of rows of test point coordinates has to agree between rho' (%d) and z' (%d)\n", rowsRp, rowsZp);
		return 1;
	}
	if (colsRp != colsZp) {
		printf("error: number of columns of test point coordinates has to agree between rho' (%d) and z' (%d)\n", colsRp, colsZp);
		return 1;
	}

	int numCases = rowsRp;

	int rowsAZRef = 0;
	int colsAZRef = 0;
	double **aZRef = loadColumnsFromFile("../resources/StraightWireSegment_A_z_ref.dat", &rowsAZRef, &colsAZRef);
	if (rowsAZRef < 1) {
		printf("error: need at least one row of reference values for A_z\n");
		return 1;
	}
	if (colsAZRef < 1) {
		printf("error: need at least one column of reference values for A_z\n");
		return 1;
	}

	if (numCases != rowsAZRef) {
		printf("error: number of reference values for A_z (%d) has to match number of test cases (%d)\n", rowsAZRef, rowsRp);
		return 1;
	}

	int rowsBPhiRef = 0;
	int colsBPhiRef = 0;
	double **bPhiRef = loadColumnsFromFile("../resources/StraightWireSegment_B_phi_ref.dat", &rowsBPhiRef, &colsBPhiRef);
	if (rowsBPhiRef < 1) {
		printf("error: need at least one row of reference values for B_phi\n");
		return 1;
	}
	if (colsBPhiRef < 1) {
		printf("error: need at least one column of reference values for B_phi\n");
		return 1;
	}

	if (numCases != rowsBPhiRef) {
		printf("error: number of reference values for B_phi (%d) has to match number of test cases (%d)\n", rowsBPhiRef, rowsRp);
		return 1;
	}

	double toleranceAZ   = 1.0e-15;
	double toleranceBPhi = 1.0e-15;

	int status = 0;
	for (int i = 0; i < numCases && !status; ++i) {

		double rp = test_points_rp[0][i];
		double zp = test_points_zp[0][i];

		// compute values using C implementation to test
		double aZ   = straightWireSegment_A_z(rp, zp);
		double bPhi = straightWireSegment_B_phi(rp, zp);

		int aZStatus = assertRelAbsEquals(aZRef[0][i], aZ, toleranceAZ);
		if (aZStatus) {
			printf("error: mismatch at Straight Wire Segment A_z test case %d\n", i);
			printf("  rho' = %.17e\n", rp);
			printf("    z' = %.17e\n", zp);
			printf("  ref A_z = %+.17e\n", aZRef[0][i]);
			printf("  act A_z = %+.17e\n", aZ);
		}
		status |= aZStatus;

		int bPhiStatus = assertRelAbsEquals(bPhiRef[0][i], bPhi, toleranceBPhi);
		if (bPhiStatus) {
			printf("error: mismatch at Straight Wire Segment B_phi test case %d\n", i);
			printf("  rho' = %.17e\n", rp);
			printf("    z' = %.17e\n", zp);
			printf("  ref B_phi = %+.17e\n", bPhiRef[0][i]);
			printf("  act B_phi = %+.17e\n", bPhi);
		}
		status |= bPhiStatus;

		if (status) {
			break;
		}
	}

	// free data loaded from text files
	free(test_points_rp[0]); free(test_points_rp);
	free(test_points_zp[0]); free(test_points_zp);
	free(aZRef[0]);          free(aZRef);
	free(bPhiRef[0]);        free(bPhiRef);

	return status;
}

/** test method for Circular Wire Loop methods */
int testCircularWireLoop() {

	int rowsRp = 0;
	int colsRp = 0;
	double **test_points_rp = loadColumnsFromFile("../resources/testPointsRpCircularWireLoop.dat", &rowsRp, &colsRp);
	if (rowsRp < 1) {
		printf("error: need at least one row of test point coordinates for rho'\n");
		return 1;
	}
	if (colsRp < 1) {
		printf("error: need at least one column of test point coordinates for rho'\n");
		return 1;
	}

	int rowsZp = 0;
	int colsZp = 0;
	double **test_points_zp = loadColumnsFromFile("../resources/testPointsZpCircularWireLoop.dat", &rowsZp, &colsZp);
	if (rowsZp < 1) {
		printf("error: need at least one row of test point coordinates for z'\n");
		return 1;
	}
	if (colsZp < 1) {
		printf("error: need at least one column of test point coordinates for z'\n");
		return 1;
	}

	if (rowsRp != rowsZp) {
		printf("error: number of rows of test point coordinates has to agree between rho' (%d) and z' (%d)\n", rowsRp, rowsZp);
		return 1;
	}
	if (colsRp != colsZp) {
		printf("error: number of columns of test point coordinates has to agree between rho' (%d) and z' (%d)\n", colsRp, colsZp);
		return 1;
	}

	int numCases = rowsRp;

	int rowsAPhiRef = 0;
	int colsAPhiRef = 0;
	double **aPhiRef = loadColumnsFromFile("../resources/CircularWireLoop_A_phi_ref.dat", &rowsAPhiRef, &colsAPhiRef);
	if (rowsAPhiRef < 1) {
		printf("error: need at least one row of reference values for A_phi\n");
		return 1;
	}
	if (colsAPhiRef < 1) {
		printf("error: need at least one column of reference values for A_phi\n");
		return 1;
	}

	if (numCases != rowsAPhiRef) {
		printf("error: number of reference values for A_phi (%d) has to match number of test cases (%d)\n", rowsAPhiRef, rowsRp);
		return 1;
	}

	int rowsBRhoRef = 0;
	int colsBRhoRef = 0;
	double **bRhoRef = loadColumnsFromFile("../resources/CircularWireLoop_B_rho_ref.dat", &rowsBRhoRef, &colsBRhoRef);
	if (rowsBRhoRef < 1) {
		printf("error: need at least one row of reference values for B_rho\n");
		return 1;
	}
	if (colsBRhoRef < 1) {
		printf("error: need at least one column of reference values for B_rho\n");
		return 1;
	}

	if (numCases != rowsBRhoRef) {
		printf("error: number of reference values for B_rho (%d) has to match number of test cases (%d)\n", rowsBRhoRef, rowsRp);
		return 1;
	}

	int rowsBZRef = 0;
	int colsBZRef = 0;
	double **bZRef = loadColumnsFromFile("../resources/CircularWireLoop_B_z_ref.dat", &rowsBZRef, &colsBZRef);
	if (rowsBZRef < 1) {
		printf("error: need at least one row of reference values for B_z\n");
		return 1;
	}
	if (colsBZRef < 1) {
		printf("error: need at least one column of reference values for B_z\n");
		return 1;
	}

	if (numCases != rowsBZRef) {
		printf("error: number of reference values for B_z (%d) has to match number of test cases (%d)\n", rowsBZRef, rowsRp);
		return 1;
	}

	double toleranceAPhi = 1.0e-15;
	double toleranceBRho = 1.0e-13;
	double toleranceBZ   = 1.0e-14;

	int status = 0;
	for (int i = 0; i < numCases && !status; ++i) {

		double rp = test_points_rp[0][i];
		double zp = test_points_zp[0][i];

		// compute values using C implementation to test
		double aPhi = circularWireLoop_A_phi(rp, zp);
		double bRho = circularWireLoop_B_rho(rp, zp);
		double bZ   = circularWireLoop_B_z(rp, zp);

		int aPhiStatus = assertRelAbsEquals(aPhiRef[0][i], aPhi, toleranceAPhi);
		if (aPhiStatus) {
			printf("error: mismatch at Circular Wire Loop A_phi test case %d\n", i);
			printf("  rho' = %.17e\n", rp);
			printf("    z' = %.17e\n", zp);
			printf("  ref A_phi = %+.17e\n", aPhiRef[0][i]);
			printf("  act A_phi = %+.17e\n", aPhi);
		}
		status |= aPhiStatus;

		int bRhoStatus = assertRelAbsEquals(bRhoRef[0][i], bRho, toleranceBRho);
		if (bRhoStatus) {
			printf("error: mismatch at Circular Wire Loop B_rho test case %d\n", i);
			printf("  rho' = %.17e\n", rp);
			printf("    z' = %.17e\n", zp);
			printf("  ref B_rho = %+.17e\n", bRhoRef[0][i]);
			printf("  act B_rho = %+.17e\n", bRho);
		}
		status |= bRhoStatus;

		int bZStatus = assertRelAbsEquals(bZRef[0][i], bZ, toleranceBZ);
		if (bZStatus) {
			printf("error: mismatch at Circular Wire Loop B_z test case %d\n", i);
			printf("  rho' = %.17e\n", rp);
			printf("    z' = %.17e\n", zp);
			printf("  ref B_z = %+.17e\n", bZRef[0][i]);
			printf("  act B_z = %+.17e\n", bZ);
		}
		status |= bZStatus;
	}

	// free data loaded from text files
	free(test_points_rp[0]); free(test_points_rp);
	free(test_points_zp[0]); free(test_points_zp);
	free(aPhiRef[0]);        free(aPhiRef);
	free(bRhoRef[0]);        free(bRhoRef);
	free(bZRef[0]);          free(bZRef);

	return status;
}

int testMagneticFieldInfiniteLineFilament() {
	double tolerance = 1.0e-15;

	// Demtroeder 2, Sec. 3.2.2 ("Magnetic field of a straight wire")
	// B(r) = mu_0 * I / (2 pi r)
	// Test this here with:
	// I = 123.0 A
	// r = 0.132 m
	// => B = 0.186 mT
	double current = 123.0;
	double r = 0.132;
	double bPhiRef = MU_0 * current / (2.0 * M_PI * r);
//	printf("ref bPhi = %.5e\n", bPhiRef);

	double vertices[] = {
			0.0, 0.0, -1.0e6,
			0.0, 0.0,  1.0e6
	};

	double evalPos[] = {
			r, 0.0, 0.0
	};

	// y component is B_phi
	double magneticField[3];
	magneticFieldPolygonFilament(2, vertices, current, 1, evalPos, magneticField);
	double bPhi = magneticField[1];
//	printf("act bPhi = %.5e\n", bPhi);

	double relAbsErr = fabs(bPhi - bPhiRef) / (1.0 + fabs(bPhiRef));
//	printf("raErr = %.5e\n", relAbsErr);

	return assertRelAbsEquals(bPhiRef, bPhi, tolerance);
}

int testBPhiInfiniteLineFilament() {
	double tolerance = 1.0e-15;

	// Demtroeder 2, Sec. 3.2.2 ("Magnetic field of a straight wire")
	// B(r) = mu_0 * I / (2 pi r)
	// Test this here with:
	// I = 123.0 A
	// r = 0.132 m
	// => B = 0.186 mT
	double current = 123.0;
	double r = 0.132;
	double bPhiRef = MU_0 * current / (2.0 * M_PI * r);
//	printf("ref bPhi = %.5e\n", bPhiRef);

	// half the length of the wire segment
	double halfL = 1e6;
	double L = 2*halfL;
	double rhoP = r / L;
	double zP = halfL / L;
	double bPhi = MU_0 * current / (4.0 * M_PI * L) * straightWireSegment_B_phi(rhoP, zP);
//	printf("act bPhi = %.5e\n", bPhi);

	double relAbsErr = fabs(bPhi - bPhiRef) / (1.0 + fabs(bPhiRef));
//	printf("raErr = %.5e\n", relAbsErr);

	return assertRelAbsEquals(bPhiRef, bPhi, tolerance);
}

int testMagneticFieldInsideLongCoil() {
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

	int N = 50000; // windings
	double L = 50.0; // total length of coil in m
	double n = N/L;

	double current = 10.0; // A
	double radius = 1.0; // m

	double bZ = 0.0;
	for (int i = 0; i < N; ++i) {

		// axial position of coil
		double z0 = -L/2.0 + (i + 0.5) / n;

		// compute magnetic field
		double prefac = MU_0 * current / (M_PI * radius);
		double bZContrib = prefac * circularWireLoop_B_z(0.0, z0);

//		printf("coil %d at z0 = % .3e => contrib = %.3e\n", i, z0, bZContrib);

		bZ += bZContrib;
	}
//	printf("B_z = %.5e\n", bZ);

	double relAbsErr = fabs(bZ - bZRef) / (1.0 + fabs(bZRef));
//	printf("raErr = %.5e\n", relAbsErr);

	return assertRelAbsEquals(bZRef, bZ, tolerance);
}

int main(int argc, char **argv) {

	int status = 0;

	status |= testStraightWireSegment();
	status |= testCircularWireLoop();
	status |= testMagneticFieldInfiniteLineFilament();
	status |= testBPhiInfiniteLineFilament();
	status |= testMagneticFieldInsideLongCoil();

	if (status != 0) {
		printf("%s: some test(s) failed :-(\n", argv[0]);
	} else {
		printf("%s: all test(s) passed :-)\n", argv[0]);
	}
}
