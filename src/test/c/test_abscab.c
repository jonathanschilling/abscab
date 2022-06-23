
#include <stdio.h>

#include "util.h"
#include "abscab.h"

/** test method for Straight Wire Segment methods */
int testStraightWireSegment() {

	int rowsRp = 0;
	int colsRp = 0;
	double *test_points_rp = loadColumnsFromFile("../resources/testPointsRpStraightWireSegment.dat", &rowsRp, &colsRp)[0];
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
	double *test_points_zp = loadColumnsFromFile("../resources/testPointsZpStraightWireSegment.dat", &rowsZp, &colsZp)[0];
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
	double *aZRef = loadColumnsFromFile("../resources/StraightWireSegment_A_z_ref.dat", &rowsAZRef, &colsAZRef)[0];
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
	double *bPhiRef = loadColumnsFromFile("../resources/StraightWireSegment_B_phi_ref.dat", &rowsBPhiRef, &colsBPhiRef)[0];
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

		double rp = test_points_rp[i];
		double zp = test_points_zp[i];

		// compute values using C implementation to test
		double aZ   = straightWireSegment_A_z(rp, zp);
		double bPhi = straightWireSegment_B_phi(rp, zp);

		int aZStatus = assertRelAbsEquals(aZRef[i], aZ, toleranceAZ);
		if (aZStatus) {
			printf("error: mismatch at Straight Wire Segment A_z test case %d\n", i);
			printf("  rho' = %.17e\n", rp);
			printf("    z' = %.17e\n", zp);
			printf("  ref A_z = %+.17e\n", aZRef[i]);
			printf("  act A_z = %+.17e\n", aZ);
		}
		status |= aZStatus;

		int bPhiStatus = assertRelAbsEquals(bPhiRef[i], bPhi, toleranceBPhi);
		if (bPhiStatus) {
			printf("error: mismatch at Straight Wire Segment B_phi test case %d\n", i);
			printf("  rho' = %.17e\n", rp);
			printf("    z' = %.17e\n", zp);
			printf("  ref B_phi = %+.17e\n", bPhiRef[i]);
			printf("  act B_phi = %+.17e\n", bPhi);
		}
		status |= bPhiStatus;
	}
	return status;
}

/** test method for Circular Wire Loop methods */
int testCircularWireLoop() {

	int rowsRp = 0;
	int colsRp = 0;
	double *test_points_rp = loadColumnsFromFile("../resources/testPointsRpCircularWireLoop.dat", &rowsRp, &colsRp)[0];
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
	double *test_points_zp = loadColumnsFromFile("../resources/testPointsZpCircularWireLoop.dat", &rowsZp, &colsZp)[0];
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
	double *aPhiRef = loadColumnsFromFile("../resources/CircularWireLoop_A_phi_ref.dat", &rowsAPhiRef, &colsAPhiRef)[0];
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
	double *bRhoRef = loadColumnsFromFile("../resources/CircularWireLoop_B_rho_ref.dat", &rowsBRhoRef, &colsBRhoRef)[0];
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
	double *bZRef = loadColumnsFromFile("../resources/CircularWireLoop_B_z_ref.dat", &rowsBZRef, &colsBZRef)[0];
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

		double rp = test_points_rp[i];
		double zp = test_points_zp[i];

		// compute values using C implementation to test
		double aPhi = circularWireLoop_A_phi(rp, zp);
		double bRho = circularWireLoop_B_rho(rp, zp);
		double bZ   = circularWireLoop_B_z(rp, zp);

		int aPhiStatus = assertRelAbsEquals(aPhiRef[i], aPhi, toleranceAPhi);
		if (aPhiStatus) {
			printf("error: mismatch at Circular Wire Loop A_phi test case %d\n", i);
			printf("  rho' = %.17e\n", rp);
			printf("    z' = %.17e\n", zp);
			printf("  ref A_phi = %+.17e\n", aPhiRef[i]);
			printf("  act A_phi = %+.17e\n", aPhi);
		}
		status |= aPhiStatus;

		int bRhoStatus = assertRelAbsEquals(bRhoRef[i], bRho, toleranceBRho);
		if (bRhoStatus) {
			printf("error: mismatch at Circular Wire Loop B_rho test case %d\n", i);
			printf("  rho' = %.17e\n", rp);
			printf("    z' = %.17e\n", zp);
			printf("  ref B_rho = %+.17e\n", bRhoRef[i]);
			printf("  act B_rho = %+.17e\n", bRho);
		}
		status |= bRhoStatus;

		int bZStatus = assertRelAbsEquals(bZRef[i], bZ, toleranceBZ);
		if (bZStatus) {
			printf("error: mismatch at Circular Wire Loop B_z test case %d\n", i);
			printf("  rho' = %.17e\n", rp);
			printf("    z' = %.17e\n", zp);
			printf("  ref B_z = %+.17e\n", bZRef[i]);
			printf("  act B_z = %+.17e\n", bZ);
		}
		status |= bZStatus;
	}
	return status;
}

int main(int argc, char **argv) {

	int status = 0;

	status |= testStraightWireSegment();
	status |= testCircularWireLoop();

	if (status != 0) {
		printf("%s: some test(s) failed :-(\n", argv[0]);
	} else {
		printf("%s: all test(s) passed :-)\n", argv[0]);
	}
}
