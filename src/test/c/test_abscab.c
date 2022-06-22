
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
		printf("error: need at least one row of reference values for A_z\n");
		return 1;
	}
	if (colsBPhiRef < 1) {
		printf("error: need at least one column of reference values for A_z\n");
		return 1;
	}

	if (numCases != rowsBPhiRef) {
		printf("error: number of reference values for B_phi (%d) has to match number of test cases (%d)\n", rowsBPhiRef, rowsRp);
		return 1;
	}


	double tolerance = 1.0e-15;
	int status = 0;

	for (int i = 0; i < numCases && !status; ++i) {

		double rp = test_points_rp[i];
		double zp = test_points_zp[i];

		// compute values using C implementation to test
		double aZ   = straightWireSegment_A_z(rp, zp);
		double bPhi = straightWireSegment_B_phi(rp, zp);

		int aZStatus = assertRelAbsEquals(aZRef[i],   aZ,   tolerance);
		if (aZStatus) {
			printf("error: mismatch at Straight Wire Segment A_z test case %d\n", i);
			printf("  rho' = %.17e\n", rp);
			printf("    z' = %.17e\n", zp);
			printf("  ref A_z = %+.17e\n", aZRef[i]);
			printf("  act A_z = %+.17e\n", aZ);
		}
		status |= aZStatus;

		int bPhiStatus = assertRelAbsEquals(bPhiRef[i], bPhi, tolerance);
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


	return 0;
}

int main(int argc, char** argv) {

	int status = 0;

	status |= testStraightWireSegment();
	status |= testCircularWireLoop();

	if (status != 0) {
		printf("%s: some test(s) failed :-(\n", argv[0]);
	} else {
		printf("%s: all test(s) passed :-)\n", argv[0]);
	}
}
