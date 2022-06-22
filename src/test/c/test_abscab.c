
#include <stdio.h>

#include "util.h"

#include "abscab.h"

/** test case for cel() implementation as described in section 4.2 of the 1969 Bulirsch article */
int test_abscab() {

	int rowsRp = 0;
	int colsRp = 0;
	double **test_points_rp = loadColumnsFromFile("../resources/testPointsRpStraightWireSegment.dat", &rowsRp, &colsRp);

	int rowsZp = 0;
	int colsZp = 0;
	double **test_points_zp = loadColumnsFromFile("../resources/testPointsZpStraightWireSegment.dat", &rowsZp, &colsZp);

	double aPhi = circularWireLoop_A_phi(1.1, 3.0);
	printf("A_phi = %g\n", aPhi);

	int status = 0;
//	status |= assertRelAbsEquals(cel1, c1, tolerance);
//	status |= assertRelAbsEquals(cel2, c2, tolerance);
	return status;
}

int main(int argc, char** argv) {

	int status = 0;

	status |= test_abscab();

	if (status != 0) {
		printf("%s: some test(s) failed :-(\n", argv[0]);
	} else {
		printf("%s: all test(s) passed :-)\n", argv[0]);
	}
}
