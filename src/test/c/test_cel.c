#include "util.h"

#include "cel.h"

/** test case for cel() implementation as described in section 4.2 of the 1969 Bulirsch article */
int test_cel() {
	double tolerance = 1.0e-15;

	double k_c = 0.1;
	double p1 =  4.1;
	double p2 = -4.1;
	double a = 1.2;
	double b = 1.1;

	double cel1 =  1.5464442694017956;
	double cel2 = -6.7687378198360556e-1;

	double c1 = cel(k_c, p1, a, b);
	double c2 = cel(k_c, p2, a, b);

//	double ra1 = fabs(cel1 - c1)/(1.0 + fabs(cel1));
//	double ra2 = fabs(cel2 - c2)/(1.0 + fabs(cel2));
//	printf("case 1: rel/abs deviation = %g\n", ra1);
//	printf("case 2: rel/abs deviation = %g\n", ra2);

	int status = 0;
	status |= assertRelAbsEquals(cel1, c1, tolerance);
	status |= assertRelAbsEquals(cel2, c2, tolerance);
	return status;
}

int main(int argc, char **argv) {

	int status = 0;

	status |= test_cel();

	if (status != 0) {
		printf("%s: some test(s) failed :-(\n", argv[0]);
	} else {
		printf("%s: all test(s) passed :-)\n", argv[0]);
	}
}
