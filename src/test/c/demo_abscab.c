
#include <stdio.h>

#include "util.h"
#include "abscab.h"









int main(int argc, char **argv) {

	int status = 0;

//	status |= testStraightWireSegment();

	if (status != 0) {
		printf("%s: some test(s) failed :-(\n", argv[0]);
	} else {
		printf("%s: all test(s) passed :-)\n", argv[0]);
	}
}

