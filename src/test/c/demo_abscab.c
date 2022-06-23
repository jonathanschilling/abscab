
#include <stdio.h>
#include <stdbool.h>

#ifdef _OPENMP
#include <omp.h>
#else // _OPENMP
#define omp_get_max_threads() 1
#endif // _OPENMP

#include "util.h"
#include "abscab.h"

double omega = 0.0;

double radius = 0.0;
void vertexSupplierStd(int idxVertex, double *point) {
	double phi = omega * idxVertex;
	point[0] = radius * cos(phi);
	point[1] = radius * sin(phi);
	point[2] = 0.0;
}

double rCorr = 0.0;
void vertexSupplierMcG(int idxVertex, double *point) {
	double phi = omega * idxVertex;
	point[0] = rCorr * cos(phi);
	point[1] = rCorr * sin(phi);
	point[2] = 0.0;
}

// GCC 12.1.0 on HP Z820, Arch Linux, 2022-06-23
// no compiler optimization flags specified:
//  1 thread : 9 min, 21 s (reference)
// 16 threads: 0 min, 38 s --> 15x speedup
// 32 threads: 0 min, 26 s --> 22x speedup
// -O3 (full optimizations, no loss is accuracy)
//  1 thread : 5 min, 56 s --> 60% speedup
// 32 threads: 0 min, 18 s --> 31x speedup on 16 physical cores with HT !!!
// -Ofast (at the expense of accuracy, but still not extremely off: 10 correct digits for binary64)
//  1 thread : 4 min, 28 s -->  2x speedup
// 32 threads: 0 min, 13 s --> 43x speedup
void demoMcGreivy() {

	radius = 1.23; // m
	double current = 17.0; // A

	double center[] = { 0.0, 0.0, 0.0 };
	double normal[] = { 0.0, 0.0, 1.0 };

	double evalPos[] = {
			10.0, 5.0, 0.0
	};

	double magneticField[3];
	magneticFieldCircularFilament(center, normal, radius, current, 1, evalPos, magneticField);
	double bZRef = magneticField[2];
	printf("ref B_z = %.3e\n", bZRef);

	// mimic circular wire loop as:
	// a) Polygon with points on the circule to be mimiced
	// b) Polygon with points slightly offset radially outward (McGreivy correction)
	// --> a) should have 2nd-order convergence;
	//     b) should have 4th-order convergence wrt. number of Polygon points

	int allNumPhi[] = {
			10, 30, 100, 300, 1000, 3000,
			10000, 30000, 100000, 300000,
			1000000, 3000000,
			10000000, 30000000,
			100000000, 300000000, 1000000000
	};
//		int[] allNumPhi = {10, 30, 100, 300, 1000, 3000};
	int numCases = 17;

	double *allBzStdErr = (double *) malloc(numCases * sizeof(double));
	double *allBzMcGErr = (double *) malloc(numCases * sizeof(double));

	double *resultTable = (double *) malloc(3 * numCases * sizeof(double));

	int numProcessors = omp_get_max_threads();
	printf("number of processors: %d\n", numProcessors);

	bool useCompensatedSummation = true;

	for (int i = 0; i < numCases; ++i) {

		int numPhi = allNumPhi[i];
		printf("case %2d/%2d: numPhi = %d\n", i+1, numCases, numPhi);

		omega = 2.0 * M_PI / (numPhi-1);

		memset(magneticField, 0, 3*sizeof(double));
		magneticFieldPolygonFilamentVertexSupplier_specPar_specSum(numPhi, vertexSupplierStd, current, 1, evalPos, magneticField, numProcessors, useCompensatedSummation);
		double bZStd = magneticField[2];

//			double[][] verticesStd = polygonCircleAround0(radius, numPhi);
//			double bZStd = ABSCAB.magneticFieldPolygonFilament(verticesStd, current, evalPos, numProcessors, useCompensatedSummation)[2][0];

		allBzStdErr[i] = errorMetric(bZRef, bZStd);
		printf("ABSCAB B_z = %.3e (err %g)\n", bZStd, allBzStdErr[i]);

		// McGreivy radius correction
		double dPhi = 2.0 * M_PI / (numPhi - 1); // spacing between points

		// TODO: understand derivation of alpha for special case of closed circle
		// |dr/ds| = 2*pi
		// --> alpha = 1/R * (dr)^2 / 12
		// == 4 pi^2 / (12 R)
		rCorr = radius * (1.0 + dPhi*dPhi/ 12);

		memset(magneticField, 0, 3*sizeof(double));
		magneticFieldPolygonFilamentVertexSupplier_specPar_specSum(numPhi, vertexSupplierMcG, current, 1, evalPos, magneticField, numProcessors, useCompensatedSummation);
		double bZMcG = magneticField[2];

//			double[][] verticesMcG = polygonCircleAround0(rCorr, numPhi);
//			double bZMcG = ABSCAB.magneticFieldPolygonFilament(verticesMcG, current, evalPos, numProcessors, useCompensatedSummation)[2][0];

		allBzMcGErr[i] = errorMetric(bZRef, bZMcG);
		printf("McGrvy B_z = %.3e (err %g)\n", bZMcG, allBzMcGErr[i]);

		resultTable[3 * i + 0] = numPhi;
		resultTable[3 * i + 1] = allBzStdErr[i];
		resultTable[3 * i + 2] = allBzMcGErr[i];
	}

	if (useCompensatedSummation) {
		dumpToFile(3, numCases, resultTable, "convergenceMcGreivy_CompensatedSummation.dat");
	} else {
		dumpToFile(3, numCases, resultTable, "convergenceMcGreivy_StandardSummation.dat");
	}

	free(allBzStdErr);
	free(allBzMcGErr);
	free(resultTable);
}

//static double[][] polygonCircleAround0(double radius, int numPhi) {
//	double[][] ret = new double[3][numPhi];
//	double omega = 2.0*Math.PI / (numPhi-1);
//	for (int i=0; i<numPhi; ++i) {
//		double phi = omega * i;
//		ret[0][i] = radius * Math.cos(phi);
//		ret[1][i] = radius * Math.sin(phi);
//	}
//
//	return ret;
//}

int main(int argc, char **argv) {


	demoMcGreivy();

	return 0;
}

