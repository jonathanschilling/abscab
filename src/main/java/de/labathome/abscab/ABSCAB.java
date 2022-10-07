package de.labathome.abscab;

import java.util.Locale;
import java.util.Objects;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.function.IntFunction;

import org.apache.commons.math3.util.FastMath;

/** Accurate Biot-Savart routines with Correct Asymptotic Behavior */
public class ABSCAB {

	/** vacuum magnetic permeability in Vs/Am (CODATA-2018) */
	public static final double MU_0 = 1.256_637_062_12e-6;

	/** vacuum magnetic permeability, divided by pi */
	private static final double MU_0_BY_PI = MU_0 / Math.PI;

	/** vacuum magnetic permeability, divided by 2 pi */
	private static final double MU_0_BY_2_PI = MU_0 / (2.0 * Math.PI);

	/** vacuum magnetic permeability, divided by 4 pi */
	private static final double MU_0_BY_4_PI = MU_0 / (4.0 * Math.PI);

	// --------------------------------------------------

	/**
	 * Compute the magnetic vector potential of a polygon filament
	 * at a number of evaluation locations.
	 * The computation is parallelized over all available processors.
	 * Kahan-Babuska compensated summation is used to compute the superposition
	 * of the contributions from the polygon vertices.
	 *
	 * @param vertices [3: x, y, z][numVertices] points along polygon; in m
	 * @param current current along polygon; in A
	 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
	 * @return [3: x, y, z][numEvalPos] magnetic vector potential at evaluation locations; in Tm
	 */
	public static double[][] vectorPotentialPolygonFilament(
			double[][] vertices,
			double current,
			double[][] evalPos) {
		int numProcessors = Runtime.getRuntime().availableProcessors();
		return vectorPotentialPolygonFilament(vertices, current, evalPos, numProcessors);
	}

	/**
	 * Compute the magnetic vector potential of a polygon filament
	 * at a number of evaluation locations.
	 * The computation is parallelized over all available processors.
	 * Kahan-Babuska compensated summation is used to compute the superposition
	 * of the contributions from the polygon vertices.
	 *
	 * @param numVertices number of points along polygon
	 * @param vertexSupplier callback which provides the points along the polygon as {x, y, z}; in m
	 * @param current current along polygon; in A
	 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
	 * @return [3: x, y, z][numEvalPos] magnetic vector potential at evaluation locations; in Tm
	 */
	public static double[][] vectorPotentialPolygonFilament(
			int numVertices,
			IntFunction<double[]> vertexSupplier,
			double current,
			double[][] evalPos) {
		int numProcessors = Runtime.getRuntime().availableProcessors();
		return vectorPotentialPolygonFilament(numVertices, vertexSupplier, current, evalPos, numProcessors);
	}

	/**
	 * Compute the magnetic field of a polygon filament
	 * at a number of evaluation locations.
	 * The computation is parallelized over all available processors.
	 * Kahan-Babuska compensated summation is used to compute the superposition
	 * of the contributions from the polygon vertices.
	 *
	 * @param vertices [3: x, y, z][numVertices] points along polygon; in m
	 * @param current current along polygon; in A
	 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
	 * @return [3: x, y, z][numEvalPos] magnetic field at evaluation locations; in T
	 */
	public static double[][] magneticFieldPolygonFilament(
			double[][] vertices,
			double current,
			double[][] evalPos) {
		int numProcessors = Runtime.getRuntime().availableProcessors();
		return magneticFieldPolygonFilament(vertices, current, evalPos, numProcessors);
	}

	/**
	 * Compute the magnetic field of a polygon filament
	 * at a number of evaluation locations.
	 * The computation is parallelized over all available processors.
	 * Kahan-Babuska compensated summation is used to compute the superposition
	 * of the contributions from the polygon vertices.
	 *
	 * @param numVertices number of points along polygon
	 * @param vertexSupplier callback which provides the points along the polygon as {x, y, z}; in m
	 * @param current current along polygon; in A
	 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
	 * @return [3: x, y, z][numEvalPos] magnetic field at evaluation locations; in T
	 */
	public static double[][] magneticFieldPolygonFilament(
			int numVertices,
			IntFunction<double[]> vertexSupplier,
			double current,
			double[][] evalPos) {
		int numProcessors = Runtime.getRuntime().availableProcessors();
		return magneticFieldPolygonFilament(numVertices, vertexSupplier, current, evalPos, numProcessors);
	}

	// --------------------------------------------------

	/**
	 * Compute the magnetic vector potential of a polygon filament
	 * at a number of evaluation locations.
	 * Kahan-Babuska compensated summation is used to compute the superposition
	 * of the contributions from the polygon vertices.
	 *
	 * @param vertices [3: x, y, z][numVertices] points along polygon; in m
	 * @param current current along polygon; in A
	 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
	 * @param numProcessors number of processors to use for parallelization
	 * @return [3: x, y, z][numEvalPos] magnetic vector potential at evaluation locations; in Tm
	 */
	public static double[][] vectorPotentialPolygonFilament(
			double[][] vertices,
			double current,
			double[][] evalPos,
			int numProcessors) {
		boolean useCompensatedSummation = true;
		return vectorPotentialPolygonFilament(vertices, current, evalPos, numProcessors, useCompensatedSummation);
	}

	/**
	 * Compute the magnetic vector potential of a polygon filament
	 * at a number of evaluation locations.
	 * Kahan-Babuska compensated summation is used to compute the superposition
	 * of the contributions from the polygon vertices.
	 *
	 * @param numVertices number of points along polygon
	 * @param vertexSupplier callback which provides the points along the polygon as {x, y, z}; in m
	 * @param current current along polygon; in A
	 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
	 * @param numProcessors number of processors to use for parallelization
	 * @return [3: x, y, z][numEvalPos] magnetic vector potential at evaluation locations; in Tm
	 */
	public static double[][] vectorPotentialPolygonFilament(
			int numVertices,
			IntFunction<double[]> vertexSupplier,
			double current,
			double[][] evalPos,
			int numProcessors) {
		boolean useCompensatedSummation = true;
		return vectorPotentialPolygonFilament(numVertices, vertexSupplier, current, evalPos, numProcessors, useCompensatedSummation);
	}

	/**
	 * Compute the magnetic field of a polygon filament
	 * at a number of evaluation locations.
	 * Kahan-Babuska compensated summation is used to compute the superposition
	 * of the contributions from the polygon vertices.
	 *
	 * @param vertices [3: x, y, z][numVertices] points along polygon; in m
	 * @param current current along polygon; in A
	 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
	 * @param numProcessors number of processors to use for parallelization
	 * @return [3: x, y, z][numEvalPos] magnetic field at evaluation locations; in T
	 */
	public static double[][] magneticFieldPolygonFilament(
			double[][] vertices,
			double current,
			double[][] evalPos,
			int numProcessors) {
		boolean useCompensatedSummation = true;
		return magneticFieldPolygonFilament(vertices, current, evalPos, numProcessors, useCompensatedSummation);
	}

	/**
	 * Compute the magnetic field of a polygon filament
	 * at a number of evaluation locations.
	 * Kahan-Babuska compensated summation is used to compute the superposition
	 * of the contributions from the polygon vertices.
	 *
	 * @param numVertices number of points along polygon
	 * @param vertexSupplier callback which provides the points along the polygon as {x, y, z}; in m
	 * @param current current along polygon; in A
	 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
	 * @param numProcessors number of processors to use for parallelization
	 * @return [3: x, y, z][numEvalPos] magnetic field at evaluation locations; in T
	 */
	public static double[][] magneticFieldPolygonFilament(
			int numVertices,
			IntFunction<double[]> vertexSupplier,
			double current,
			double[][] evalPos,
			int numProcessors) {
		boolean useCompensatedSummation = true;
		return magneticFieldPolygonFilament(numVertices, vertexSupplier, current, evalPos, numProcessors, useCompensatedSummation);
	}

	// --------------------------------------------------

	/**
	 * Compute the magnetic vector potential of a polygon filament
	 * at a number of evaluation locations.
	 *
	 * @param vertices [3: x, y, z][numVertices] points along polygon; in m
	 * @param current current along polygon; in A
	 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
	 * @param numProcessors number of processors to use for parallelization
	 * @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
	 *                                of the contributions from the polygon vertices; otherwise, use standard += summation
	 * @return [3: x, y, z][numEvalPos] magnetic vector potential at evaluation locations; in Tm
	 */
	public static double[][] vectorPotentialPolygonFilament(
			double[][] vertices,
			double current,
			double[][] evalPos,
			int numProcessors,
			boolean useCompensatedSummation) {
		final int numEvalPos = validateCartesianVectorInput(evalPos);
		double[][] vectorPotential = new double[3][numEvalPos];
		vectorPotentialPolygonFilament(vertices, current, evalPos, vectorPotential, numProcessors, useCompensatedSummation);
		return vectorPotential;
	}

	/**
	 * Compute the magnetic vector potential of a polygon filament
	 * at a number of evaluation locations.
	 *
	 * @param numVertices number of points along polygon
	 * @param vertexSupplier callback which provides the points along the polygon as {x, y, z}; in m
	 * @param current current along polygon; in A
	 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
	 * @param numProcessors number of processors to use for parallelization
	 * @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
	 *                                of the contributions from the polygon vertices; otherwise, use standard += summation
	 * @return [3: x, y, z][numEvalPos] magnetic vector potential at evaluation locations; in Tm
	 */
	public static double[][] vectorPotentialPolygonFilament(
			int numVertices,
			IntFunction<double[]> vertexSupplier,
			double current,
			double[][] evalPos,
			int numProcessors,
			boolean useCompensatedSummation) {
		final int numEvalPos = validateCartesianVectorInput(evalPos);
		double[][] vectorPotential = new double[3][numEvalPos];
		vectorPotentialPolygonFilament(numVertices, vertexSupplier, current, evalPos, vectorPotential, numProcessors, useCompensatedSummation);
		return vectorPotential;
	}

	/**
	 * Compute the magnetic field of a polygon filament
	 * at a number of evaluation locations.
	 *
	 * @param vertices [3: x, y, z][numVertices] points along polygon; in m
	 * @param current current along polygon; in A
	 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
	 * @param numProcessors number of processors to use for parallelization
	 * @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
	 *                                of the contributions from the polygon vertices; otherwise, use standard += summation
	 * @return [3: x, y, z][numEvalPos] magnetic field at evaluation locations; in T
	 */
	public static double[][] magneticFieldPolygonFilament(
			double[][] vertices,
			double current,
			double[][] evalPos,
			int numProcessors,
			boolean useCompensatedSummation) {
		final int numEvalPos = validateCartesianVectorInput(evalPos);
		double[][] magneticField = new double[3][numEvalPos];
		magneticFieldPolygonFilament(vertices, current, evalPos, magneticField, numProcessors, useCompensatedSummation);
		return magneticField;
	}

	/**
	 * Compute the magnetic field of a polygon filament
	 * at a number of evaluation locations.
	 *
	 * @param numVertices number of points along polygon
	 * @param vertexSupplier callback which provides the points along the polygon as {x, y, z}; in m
	 * @param current current along polygon; in A
	 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
	 * @param numProcessors number of processors to use for parallelization
	 * @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
	 *                                of the contributions from the polygon vertices; otherwise, use standard += summation
	 * @return [3: x, y, z][numEvalPos] magnetic field at evaluation locations; in T
	 */
	public static double[][] magneticFieldPolygonFilament(
			int numVertices,
			IntFunction<double[]> vertexSupplier,
			double current,
			double[][] evalPos,
			int numProcessors,
			boolean useCompensatedSummation) {
		final int numEvalPos = validateCartesianVectorInput(evalPos);
		double[][] magneticField = new double[3][numEvalPos];
		magneticFieldPolygonFilament(numVertices, vertexSupplier, current, evalPos, magneticField, numProcessors, useCompensatedSummation);
		return magneticField;
	}

	// --------------------------------------------------

	/**
	 * Compute the magnetic vector potential of a polygon filament
	 * at a number of evaluation locations.
	 *
	 * @param vertices [3: x, y, z][numVertices] points along polygon; in m
	 * @param current current along polygon; in A
	 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
	 * @param vectorPotential [3: x, y, z][numEvalPos] target array for magnetic vector potential at evaluation locations; in Tm
	 * @param numProcessors number of processors to use for parallelization
	 * @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
	 *                                of the contributions from the polygon vertices; otherwise, use standard += summation
	 */
	public static void vectorPotentialPolygonFilament(
			double[][] vertices,
			double current,
			double[][] evalPos,
			double[][] vectorPotential,
			int numProcessors,
			boolean useCompensatedSummation) {

		final int numVertices = validateCartesianVectorInput(vertices);

		if (numVertices < 2) {
			throw new RuntimeException("need at least 2 vertices, but only got " + numVertices);
		}

		final int numEvalPos = validateCartesianVectorInput(evalPos);

		if (numProcessors < 1) {
			throw new RuntimeException("need at least 1 processor, but only got " + numProcessors);
		}

		if (current == 0.0) {
			// TODO: Arrays.fill(vectorPotential, 0.0)
			return;
		}

		if (numProcessors == 1) {
			// single-threaded call
			final int idxSourceStart = 0;
			final int idxSourceEnd   = numVertices-1;
			final int idxEvalStart   = 0;
			final int idxEvalEnd     = numEvalPos;
			kernelVectorPotentialPolygonFilament(
					vertices, current,
					evalPos,
					vectorPotential,
					idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
					useCompensatedSummation);
		} else {
			// use multithreading

			if (numVertices-1 > numEvalPos) {
				// parallelize over nSource-1

				// Note that each thread needs its own copy of the vectorPotential array,
				// so this approach might need quite some memory in case the number of
				// threads and the number of evaluation points is large.

				final int nThreads;
				final int nSourcePerThread;
				if (numVertices-1 < numProcessors) {
					nThreads = numVertices-1;
					nSourcePerThread = 1;
				} else {
					nThreads = numProcessors;

					// It is better that many threads do more
					// than one thread needs to do more.
					nSourcePerThread = (int) Math.ceil( (numVertices-1.0) / nThreads);
				}

				final double[][][] vectorPotentialContributions = new double[nThreads][3][numEvalPos];

				ExecutorService service = Executors.newFixedThreadPool(nThreads);

				// submit jobs
				for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
					service.submit(new Runnable() {
						private int idxThread;

						public Runnable init(final int idxThread) {
							this.idxThread = idxThread;
							return this;
						}

						@Override
						public void run() {
							try {
								final int idxSourceStart =           idxThread    * nSourcePerThread;
								final int idxSourceEnd   = Math.min((idxThread+1) * nSourcePerThread, numVertices-1);
								final int idxEvalStart   = 0;
								final int idxEvalEnd     = numEvalPos;

								kernelVectorPotentialPolygonFilament(
										vertices, current,
										evalPos,
										vectorPotentialContributions[idxThread],
										idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
										useCompensatedSummation);

							} catch (Exception e) {
								e.printStackTrace();
							}
						}
					}.init(idxThread));
				}

				// accept no more new threads and start execution
				service.shutdown();

				// wait until all threads are done
				try {
					service.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}

				// sum up contributions from source chunks
				if (useCompensatedSummation) {
					for (int i=0; i<numEvalPos; ++i) {
						CompensatedSummation sumX = new CompensatedSummation();
						CompensatedSummation sumY = new CompensatedSummation();
						CompensatedSummation sumZ = new CompensatedSummation();
						for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
							sumX.add(vectorPotentialContributions[idxThread][0][i]);
							sumY.add(vectorPotentialContributions[idxThread][1][i]);
							sumZ.add(vectorPotentialContributions[idxThread][2][i]);
						}
						vectorPotential[0][i] = sumX.getSum();
						vectorPotential[1][i] = sumY.getSum();
						vectorPotential[2][i] = sumZ.getSum();
					}
				} else {
					for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
						for (int i=0; i<numEvalPos; ++i) {
							vectorPotential[0][i] += vectorPotentialContributions[idxThread][0][i];
							vectorPotential[1][i] += vectorPotentialContributions[idxThread][1][i];
							vectorPotential[2][i] += vectorPotentialContributions[idxThread][2][i];
						}
					}
				}
			} else { // nEval > nSource
				// parallelize over nEval

				final int nThreads;
				final int nEvalPerThread;
				if (numEvalPos < numProcessors) {
					nThreads = numEvalPos;
					nEvalPerThread = 1;
				} else {
					nThreads = numProcessors;
					nEvalPerThread = (int) Math.ceil( ((double) numEvalPos) / nThreads );
				}

				ExecutorService service = Executors.newFixedThreadPool(nThreads);

				// submit jobs
				for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
					service.submit(new Runnable() {
						private int idxThread;

						public Runnable init(final int idxThread) {
							this.idxThread = idxThread;
							return this;
						}

						@Override
						public void run() {
							try {
								final int idxSourceStart = 0;
								final int idxSourceEnd   = numVertices-1;
								final int idxEvalStart   =           idxThread    * nEvalPerThread;
								final int idxEvalEnd     = Math.min((idxThread+1) * nEvalPerThread, numEvalPos);

								kernelVectorPotentialPolygonFilament(
										vertices, current,
										evalPos,
										vectorPotential,
										idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
										useCompensatedSummation);

							} catch (Exception e) {
								e.printStackTrace();
							}
						}
					}.init(idxThread));
				}

				// accept no more new threads and start execution
				service.shutdown();

				// wait until all threads are done
				try {
					service.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			} // parallelize over nSource or nEval
		} // parallelization
	}

	/**
	 * Compute the magnetic vector potential of a polygon filament
	 * at a number of evaluation locations.
	 *
	 * @param numVertices number of points along polygon
	 * @param vertexSupplier callback which provides the points along the polygon as {x, y, z}; in m
	 * @param current current along polygon; in A
	 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
	 * @param vectorPotential [3: x, y, z][numEvalPos] target array for magnetic vector potential at evaluation locations; in Tm
	 * @param numProcessors number of processors to use for parallelization
	 * @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
	 *                                of the contributions from the polygon vertices; otherwise, use standard += summation
	 */
	public static void vectorPotentialPolygonFilament(
			int numVertices,
			IntFunction<double[]> vertexSupplier,
			double current,
			double[][] evalPos,
			double[][] vectorPotential,
			int numProcessors,
			boolean useCompensatedSummation) {

		if (numVertices < 2) {
			throw new RuntimeException("need at least 2 vertices, but only got " + numVertices);
		}

		final int numEvalPos = validateCartesianVectorInput(evalPos);

		if (numProcessors < 1) {
			throw new RuntimeException("need at least 1 processor, but only got " + numProcessors);
		}

		if (current == 0.0) {
			// TODO: Arrays.fill(vectorPotential, 0.0)
			return;
		}

		if (numProcessors == 1) {
			// single-threaded call
			final int idxSourceStart = 0;
			final int idxSourceEnd   = numVertices-1;
			final int idxEvalStart   = 0;
			final int idxEvalEnd     = numEvalPos;
			kernelVectorPotentialPolygonFilament(
					vertexSupplier, current,
					evalPos,
					vectorPotential,
					idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
					useCompensatedSummation);
		} else {
			// use multithreading

			if (numVertices-1 > numEvalPos) {
				// parallelize over nSource-1

				// Note that each thread needs its own copy of the vectorPotential array,
				// so this approach might need quite some memory in case the number of
				// threads and the number of evaluation points is large.

				final int nThreads;
				final int nSourcePerThread;
				if (numVertices-1 < numProcessors) {
					nThreads = numVertices-1;
					nSourcePerThread = 1;
				} else {
					nThreads = numProcessors;

					// It is better that many threads do more
					// than one thread needs to do more.
					nSourcePerThread = (int) Math.ceil( ((double)(numVertices - 1)) / nThreads);
				}

				final double[][][] vectorPotentialContributions = new double[nThreads][3][numEvalPos];

				ExecutorService service = Executors.newFixedThreadPool(nThreads);

				// submit jobs
				for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
					service.submit(new Runnable() {
						private int idxThread;

						public Runnable init(final int idxThread) {
							this.idxThread = idxThread;
							return this;
						}

						@Override
						public void run() {
							try {
								final int idxSourceStart =           idxThread    * nSourcePerThread;
								final int idxSourceEnd   = Math.min((idxThread+1) * nSourcePerThread, numVertices-1);
								final int idxEvalStart   = 0;
								final int idxEvalEnd     = numEvalPos;

								kernelVectorPotentialPolygonFilament(
										vertexSupplier, current,
										evalPos,
										vectorPotentialContributions[idxThread],
										idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
										useCompensatedSummation);

							} catch (Exception e) {
								e.printStackTrace();
							}
						}
					}.init(idxThread));
				}

				// accept no more new threads and start execution
				service.shutdown();

				// wait until all threads are done
				try {
					service.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}

				// sum up contributions from source chunks
				if (useCompensatedSummation) {
					for (int i=0; i<numEvalPos; ++i) {
						CompensatedSummation sumX = new CompensatedSummation();
						CompensatedSummation sumY = new CompensatedSummation();
						CompensatedSummation sumZ = new CompensatedSummation();
						for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
							sumX.add(vectorPotentialContributions[idxThread][0][i]);
							sumY.add(vectorPotentialContributions[idxThread][1][i]);
							sumZ.add(vectorPotentialContributions[idxThread][2][i]);
						}
						vectorPotential[0][i] = sumX.getSum();
						vectorPotential[1][i] = sumY.getSum();
						vectorPotential[2][i] = sumZ.getSum();
					}
				} else {
					for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
						for (int i=0; i<numEvalPos; ++i) {
							vectorPotential[0][i] += vectorPotentialContributions[idxThread][0][i];
							vectorPotential[1][i] += vectorPotentialContributions[idxThread][1][i];
							vectorPotential[2][i] += vectorPotentialContributions[idxThread][2][i];
						}
					}
				}
			} else { // nEval > nSource
				// parallelize over nEval

				final int nThreads;
				final int nEvalPerThread;
				if (numEvalPos < numProcessors) {
					nThreads = numEvalPos;
					nEvalPerThread = 1;
				} else {
					nThreads = numProcessors;

					// It is better that many threads do more
					// than one thread needs to do more.
					nEvalPerThread = (int) Math.ceil( ((double) numEvalPos) / nThreads );
				}

				ExecutorService service = Executors.newFixedThreadPool(nThreads);

				// submit jobs
				for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
					service.submit(new Runnable() {
						private int idxThread;

						public Runnable init(final int idxThread) {
							this.idxThread = idxThread;
							return this;
						}

						@Override
						public void run() {
							try {
								final int idxSourceStart = 0;
								final int idxSourceEnd   = numVertices-1;
								final int idxEvalStart   =           idxThread    * nEvalPerThread;
								final int idxEvalEnd     = Math.min((idxThread+1) * nEvalPerThread, numEvalPos);

								kernelVectorPotentialPolygonFilament(
										vertexSupplier, current,
										evalPos,
										vectorPotential,
										idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
										useCompensatedSummation);

							} catch (Exception e) {
								e.printStackTrace();
							}
						}
					}.init(idxThread));
				}

				// accept no more new threads and start execution
				service.shutdown();

				// wait until all threads are done
				try {
					service.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			} // parallelize over nSource or nEval
		} // parallelization
	}

	/**
	 * Compute the magnetic field of a polygon filament
	 * at a number of evaluation locations.
	 *
	 * @param vertices [3: x, y, z][numVertices] points along polygon; in m
	 * @param current current along polygon; in A
	 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
	 * @param magneticField [3: x, y, z][numEvalPos] target array for magnetic field at evaluation locations; in T
	 * @param numProcessors number of processors to use for parallelization
	 * @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
	 *                                of the contributions from the polygon vertices; otherwise, use standard += summation
	 */
	public static void magneticFieldPolygonFilament(
			double[][] vertices,
			double current,
			double[][] evalPos,
			double[][] magneticField,
			int numProcessors,
			boolean useCompensatedSummation) {

		final int numVertices = validateCartesianVectorInput(vertices);

		if (numVertices < 2) {
			throw new RuntimeException("need at least 2 vertices, but only got " + numVertices);
		}

		final int numEvalPos = validateCartesianVectorInput(evalPos);

		if (numProcessors < 1) {
			throw new RuntimeException("need at least 1 processor, but only got " + numProcessors);
		}

		if (current == 0.0) {
			// TODO: Arrays.fill(magneticField, 0.0)
			return;
		}

		if (numProcessors == 1) {
			// single-threaded call
			final int idxSourceStart = 0;
			final int idxSourceEnd   = numVertices-1;
			final int idxEvalStart   = 0;
			final int idxEvalEnd     = numEvalPos;
			kernelMagneticFieldPolygonFilament(
					vertices, current,
					evalPos,
					magneticField,
					idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
					useCompensatedSummation);
		} else {
			// use multithreading

			if (numVertices-1 > numEvalPos) {
				// parallelize over nSource-1

				// Note that each thread needs its own copy of the vectorPotential array,
				// so this approach might need quite some memory in case the number of
				// threads and the number of evaluation points is large.

				final int nThreads;
				final int nSourcePerThread;
				if (numVertices-1 < numProcessors) {
					nThreads = numVertices-1;
					nSourcePerThread = 1;
				} else {
					nThreads = numProcessors;

					// It is better that many threads do more
					// than one thread needs to do more.
					nSourcePerThread = (int) Math.ceil( (numVertices-1.0) / nThreads);
				}

				final double[][][] magneticFieldContributions = new double[nThreads][3][numEvalPos];

				ExecutorService service = Executors.newFixedThreadPool(nThreads);

				// submit jobs
				for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
					service.submit(new Runnable() {
						private int idxThread;

						public Runnable init(final int idxThread) {
							this.idxThread = idxThread;
							return this;
						}

						@Override
						public void run() {
							try {
								final int idxSourceStart =           idxThread    * nSourcePerThread;
								final int idxSourceEnd   = Math.min((idxThread+1) * nSourcePerThread, numVertices-1);
								final int idxEvalStart   = 0;
								final int idxEvalEnd     = numEvalPos;

								kernelMagneticFieldPolygonFilament(
										vertices, current,
										evalPos,
										magneticFieldContributions[idxThread],
										idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
										useCompensatedSummation);

							} catch (Exception e) {
								e.printStackTrace();
							}
						}
					}.init(idxThread));
				}

				// accept no more new threads and start execution
				service.shutdown();

				// wait until all threads are done
				try {
					service.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}

				// sum up contributions from source chunks
				if (useCompensatedSummation) {
					for (int i=0; i<numEvalPos; ++i) {
						CompensatedSummation sumX = new CompensatedSummation();
						CompensatedSummation sumY = new CompensatedSummation();
						CompensatedSummation sumZ = new CompensatedSummation();
						for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
							sumX.add(magneticFieldContributions[idxThread][0][i]);
							sumY.add(magneticFieldContributions[idxThread][1][i]);
							sumZ.add(magneticFieldContributions[idxThread][2][i]);
						}
						magneticField[0][i] = sumX.getSum();
						magneticField[1][i] = sumY.getSum();
						magneticField[2][i] = sumZ.getSum();
					}
				} else {
					for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
						for (int i=0; i<numEvalPos; ++i) {
							magneticField[0][i] += magneticFieldContributions[idxThread][0][i];
							magneticField[1][i] += magneticFieldContributions[idxThread][1][i];
							magneticField[2][i] += magneticFieldContributions[idxThread][2][i];
						}
					}
				}
			} else { // nEval > nSource
				// parallelize over nEval

				final int nThreads;
				final int nEvalPerThread;
				if (numEvalPos < numProcessors) {
					nThreads = numEvalPos;
					nEvalPerThread = 1;
				} else {
					nThreads = numProcessors;
					nEvalPerThread = (int) Math.ceil( ((double) numEvalPos) / nThreads );
				}

				ExecutorService service = Executors.newFixedThreadPool(nThreads);

				// submit jobs
				for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
					service.submit(new Runnable() {
						private int idxThread;

						public Runnable init(final int idxThread) {
							this.idxThread = idxThread;
							return this;
						}

						@Override
						public void run() {
							try {
								final int idxSourceStart = 0;
								final int idxSourceEnd   = numVertices-1;
								final int idxEvalStart   =           idxThread    * nEvalPerThread;
								final int idxEvalEnd     = Math.min((idxThread+1) * nEvalPerThread, numEvalPos);

								kernelMagneticFieldPolygonFilament(
										vertices, current,
										evalPos,
										magneticField,
										idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
										useCompensatedSummation);

							} catch (Exception e) {
								e.printStackTrace();
							}
						}
					}.init(idxThread));
				}

				// accept no more new threads and start execution
				service.shutdown();

				// wait until all threads are done
				try {
					service.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			} // parallelize over nSource or nEval
		} // parallelization
	}

	/**
	 * Compute the magnetic field of a polygon filament
	 * at a number of evaluation locations.
	 *
	 * @param numVertices number of points along polygon
	 * @param vertexSupplier callback which provides the points along the polygon as {x, y, z}; in m
	 * @param current current along polygon; in A
	 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
	 * @param magneticField [3: x, y, z][numEvalPos] target array for magnetic field at evaluation locations; in T
	 * @param numProcessors number of processors to use for parallelization
	 * @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
	 *                                of the contributions from the polygon vertices; otherwise, use standard += summation
	 */
	public static void magneticFieldPolygonFilament(
			int numVertices,
			IntFunction<double[]> vertexSupplier,
			double current,
			double[][] evalPos,
			double[][] magneticField,
			int numProcessors,
			boolean useCompensatedSummation) {

		if (numVertices < 2) {
			throw new RuntimeException("need at least 2 vertices, but only got " + numVertices);
		}

		final int numEvalPos = validateCartesianVectorInput(evalPos);

		if (numProcessors < 1) {
			throw new RuntimeException("need at least 1 processor, but only got " + numProcessors);
		}

		if (current == 0.0) {
			// TODO: Arrays.fill(magneticField, 0.0)
			return;
		}

		if (numProcessors == 1) {
			// single-threaded call
			final int idxSourceStart = 0;
			final int idxSourceEnd   = numVertices-1;
			final int idxEvalStart   = 0;
			final int idxEvalEnd     = numEvalPos;
			kernelMagneticFieldPolygonFilament(
					vertexSupplier, current,
					evalPos,
					magneticField,
					idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
					useCompensatedSummation);
		} else {
			// use multithreading

			if (numVertices-1 > numEvalPos) {
				// parallelize over nSource-1

				// Note that each thread needs its own copy of the vectorPotential array,
				// so this approach might need quite some memory in case the number of
				// threads and the number of evaluation points is large.

				final int nThreads;
				final int nSourcePerThread;
				if (numVertices-1 < numProcessors) {
					nThreads = numVertices-1;
					nSourcePerThread = 1;
				} else {
					nThreads = numProcessors;

					// It is better that many threads do more
					// than one thread needs to do more.
					nSourcePerThread = (int) Math.ceil( ((double)(numVertices - 1)) / nThreads);
				}

				final double[][][] magneticFieldContributions = new double[nThreads][3][numEvalPos];

				ExecutorService service = Executors.newFixedThreadPool(nThreads);

				// submit jobs
				for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
					service.submit(new Runnable() {
						private int idxThread;

						public Runnable init(final int idxThread) {
							this.idxThread = idxThread;
							return this;
						}

						@Override
						public void run() {
							try {
								final int idxSourceStart =           idxThread    * nSourcePerThread;
								final int idxSourceEnd   = Math.min((idxThread+1) * nSourcePerThread, numVertices-1);
								final int idxEvalStart   = 0;
								final int idxEvalEnd     = numEvalPos;

								kernelMagneticFieldPolygonFilament(
										vertexSupplier, current,
										evalPos,
										magneticFieldContributions[idxThread],
										idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
										useCompensatedSummation);

							} catch (Exception e) {
								e.printStackTrace();
							}
						}
					}.init(idxThread));
				}

				// accept no more new threads and start execution
				service.shutdown();

				// wait until all threads are done
				try {
					service.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}

				// sum up contributions from source chunks
				if (useCompensatedSummation) {
					for (int i=0; i<numEvalPos; ++i) {
						CompensatedSummation sumX = new CompensatedSummation();
						CompensatedSummation sumY = new CompensatedSummation();
						CompensatedSummation sumZ = new CompensatedSummation();
						for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
							sumX.add(magneticFieldContributions[idxThread][0][i]);
							sumY.add(magneticFieldContributions[idxThread][1][i]);
							sumZ.add(magneticFieldContributions[idxThread][2][i]);
						}
						magneticField[0][i] = sumX.getSum();
						magneticField[1][i] = sumY.getSum();
						magneticField[2][i] = sumZ.getSum();
					}
				} else {
					for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
						for (int i=0; i<numEvalPos; ++i) {
							magneticField[0][i] += magneticFieldContributions[idxThread][0][i];
							magneticField[1][i] += magneticFieldContributions[idxThread][1][i];
							magneticField[2][i] += magneticFieldContributions[idxThread][2][i];
						}
					}
				}
			} else { // nEval > nSource
				// parallelize over nEval

				final int nThreads;
				final int nEvalPerThread;
				if (numEvalPos < numProcessors) {
					nThreads = numEvalPos;
					nEvalPerThread = 1;
				} else {
					nThreads = numProcessors;

					// It is better that many threads do more
					// than one thread needs to do more.
					nEvalPerThread = (int) Math.ceil( ((double) numEvalPos) / nThreads );
				}

				ExecutorService service = Executors.newFixedThreadPool(nThreads);

				// submit jobs
				for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
					service.submit(new Runnable() {
						private int idxThread;

						public Runnable init(final int idxThread) {
							this.idxThread = idxThread;
							return this;
						}

						@Override
						public void run() {
							try {
								final int idxSourceStart = 0;
								final int idxSourceEnd   = numVertices-1;
								final int idxEvalStart   =           idxThread    * nEvalPerThread;
								final int idxEvalEnd     = Math.min((idxThread+1) * nEvalPerThread, numEvalPos);

								kernelMagneticFieldPolygonFilament(
										vertexSupplier, current,
										evalPos,
										magneticField,
										idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
										useCompensatedSummation);

							} catch (Exception e) {
								e.printStackTrace();
							}
						}
					}.init(idxThread));
				}

				// accept no more new threads and start execution
				service.shutdown();

				// wait until all threads are done
				try {
					service.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			} // parallelize over nSource or nEval
		} // parallelization
	}

	// --------------------------------------------------

	/**
	 * Compute the magnetic vector potential of a polygon filament
	 * at a number of evaluation locations.
	 *
	 * @param vertices [3: x, y, z][numVertices] points along polygon; in m
	 * @param current current along polygon; in A
	 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
	 * @param vectorPotential [3: x, y, z][numEvalPos] target array for magnetic vector potential at evaluation locations; in Tm
	 * @param idxSourceStart first index in {@code vertices} to take into account
	 * @param idxSourceEnd (last+1) index in {@code vertices} to take into account
	 * @param idxEvalStart first index in {@code evalPos} to take into account
	 * @param idxEvalEnd (last+1) index in {@code evalPos} to take into account
	 * @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
	 *                                of the contributions from the polygon vertices; otherwise, use standard += summation
	 */
	private static void kernelVectorPotentialPolygonFilament(
			double[][] vertices,
			double current,
			double[][] evalPos,
			double[][] vectorPotential,
			int idxSourceStart,
			int idxSourceEnd,
			int idxEvalStart,
			int idxEvalEnd,
			boolean useCompensatedSummation) {

		final double aPrefactor = MU_0_BY_2_PI * current;

		// setup compensated summation objects
		final CompensatedSummation[] aXSum;
		final CompensatedSummation[] aYSum;
		final CompensatedSummation[] aZSum;
		if (useCompensatedSummation) {
			int numEvalPos = idxEvalEnd - idxEvalStart;
			aXSum = new CompensatedSummation[numEvalPos];
			aYSum = new CompensatedSummation[numEvalPos];
			aZSum = new CompensatedSummation[numEvalPos];
			for (int idxEval = 0; idxEval < numEvalPos; ++idxEval) {
				aXSum[idxEval] = new CompensatedSummation();
				aYSum[idxEval] = new CompensatedSummation();
				aZSum[idxEval] = new CompensatedSummation();
			}
		} else {
			aXSum = null;
			aYSum = null;
			aZSum = null;
		}

		double x_i = vertices[0][idxSourceStart];
		double y_i = vertices[1][idxSourceStart];
		double z_i = vertices[2][idxSourceStart];

		for (int idxSource = idxSourceStart; idxSource < idxSourceEnd; ++idxSource) {

			final double x_f = vertices[0][idxSource+1];
			final double y_f = vertices[1][idxSource+1];
			final double z_f = vertices[2][idxSource+1];

			// vector from start to end of i:th wire segment
			final double dx = x_f - x_i;
			final double dy = y_f - y_i;
			final double dz = z_f - z_i;

			// squared length of wire segment
			final double l2 = dx * dx + dy * dy + dz * dz;
			if (l2 == 0.0) {
				// skip zero-length segments: no contribution
				continue;
			}

			// length of wire segment
			final double l = Math.sqrt(l2);

			// unit vector parallel to wire segment
			final double eX = dx / l;
			final double eY = dy / l;
			final double eZ = dz / l;

			for (int idxEval = idxEvalStart; idxEval < idxEvalEnd; ++idxEval) {

				// vector from start of wire segment to eval pos
				final double r0x = evalPos[0][idxEval] - x_i;
				final double r0y = evalPos[1][idxEval] - y_i;
				final double r0z = evalPos[2][idxEval] - z_i;

				// z position along axis of wire segment
				final double alignedZ = eX * r0x + eY * r0y + eZ * r0z;

				// normalized z component of evaluation location in coordinate system of wire segment
				final double zP = alignedZ / l;

				// vector perpendicular to axis of wire segment, pointing at evaluation pos
				final double rPerpX = r0x - alignedZ * eX;
				final double rPerpY = r0y - alignedZ * eY;
				final double rPerpZ = r0z - alignedZ * eZ;

				// perpendicular distance between evalPos and axis of wire segment
				final double alignedR = Math.sqrt(rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ);

				// normalized rho component of evaluation location in coordinate system of wire segment
				final double rhoP = alignedR / l;

				// compute parallel component of magnetic vector potential, including current and mu_0
				final double aParallel = aPrefactor * straightWireSegment_A_z(rhoP, zP);

				// add contribution from wire segment to result
				if (useCompensatedSummation) {
					aXSum[idxEval - idxEvalStart].add(aParallel * eX);
					aYSum[idxEval - idxEvalStart].add(aParallel * eY);
					aZSum[idxEval - idxEvalStart].add(aParallel * eZ);
				} else {
					vectorPotential[0][idxEval] += aParallel * eX;
					vectorPotential[1][idxEval] += aParallel * eY;
					vectorPotential[2][idxEval] += aParallel * eZ;
				}
			}

			// shift to next point
			x_i = x_f;
			y_i = y_f;
			z_i = z_f;
		}

		if (useCompensatedSummation) {
			// obtain compensated sums from summation objects
			for (int idxEval = idxEvalStart; idxEval < idxEvalEnd; ++idxEval) {
				vectorPotential[0][idxEval] = aXSum[idxEval - idxEvalStart].getSum();
				vectorPotential[1][idxEval] = aYSum[idxEval - idxEvalStart].getSum();
				vectorPotential[2][idxEval] = aZSum[idxEval - idxEvalStart].getSum();
			}
		}
	}

	/**
	 * Compute the magnetic vector potential of a polygon filament
	 * at a number of evaluation locations.
	 *
	 * @param vertexSupplier callback which provides the points along the polygon as {x, y, z}; in m
	 * @param current current along polygon; in A
	 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
	 * @param vectorPotential [3: x, y, z][numEvalPos] target array for magnetic vector potential at evaluation locations; in Tm
	 * @param idxSourceStart first index in {@code vertices} to take into account
	 * @param idxSourceEnd (last+1) index in {@code vertices} to take into account
	 * @param idxEvalStart first index in {@code evalPos} to take into account
	 * @param idxEvalEnd (last+1) index in {@code evalPos} to take into account
	 * @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
	 *                                of the contributions from the polygon vertices; otherwise, use standard += summation
	 */
	private static void kernelVectorPotentialPolygonFilament(
			IntFunction<double[]> vertexSupplier,
			double current,
			double[][] evalPos,
			double[][] vectorPotential,
			int idxSourceStart,
			int idxSourceEnd,
			int idxEvalStart,
			int idxEvalEnd,
			boolean useCompensatedSummation) {

		final double aPrefactor = MU_0_BY_2_PI * current;

		// setup compensated summation objects
		final CompensatedSummation[] aXSum;
		final CompensatedSummation[] aYSum;
		final CompensatedSummation[] aZSum;
		if (useCompensatedSummation) {
			int numEvalPos = idxEvalEnd - idxEvalStart;
			aXSum = new CompensatedSummation[numEvalPos];
			aYSum = new CompensatedSummation[numEvalPos];
			aZSum = new CompensatedSummation[numEvalPos];
			for (int idxEval = 0; idxEval < numEvalPos; ++idxEval) {
				aXSum[idxEval] = new CompensatedSummation();
				aYSum[idxEval] = new CompensatedSummation();
				aZSum[idxEval] = new CompensatedSummation();
			}
		} else {
			aXSum = null;
			aYSum = null;
			aZSum = null;
		}

		double[] firstPoint = vertexSupplier.apply(idxSourceStart);
		double x_i = firstPoint[0];
		double y_i = firstPoint[1];
		double z_i = firstPoint[2];

		for (int idxSource = idxSourceStart; idxSource < idxSourceEnd; ++idxSource) {

			final double[] nextPoint = vertexSupplier.apply(idxSource+1);
			final double x_f = nextPoint[0];
			final double y_f = nextPoint[1];
			final double z_f = nextPoint[2];

			// vector from start to end of i:th wire segment
			final double dx = x_f - x_i;
			final double dy = y_f - y_i;
			final double dz = z_f - z_i;

			// squared length of wire segment
			final double l2 = dx * dx + dy * dy + dz * dz;
			if (l2 == 0.0) {
				// skip zero-length segments: no contribution
				continue;
			}

			// length of wire segment
			final double l = Math.sqrt(l2);

			// unit vector parallel to wire segment
			final double eX = dx / l;
			final double eY = dy / l;
			final double eZ = dz / l;

			for (int idxEval = idxEvalStart; idxEval < idxEvalEnd; ++idxEval) {

				// vector from start of wire segment to eval pos
				final double r0x = evalPos[0][idxEval] - x_i;
				final double r0y = evalPos[1][idxEval] - y_i;
				final double r0z = evalPos[2][idxEval] - z_i;

				// z position along axis of wire segment
				final double alignedZ = eX * r0x + eY * r0y + eZ * r0z;

				// normalized z component of evaluation location in coordinate system of wire segment
				final double zP = alignedZ / l;

				// vector perpendicular to axis of wire segment, pointing at evaluation pos
				final double rPerpX = r0x - alignedZ * eX;
				final double rPerpY = r0y - alignedZ * eY;
				final double rPerpZ = r0z - alignedZ * eZ;

				// perpendicular distance between evalPos and axis of wire segment
				final double alignedR = Math.sqrt(rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ);

				// normalized rho component of evaluation location in coordinate system of wire segment
				final double rhoP = alignedR / l;

				// compute parallel component of magnetic vector potential, including current and mu_0
				final double aParallel = aPrefactor * straightWireSegment_A_z(rhoP, zP);

				// add contribution from wire segment to result
				if (useCompensatedSummation) {
					aXSum[idxEval - idxEvalStart].add(aParallel * eX);
					aYSum[idxEval - idxEvalStart].add(aParallel * eY);
					aZSum[idxEval - idxEvalStart].add(aParallel * eZ);
				} else {
					vectorPotential[0][idxEval] += aParallel * eX;
					vectorPotential[1][idxEval] += aParallel * eY;
					vectorPotential[2][idxEval] += aParallel * eZ;
				}
			}

			// shift to next point
			x_i = x_f;
			y_i = y_f;
			z_i = z_f;
		}

		if (useCompensatedSummation) {
			// obtain compensated sums from summation objects
			for (int idxEval = idxEvalStart; idxEval < idxEvalEnd; ++idxEval) {
				vectorPotential[0][idxEval] = aXSum[idxEval - idxEvalStart].getSum();
				vectorPotential[1][idxEval] = aYSum[idxEval - idxEvalStart].getSum();
				vectorPotential[2][idxEval] = aZSum[idxEval - idxEvalStart].getSum();
			}
		}
	}

	/**
	 * Compute the magnetic field of a polygon filament
	 * at a number of evaluation locations.
	 *
	 * @param vertices [3: x, y, z][numVertices] points along polygon; in m
	 * @param current current along polygon; in A
	 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
	 * @param magneticField [3: x, y, z][numEvalPos] target array for magnetic field at evaluation locations; in T
	 * @param idxSourceStart first index in {@code vertices} to take into account
	 * @param idxSourceEnd (last+1) index in {@code vertices} to take into account
	 * @param idxEvalStart first index in {@code evalPos} to take into account
	 * @param idxEvalEnd (last+1) index in {@code evalPos} to take into account
	 * @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
	 *                                of the contributions from the polygon vertices; otherwise, use standard += summation
	 */
	private static void kernelMagneticFieldPolygonFilament(
			double[][] vertices,
			double current,
			double[][] evalPos,
			double[][] magneticField,
			int idxSourceStart,
			int idxSourceEnd,
			int idxEvalStart,
			int idxEvalEnd,
			boolean useCompensatedSummation) {

		// needs additional division by length of wire segment!
		final double bPrefactorL = MU_0_BY_4_PI * current;

		// setup compensated summation objects
		final CompensatedSummation[] bXSum;
		final CompensatedSummation[] bYSum;
		final CompensatedSummation[] bZSum;
		if (useCompensatedSummation) {
			int numEvalPos = idxEvalEnd - idxEvalStart;
			bXSum = new CompensatedSummation[numEvalPos];
			bYSum = new CompensatedSummation[numEvalPos];
			bZSum = new CompensatedSummation[numEvalPos];
			for (int idxEval = 0; idxEval < numEvalPos; ++idxEval) {
				bXSum[idxEval] = new CompensatedSummation();
				bYSum[idxEval] = new CompensatedSummation();
				bZSum[idxEval] = new CompensatedSummation();
			}
		} else {
			bXSum = null;
			bYSum = null;
			bZSum = null;
		}

		double x_i = vertices[0][idxSourceStart];
		double y_i = vertices[1][idxSourceStart];
		double z_i = vertices[2][idxSourceStart];

		for (int idxSource = idxSourceStart; idxSource < idxSourceEnd; ++idxSource) {

			final double x_f = vertices[0][idxSource+1];
			final double y_f = vertices[1][idxSource+1];
			final double z_f = vertices[2][idxSource+1];

			// vector from start to end of i:th wire segment
			final double dx = x_f - x_i;
			final double dy = y_f - y_i;
			final double dz = z_f - z_i;

			// squared length of wire segment
			final double l2 = dx * dx + dy * dy + dz * dz;
			if (l2 == 0.0) {
				// skip zero-length segments: no contribution
				continue;
			}

			// length of wire segment
			final double l = Math.sqrt(l2);

			// assemble full prefactor for B_phi
			final double bPrefactor = bPrefactorL / l;

			// unit vector parallel to wire segment
			final double eX = dx / l;
			final double eY = dy / l;
			final double eZ = dz / l;

			for (int idxEval = idxEvalStart; idxEval < idxEvalEnd; ++idxEval) {

				// vector from start of wire segment to eval pos
				final double r0x = evalPos[0][idxEval] - x_i;
				final double r0y = evalPos[1][idxEval] - y_i;
				final double r0z = evalPos[2][idxEval] - z_i;

				// z position along axis of wire segment
				final double alignedZ = eX * r0x + eY * r0y + eZ * r0z;

				// normalized z component of evaluation location in coordinate system of wire segment
				final double zP = alignedZ / l;

				// vector perpendicular to axis of wire segment, pointing at evaluation pos
				final double rPerpX = r0x - alignedZ * eX;
				final double rPerpY = r0y - alignedZ * eY;
				final double rPerpZ = r0z - alignedZ * eZ;

				// perpendicular distance squared between evalPos and axis of wire segment
				final double alignedRSq = rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ;

				// B_phi is zero along axis of filament
				if (alignedRSq > 0.0) {

					// perpendicular distance between evalPos and axis of wire segment
					final double alignedR = Math.sqrt(alignedRSq) ;

					// normalized rho component of evaluation location in coordinate system of wire segment
					final double rhoP = alignedR / l;

					// compute tangential component of magnetic vector potential, including current and mu_0
					final double bPhi = bPrefactor * straightWireSegment_B_phi(rhoP, zP);

					// unit vector in radial direction
					final double eRX = rPerpX / alignedR;
					final double eRY = rPerpY / alignedR;
					final double eRZ = rPerpZ / alignedR;

					// compute cross product between e_z and e_rho to get e_phi
					final double ePhiX = eY * eRZ - eZ * eRY;
					final double ePhiY = eZ * eRX - eX * eRZ;
					final double ePhiZ = eX * eRY - eY * eRX;

					// add contribution from wire segment to result
					if (useCompensatedSummation) {
						bXSum[idxEval - idxEvalStart].add(bPhi * ePhiX);
						bYSum[idxEval - idxEvalStart].add(bPhi * ePhiY);
						bZSum[idxEval - idxEvalStart].add(bPhi * ePhiZ);
					} else {
						magneticField[0][idxEval] += bPhi * ePhiX;
						magneticField[1][idxEval] += bPhi * ePhiY;
						magneticField[2][idxEval] += bPhi * ePhiZ;
					}
				}
			}

			// shift to next point
			x_i = x_f;
			y_i = y_f;
			z_i = z_f;
		}

		if (useCompensatedSummation) {
			// obtain compensated sums from summation objects
			for (int idxEval = idxEvalStart; idxEval < idxEvalEnd; ++idxEval) {
				magneticField[0][idxEval] = bXSum[idxEval - idxEvalStart].getSum();
				magneticField[1][idxEval] = bYSum[idxEval - idxEvalStart].getSum();
				magneticField[2][idxEval] = bZSum[idxEval - idxEvalStart].getSum();
			}
		}
	}

	/**
	 * Compute the magnetic field of a polygon filament
	 * at a number of evaluation locations.
	 *
	 * @param vertexSupplier callback which provides the points along the polygon as {x, y, z}; in m
	 * @param current current along polygon; in A
	 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
	 * @param magneticField [3: x, y, z][numEvalPos] target array for magnetic field at evaluation locations; in T
	 * @param idxSourceStart first index in {@code vertices} to take into account
	 * @param idxSourceEnd (last+1) index in {@code vertices} to take into account
	 * @param idxEvalStart first index in {@code evalPos} to take into account
	 * @param idxEvalEnd (last+1) index in {@code evalPos} to take into account
	 * @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
	 *                                of the contributions from the polygon vertices; otherwise, use standard += summation
	 */
	private static void kernelMagneticFieldPolygonFilament(
			IntFunction<double[]> vertexSupplier,
			double current,
			double[][] evalPos,
			double[][] magneticField,
			int idxSourceStart,
			int idxSourceEnd,
			int idxEvalStart,
			int idxEvalEnd,
			boolean useCompensatedSummation) {

		// setup compensated summation objects
		final CompensatedSummation[] bXSum;
		final CompensatedSummation[] bYSum;
		final CompensatedSummation[] bZSum;
		if (useCompensatedSummation) {
			int numEvalPos = idxEvalEnd - idxEvalStart;
			bXSum = new CompensatedSummation[numEvalPos];
			bYSum = new CompensatedSummation[numEvalPos];
			bZSum = new CompensatedSummation[numEvalPos];
			for (int idxEval = 0; idxEval < numEvalPos; ++idxEval) {
				bXSum[idxEval] = new CompensatedSummation();
				bYSum[idxEval] = new CompensatedSummation();
				bZSum[idxEval] = new CompensatedSummation();
			}
		} else {
			bXSum = null;
			bYSum = null;
			bZSum = null;
		}

		// needs additional division by length of wire segment!
		final double bPrefactorL = MU_0_BY_4_PI * current;

		double[] firstPoint = vertexSupplier.apply(idxSourceStart);
		double x_i = firstPoint[0];
		double y_i = firstPoint[1];
		double z_i = firstPoint[2];

		for (int idxSource = idxSourceStart; idxSource < idxSourceEnd; ++idxSource) {

			final double[] nextPoint = vertexSupplier.apply(idxSource+1);
			final double x_f = nextPoint[0];
			final double y_f = nextPoint[1];
			final double z_f = nextPoint[2];

			// vector from start to end of i:th wire segment
			final double dx = x_f - x_i;
			final double dy = y_f - y_i;
			final double dz = z_f - z_i;

			// squared length of wire segment
			final double l2 = dx * dx + dy * dy + dz * dz;
			if (l2 == 0.0) {
				// skip zero-length segments: no contribution
				continue;
			}

			// length of wire segment
			final double l = Math.sqrt(l2);

			// assemble full prefactor for B_phi
			final double bPrefactor = bPrefactorL / l;

			// unit vector parallel to wire segment
			final double eX = dx / l;
			final double eY = dy / l;
			final double eZ = dz / l;

			for (int idxEval = idxEvalStart; idxEval < idxEvalEnd; ++idxEval) {

				// vector from start of wire segment to eval pos
				final double r0x = evalPos[0][idxEval] - x_i;
				final double r0y = evalPos[1][idxEval] - y_i;
				final double r0z = evalPos[2][idxEval] - z_i;

				// z position along axis of wire segment
				final double alignedZ = eX * r0x + eY * r0y + eZ * r0z;

				// normalized z component of evaluation location in coordinate system of wire segment
				final double zP = alignedZ / l;

				// vector perpendicular to axis of wire segment, pointing at evaluation pos
				final double rPerpX = r0x - alignedZ * eX;
				final double rPerpY = r0y - alignedZ * eY;
				final double rPerpZ = r0z - alignedZ * eZ;

				// perpendicular distance squared between evalPos and axis of wire segment
				final double alignedRSq = rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ;

				// B_phi is zero along axis of filament
				if (alignedRSq > 0.0) {

					// perpendicular distance between evalPos and axis of wire segment
					final double alignedR = Math.sqrt(alignedRSq);

					// normalized rho component of evaluation location in coordinate system of wire segment
					final double rhoP = alignedR / l;

					// compute tangential component of magnetic vector potential, including current and mu_0
					final double bPhi = bPrefactor * straightWireSegment_B_phi(rhoP, zP);

					// unit vector in radial direction
					final double eRX = rPerpX / alignedR;
					final double eRY = rPerpY / alignedR;
					final double eRZ = rPerpZ / alignedR;

					// compute cross product between e_z and e_rho to get e_phi
					final double ePhiX = eY * eRZ - eZ * eRY;
					final double ePhiY = eZ * eRX - eX * eRZ;
					final double ePhiZ = eX * eRY - eY * eRX;

					// add contribution from wire segment to result
					if (useCompensatedSummation) {
						bXSum[idxEval - idxEvalStart].add(bPhi * ePhiX);
						bYSum[idxEval - idxEvalStart].add(bPhi * ePhiY);
						bZSum[idxEval - idxEvalStart].add(bPhi * ePhiZ);
					} else {
						magneticField[0][idxEval] += bPhi * ePhiX;
						magneticField[1][idxEval] += bPhi * ePhiY;
						magneticField[2][idxEval] += bPhi * ePhiZ;
					}
				}
			}

			// shift to next point
			x_i = x_f;
			y_i = y_f;
			z_i = z_f;
		}

		if (useCompensatedSummation) {
			// obtain compensated sums from summation objects
			for (int idxEval = idxEvalStart; idxEval < idxEvalEnd; ++idxEval) {
				magneticField[0][idxEval] = bXSum[idxEval - idxEvalStart].getSum();
				magneticField[1][idxEval] = bYSum[idxEval - idxEvalStart].getSum();
				magneticField[2][idxEval] = bZSum[idxEval - idxEvalStart].getSum();
			}
		}
	}

	// --------------------------------------------------

	/**
	 * Compute the magnetic vector potential of a circular wire loop.
	 *
	 * @param center  [3: x, y, z] origin of loop (in meters)
	 * @param normal  [3: x, y, z] normal vector of loop (in meters); will be
	 *                normalized internally
	 * @param radius  radius of the wire loop (in meters)
	 * @param current loop current (in A)
	 * @param evalPos [3: x, y, z][nEvalPos] evaluation locations (in meters)
	 * @return [3: A_x, A_y, A_z][nEvalPos] Cartesian components of the magnetic
	 *         vector potential evaluated at the given locations (in Tm)
	 */
	public static double[][] vectorPotentialCircularFilament(double[] center, double[] normal, double radius,
			double current, double[][] evalPos) {

		Objects.requireNonNull(center);
		if (center.length != 3) {
			throw new RuntimeException("center needs to have 3 elements, but has " + center.length);
		}

		Objects.requireNonNull(normal);
		if (normal.length != 3) {
			throw new RuntimeException("normal needs to have 3 elements, but has " + normal.length);
		}

		if (!Double.isFinite(radius) || radius <= 0.0) {
			throw new RuntimeException("radius must be finite and positive, but is " + radius);
		}

		int nEvalPos = validateCartesianVectorInput(evalPos);

		final double aPrefactor = MU_0_BY_PI * current;

		// squared length of normal vector
		double nLen2 = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];

		if (nLen2 == 0.0) {
			throw new RuntimeException("length of normal vector must not be zero");
		}

		// length of normal vector
		double nLen = Math.sqrt(nLen2);

		// unit normal vector of wire loop
		double eX = normal[0] / nLen;
		double eY = normal[1] / nLen;
		double eZ = normal[2] / nLen;

		final double[][] ret = new double[3][nEvalPos];

		for (int idxEval=0; idxEval<nEvalPos; ++idxEval) {

			// vector from center of wire loop to eval pos
			final double r0x = evalPos[0][idxEval] - center[0];
			final double r0y = evalPos[1][idxEval] - center[1];
			final double r0z = evalPos[2][idxEval] - center[2];

			// z position along normal of wire loop
			final double alignedZ = eX * r0x + eY * r0y + eZ * r0z;

			// normalized z component of evaluation location in coordinate system of wire loop
			final double zP = alignedZ / radius;

			// r0 projected onto axis of wire loop
			final double rParallelX = alignedZ * eX;
			final double rParallelY = alignedZ * eY;
			final double rParallelZ = alignedZ * eZ;

			// vector perpendicular to axis of wire loop, pointing at evaluation pos
			final double rPerpX = r0x - rParallelX;
			final double rPerpY = r0y - rParallelY;
			final double rPerpZ = r0z - rParallelZ;

			// perpendicular distance squared between evalPos and axis of wire loop
			double alignedRSq = rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ;

			// prevent division-by-zero when computing radial unit vector
			// A_phi is zero anyway on-axis --> no contribution expected
			if (alignedRSq > 0.0) {

				// perpendicular distance between evalPos and axis of wire loop
				final double alignedR = Math.sqrt(alignedRSq);

				// unit vector in radial direction
				final double eRX = rPerpX / alignedR;
				final double eRY = rPerpY / alignedR;
				final double eRZ = rPerpZ / alignedR;

				// normalized rho component of evaluation location in coordinate system of wire loop
				final double rhoP = alignedR / radius;

				// compute tangential component of magnetic vector potential, including current and mu_0
				final double aPhi = aPrefactor * circularWireLoop_A_phi(rhoP, zP);

				// compute cross product between e_z and e_rho to get e_phi
				final double ePhiX = eRY * eZ - eRZ * eY;
				final double ePhiY = eRZ * eX - eRX * eZ;
				final double ePhiZ = eRX * eY - eRY * eX;

				// add contribution from wire loop to result
				ret[0][idxEval] = aPhi * ePhiX;
				ret[1][idxEval] = aPhi * ePhiY;
				ret[2][idxEval] = aPhi * ePhiZ;
			}
		}

		return ret;
	}

	/**
	 * Compute the magnetic field of a circular wire loop.
	 *
	 * @param center  [3: x, y, z] origin of loop (in meters)
	 * @param normal  [3: x, y, z] normal vector of loop (in meters); will be
	 *                normalized internally
	 * @param radius  radius of the wire loop (in meters)
	 * @param current loop current (in A)
	 * @param evalPos [3: x, y, z][nEvalPos] evaluation locations (in meters)
	 * @return [3: B_x, B_y, B_z][nEvalPos] Cartesian components of the magnetic
	 *         field evaluated at the given locations (in T)
	 */
	public static double[][] magneticFieldCircularFilament(double[] center, double[] normal, double radius,
			double current, double[][] evalPos) {

		Objects.requireNonNull(center);
		if (center.length != 3) {
			throw new RuntimeException("center needs to have 3 elements, but has " + center.length);
		}

		Objects.requireNonNull(normal);
		if (normal.length != 3) {
			throw new RuntimeException("normal needs to have 3 elements, but has " + normal.length);
		}

		if (!Double.isFinite(radius) || radius <= 0.0) {
			throw new RuntimeException("radius must be finite and positive, but is " + radius);
		}

		int nEvalPos = validateCartesianVectorInput(evalPos);

		final double bPrefactor = MU_0_BY_PI * current / radius;

		// squared length of normal vector
		double nLen2 = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];

		if (nLen2 == 0.0) {
			throw new RuntimeException("length of normal vector must not be zero");
		}

		// length of normal vector
		double nLen = Math.sqrt(nLen2);

		// unit normal vector of wire loop
		double eX = normal[0] / nLen;
		double eY = normal[1] / nLen;
		double eZ = normal[2] / nLen;

		final double[][] ret = new double[3][nEvalPos];

		for (int idxEval=0; idxEval<nEvalPos; ++idxEval) {

			// vector from center of wire loop to eval pos
			final double r0x = evalPos[0][idxEval] - center[0];
			final double r0y = evalPos[1][idxEval] - center[1];
			final double r0z = evalPos[2][idxEval] - center[2];

			// z position along normal of wire loop
			final double alignedZ = eX * r0x + eY * r0y + eZ * r0z;

			// normalized z component of evaluation location in coordinate system of wire loop
			final double zP = alignedZ / radius;

			// r0 projected onto axis of wire loop
			final double rParallelX = alignedZ * eX;
			final double rParallelY = alignedZ * eY;
			final double rParallelZ = alignedZ * eZ;

			// vector perpendicular to axis of wire loop, pointing at evaluation pos
			final double rPerpX = r0x - rParallelX;
			final double rPerpY = r0y - rParallelY;
			final double rPerpZ = r0z - rParallelZ;

			// perpendicular distance squared between evalPos and axis of wire loop
			final double alignedRSq = rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ;

			final double rhoP;
			if (alignedRSq > 0.0) {
				// radial unit vector is only defined if evaluation pos is off-axis

				// perpendicular distance between evalPos and axis of wire loop
				final double alignedR = Math.sqrt(alignedRSq);

				// unit vector in radial direction
				final double eRX = rPerpX / alignedR;
				final double eRY = rPerpY / alignedR;
				final double eRZ = rPerpZ / alignedR;

				// normalized rho component of evaluation location in coordinate system of wire loop
				rhoP = alignedR / radius;

				// compute radial component of normalized magnetic field
				// and scale by current and mu_0
				final double bRho = bPrefactor * circularWireLoop_B_rho(rhoP, zP);

				// add contribution from B_rho of wire loop to result
				ret[0][idxEval] = bRho * eRX;
				ret[1][idxEval] = bRho * eRY;
				ret[2][idxEval] = bRho * eRZ;
			} else {
				rhoP = 0.0;
			}

			// compute vertical component of normalized magnetic field
			// and scale by current and mu_0
			final double bZ = bPrefactor * circularWireLoop_B_z(rhoP, zP);

			// add contribution from B_z of wire loop to result
			ret[0][idxEval] += bZ * eX;
			ret[1][idxEval] += bZ * eY;
			ret[2][idxEval] += bZ * eZ;
		}

		return ret;
	}

	/**
	 * Validate a given vector of Cartesian coordinates and compute the number of
	 * elements in it. It is checked that all components are non-null and of equal
	 * length.
	 *
	 * @param a [3: x, y, z][numPos] vector of Cartesian coordinates
	 * @return numPos: number of positions in the given vector
	 */
	private static int validateCartesianVectorInput(double[][] a) {
		Objects.requireNonNull(a);

		if (a.length != 3) {
			throw new RuntimeException("first dimension must have 3 dimension, but has " + a.length);
		}

		Objects.requireNonNull(a[0]);
		Objects.requireNonNull(a[1]);
		Objects.requireNonNull(a[2]);

		int len = a[0].length;
		if (a[1].length != len) {
			throw new RuntimeException("expected " + len + " y coordinates, but only got " + a[1].length);
		}
		if (a[2].length != len) {
			throw new RuntimeException("expected " + len + " z coordinates, but only got " + a[2].length);
		}

		return len;
	}

	// --------------------------------------------------

	/**
	 * Compute the normalized axial component of the magnetic vector potential of a straight wire segment.
	 *
	 * @param rhoP normalized radial coordinate of evaluation location
	 * @param zP normalized axial coordinate of evaluation location
	 * @return normalized axial component of magnetic vector potential
	 */
	public static double straightWireSegment_A_z(double rhoP, double zP) {
		if (rhoP == 0.0) {
			if (zP < 0 || zP > 1) {
				return sws_A_z_ax(zP);
			} else {
				throw new IllegalArgumentException(
						String.format(Locale.ENGLISH, "evaluation locations on the wire segment (rho'=%g z'=%g) are not allowed", rhoP, zP));
			}
		} else if (zP == 0.0 || zP == 1.0) {
			return sws_A_z_rad(rhoP);
		} else if (rhoP >= 1.0 || zP <= -1.0 || zP > 2.0) {
			return sws_A_z_f(rhoP, zP);
		} else {
			return sws_A_z_n(rhoP, zP);
		}
	}

	/**
	 * Compute the normalized tangential component of the magnetic field of a straight wire segment.
	 *
	 * @param rhoP normalized radial coordinate of evaluation location
	 * @param zP normalized axial coordinate of evaluation location
	 * @return normalized tangential component of magnetic field
	 */
	public static double straightWireSegment_B_phi(double rhoP, double zP) {
		if (rhoP == 0.0) {
			if (zP < 0 || zP > 1) {
				return 0.0;
			} else {
				throw new IllegalArgumentException(
						String.format(Locale.ENGLISH, "evaluation locations on the wire segment (rho'=%g z'=%g) are not allowed", rhoP, zP));
			}
		} else if (zP == 0.0 || zP == 1.0) {
			return sws_B_phi_rad(rhoP);
		} else if (rhoP >= zP || rhoP >= 1 - zP || zP < 0.0 || zP > 1.0) {
			return sws_B_phi_f(rhoP, zP);
		} else {
			return sws_B_phi_n(rhoP, zP);
		}
	}

	/**
	 * Geometric part of magnetic vector potential computation for circular wire
	 * loop at rho'=1, z'=0 (normalized coordinates). This routine selects special
	 * case routines to get the most accurate formulation for given evaluation
	 * coordinates.
	 *
	 * @param rhoP normalized radial evaluation position
	 * @param zP   normalized vertical evaluation position
	 * @return A_phi: toroidal component of magnetic vector potential: geometric
	 *         part (no mu0*I/pi factor included)
	 */
	public static double circularWireLoop_A_phi(double rhoP, double zP) {
		if (rhoP == 0.0) {
			return 0.0;
		} else if (rhoP < 0.5 || rhoP > 2.0 || Math.abs(zP) >= 1.0) {
			return cwl_A_phi_f(rhoP, zP);
		} else if (rhoP != 1.0) {
			return cwl_A_phi_n(rhoP, zP);
		} else {
			if (zP != 0) {
				return cwl_A_phi_v(zP);
			} else {
				throw new IllegalArgumentException("evaluation at location of wire loop (rho' = 1, z' = 0) is not defined");
			}
		}
	}

	/**
	 * Geometric part of radial magnetic field computation for circular wire loop at
	 * rho'=1, z'=0 (normalized coordinates). This routine selects special case
	 * routines to get the most accurate formulation for given evaluation
	 * coordinates.
	 *
	 * @param rhoP normalized radial evaluation position
	 * @param zP   normalized vertical evaluation position
	 * @return B_rho: radial component of magnetic field: geometric part (no
	 *         mu0*I/(pi*a) factor included)
	 */
	public static double circularWireLoop_B_rho(double rhoP, double zP) {
		if (rhoP == 0.0 || zP == 0.0) {
			if (rhoP != 1.0) {
				return 0.0;
			} else {
				throw new IllegalArgumentException("evaluation at location of wire loop (rho' = 1, z' = 0) is not defined");
			}
		} else if (rhoP < 0.5 || rhoP > 2.0 || Math.abs(zP) >= 1.0) {
			return cwl_B_rho_f(rhoP, zP);
		} else if (rhoP != 1.0) {
			return cwl_B_rho_n(rhoP, zP);
		} else {
			return cwl_B_rho_v(zP);
		}
	}

	/**
	 * Geometric part of vertical magnetic field computation for circular wire loop
	 * at rho'=1, z'=0 (normalized coordinates). This routine selects special case
	 * routines to get the most accurate formulation for given evaluation
	 * coordinates.
	 *
	 * @param rhoP normalized radial evaluation position
	 * @param zP   normalized vertical evaluation position
	 * @return B_z: vertical component of magnetic field: geometric part (no
	 *         mu0*I/(pi*a) factor included)
	 */
	public static double circularWireLoop_B_z(double rhoP, double zP) {
		if (rhoP < 0.5 || (rhoP <= 2 && Math.abs(zP) > 1)) {
			return cwl_B_z_f1(rhoP, zP);
		} else if (rhoP > 2) {
			return cwl_B_z_f2(rhoP, zP);
		} else if (rhoP != 1.0) {
			return cwl_B_z_n(rhoP, zP);
		} else {
			if (zP != 0) {
				return cwl_B_z_v(zP);
			} else {
				throw new IllegalArgumentException("evaluation at location of wire loop (rho' = 1, z' = 0) is not defined");
			}
		}
	}

	/////// A_z of straight wire segment

	/**
	 * Compute the normalized axial component of magnetic vector potential of straight wire segment,
	 * evaluated along axis of wire segment (rho = 0).
	 *
	 * @param zP normalized axial coordinate of evaluation location; must not be in [0, 1] (on wire segment)
	 * @return normalized axial component of magnetic vector potential
	 */
	static double sws_A_z_ax(double zP) {
		if (zP < -1 || zP >= 2) {
			return sws_A_z_ax_f(zP);
		} else {
			return sws_A_z_ax_n(zP);
		}
	}

	/**
	 * Compute the normalized axial component of magnetic vector potential of straight wire segment,
	 * evaluated along axis of wire segment (rho = 0).
	 * This is a special case for points away from the wire ("far-field") for zP < -1 or zP >= 2.
	 *
	 * @param zP normalized axial coordinate of evaluation location; must not be in [0, 1] (on wire segment)
	 * @return normalized axial component of magnetic vector potential
	 */
	static double sws_A_z_ax_f(double zP) {
		return FastMath.atanh(1 / (Math.abs(zP) + Math.abs(1 - zP)));
	}

	/**
	 * Compute the normalized axial component of magnetic vector potential of straight wire segment,
	 * evaluated along axis of wire segment (rhoP = 0).
	 * This is a special case for points close to the wire ("near-field") for -1 <= zP < 2.
	 *
	 * @param zP normalized axial coordinate of evaluation location; must not be in [0, 1] (on wire segment)
	 * @return normalized axial component of magnetic vector potential
	 */
	static double sws_A_z_ax_n(double zP) {
		return Math.signum(zP) * Math.log(zP / (zP - 1)) / 2;
	}

	/**
	 * Compute the normalized axial component of the magnetic vector potential of a straight wire segment,
	 * evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
	 *
	 * @param rhoP normalized radial coordinate of evaluation location; must not be zero (on wire segment)
	 * @return normalized axial component of magnetic vector potential
	 */
	static double sws_A_z_rad(double rhoP) {
		if (rhoP > 1) {
			return sws_A_z_rad_f(rhoP);
		} else {
			return sws_A_z_rad_n(rhoP);
		}
	}

	/**
	 * Compute the normalized axial component of the magnetic vector potential of a straight wire segment,
	 * evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
	 * This is a special case for points away from the wire ("far-field") for rhoP > 1.
	 *
	 * @param rhoP normalized radial coordinate of evaluation location; must not be zero (on wire segment)
	 * @return normalized axial component of magnetic vector potential
	 */
	static double sws_A_z_rad_f(double rhoP) {
		return FastMath.atanh(1 / (rhoP + Math.hypot(rhoP, 1)));
	}

	/**
	 * Compute the normalized axial component of the magnetic vector potential of a straight wire segment,
	 * evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
	 * This is a special case for points close to the wire ("near-field") for rhoP <= 1.
	 *
	 * @param rhoP normalized radial coordinate of evaluation location; must not be zero (on wire segment)
	 * @return normalized axial component of magnetic vector potential
	 */
	static double sws_A_z_rad_n(double rhoP) {
		double cat = 1 / Math.hypot(rhoP, 1);       // cos(atan(...)  )
		double sat = Math.sin(Math.atan(rhoP) / 2); // sin(atan(...)/2)
		double rc = rhoP * cat;
		double num = rc + 1 + cat;
		double den = rc + 2 * sat * sat;
		return Math.log(num / den) / 2;
	}

	/**
	 * Compute the normalized axial component of the magnetic vector potential of a straight wire segment.
	 * This formulation is useful for points away from the wire ("far-field")
	 * at rhoP >= 1 or zP <= -1 or zP > 2.
	 *
	 * @param rhoP normalized radial coordinate of evaluation location
	 * @param zP normalized axial coordinate of evaluation location
	 * @return normalized axial component of magnetic vector potential
	 */
	static double sws_A_z_f(double rhoP, double zP) {
		double r_i = Math.hypot(rhoP, zP);
		double r_f = Math.hypot(rhoP, 1 - zP);
		return FastMath.atanh(1 / (r_i + r_f));
	}

	/**
	 * Compute the normalized axial component of the magnetic vector potential of a straight wire segment.
	 * This formulation is useful for points close to the wire ("near-field")
	 * at rhoP < 1 and -1 < zP <= 2.
	 *
	 * @param rhoP normalized radial coordinate of evaluation location
	 * @param zP normalized axial coordinate of evaluation location
	 * @return normalized axial component of magnetic vector potential
	 */
	static double sws_A_z_n(double rhoP, double zP) {
		double omz = 1 - zP;

		double r_i = Math.hypot(rhoP, zP);
		double r_f = Math.hypot(rhoP, omz);

		double alpha = Math.atan2(rhoP, zP);
		double sinAlphaHalf = Math.sin(alpha / 2);

		double beta = Math.atan2(rhoP, omz);
		double sinBetaHalf = Math.sin(beta / 2);

		double Ri_zP    = r_i * sinAlphaHalf * sinAlphaHalf; // r_i - z'
		double Rf_p_zM1 = r_f * sinBetaHalf  * sinBetaHalf;  // r_f - (1 - z')

		double n = Ri_zP + Rf_p_zM1;

		return (Math.log(1 + n) - Math.log(n)) / 2;
	}

	/////// B_phi of straight wire segment

	/**
	 * Compute the normalized tangential component of the magnetic field of a straight wire segment,
	 * evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
	 *
	 * @param rhoP normalized radial coordinate of evaluation location
	 * @return normalized tangential component of magnetic field
	 */
	static double sws_B_phi_rad(double rhoP) {
		return 1 / (rhoP * Math.hypot(rhoP, 1));
	}

	/**
	 * Compute the normalized tangential component of the magnetic field of a straight wire segment.
	 * This formulation is useful for points away from the wire ("far-field")
	 * at rhoP >= 1 or zP <= 0 or zP >= 1 or rhoP/(1-zP) >= 1 or rhoP/zP >= 1.
	 *
	 * @param rhoP normalized radial coordinate of evaluation location
	 * @param zP normalized axial coordinate of evaluation location
	 * @return normalized tangential component of magnetic field
	 */
	static double sws_B_phi_f(double rhoP, double zP) {
		double omz = 1 - zP;

		double r_i = Math.hypot(rhoP, zP);
		double r_f = Math.hypot(rhoP, omz);

		double num = rhoP * (1/r_i + 1/r_f);
		double den = rhoP * rhoP - zP * omz + r_i * r_f;

		return num / den;
	}

	/**
	 * Compute the normalized tangential component of the magnetic field of a straight wire segment.
	 * This formulation is useful for points close to the wire ("near-field")
	 * at rhoP < 1 and 0 < zP < 1 and rhoP/(1-zP) < 1 and rhoP/zP < 1.
	 *
	 * @param rhoP normalized radial coordinate of evaluation location
	 * @param zP normalized axial coordinate of evaluation location
	 * @return normalized tangential component of magnetic field
	 */
	static double sws_B_phi_n(double rhoP, double zP) {
		double omz = 1 - zP;

		double r_i = Math.hypot(rhoP, zP);
		double r_f = Math.hypot(rhoP, omz);

		double num = rhoP * (1/r_i + 1/r_f);

		double alpha = Math.atan2(rhoP, zP);
		double sinAlphaHalf = Math.sin(alpha / 2);

		double beta = Math.atan2(rhoP, omz);
		double sinBetaHalf = Math.sin(beta / 2);

		// r_f * sin^2(beta/2) + (1 - z') * sin^2(alpha/2)
		double rfb_omza = r_f * sinBetaHalf * sinBetaHalf + omz * sinAlphaHalf * sinAlphaHalf;

		//     r_i * r_f - z' * (1 - z')
		// =   r_i * r_f - r_i * (1 - z') + r_i * (1 - z') - z' * (1 - z')
		// =   r_i * r_f - r_i * r_f * cos(beta)
		//   + r_i * (1 - z') + (1 - z') * r_i * cos(alpha)
		// =   r_i *    r_f   * (1 - cos(beta))
		//   + r_i * (1 - z') * (1 - cos(alpha))
		// = 2 * r_i * [ r_f * sin^2(beta/2) + (1 - z') * sin^2(alpha/2) ]
		double den = rhoP * rhoP + 2 * r_i * rfb_omza;

		return num / den;
	}

	///// A_phi of circular wire loop

	/**
	 * Compute the normalized tangential component of the magnetic vector potential of a circular wire loop.
	 * This formulation is useful for points away from the wire ("far-field")
	 * at rhoP < 1/2 or rhoP > 2 or |zP| >= 1.
	 *
	 * @param rhoP normalized radial coordinate of evaluation location
	 * @param zP normalized axial coordinate of evaluation location
	 * @return normalized tangential component of magnetic vector potential
	 */
	static double cwl_A_phi_f(double rhoP, double zP) {
		double sqrt_kCSqNum = Math.hypot(zP, 1 - rhoP);
		double sqrt_kCSqDen = Math.hypot(zP, 1 + rhoP);

		double kC = sqrt_kCSqNum / sqrt_kCSqDen;
		double kSq = 4 * rhoP / (sqrt_kCSqDen * sqrt_kCSqDen);

		double kCp1 = 1 + kC;
		double arg1 = 2 * Math.sqrt(kC) / kCp1;
		double arg2 = 2 / (kCp1 * kCp1 * kCp1);
		double C = CompleteEllipticIntegral.cel(arg1, 1, 0, arg2);

		return kSq/sqrt_kCSqDen * C;
	}

	/**
	 * Compute the normalized tangential component of the magnetic vector potential of a circular wire loop.
	 * This formulation is useful for points close to the wire ("near-field")
	 * at 1/2 <= rhoP <= 2 and |zP| < 1.
	 *
	 * @param rhoP normalized radial coordinate of evaluation location
	 * @param zP normalized axial coordinate of evaluation location
	 * @return normalized tangential component of magnetic vector potential
	 */
	static double cwl_A_phi_n(double rhoP, double zP) {
		double rhoP_m_1 = rhoP - 1;

		double n = zP / rhoP_m_1;
		double m = 1 + 2 / rhoP_m_1;

		double num = n * n + 1;
		double den = n * n + m * m;

		double kCSq = num / den;

		double prefac = 1 / (Math.abs(rhoP - 1) * Math.sqrt(den));
		double celPart = CompleteEllipticIntegral.cel(Math.sqrt(kCSq), 1, -1, 1);
		return prefac * celPart;
	}

	/**
	 * Compute the normalized tangential component of the magnetic vector potential of a circular wire loop.
	 * This formulation is useful for points along rhoP=1 with |zP| < 1.
	 *
	 * @param zP normalized axial coordinate of evaluation location
	 * @return normalized tangential component of magnetic vector potential
	 */
	static double cwl_A_phi_v(double zP) {
		double absZp = Math.abs(zP);

		// 1/k_c
		double kCInv = Math.sqrt(4 + zP * zP) / absZp;

		return CompleteEllipticIntegral.cel(kCInv, 1, 1, -1) / absZp;
	}

	//////// B_rho of circular wire loop

	/**
	 * Compute the normalized radial component of the magnetic field of a circular wire loop.
	 * This formulation is useful for points away from the wire ("far-field")
	 * at rhoP < 1/2 or rhoP > 2 or |zP| >= 1.
	 *
	 * @param rhoP normalized radial coordinate of evaluation location
	 * @param zP normalized axial coordinate of evaluation location
	 * @return normalized radial component of magnetic field
	 */
	static double cwl_B_rho_f(double rhoP, double zP) {
		double sqrt_kCSqNum = Math.hypot(zP, 1 - rhoP);
		double sqrt_kCSqDen = Math.hypot(zP, 1 + rhoP);

		double kCSqNum = sqrt_kCSqNum * sqrt_kCSqNum;
		double kCSqDen = sqrt_kCSqDen * sqrt_kCSqDen;

		double kCSq = kCSqNum / kCSqDen;
		double kC = Math.sqrt(kCSq);

		double D = CompleteEllipticIntegral.cel(kC, 1, 0, 1);

		double kCp1 = 1 + kC;
		double arg1 = 2 * Math.sqrt(kC) / kCp1;
		double arg2 = 2 / (kCp1 * kCp1 * kCp1);
		double C = CompleteEllipticIntegral.cel(arg1, 1, 0, arg2);

		double prefac = 4 * rhoP / (kCSqDen * sqrt_kCSqDen * kCSqNum);

		return prefac * zP * (D - C);
	}

	/**
	 * Compute the normalized radial component of the magnetic field of a circular wire loop.
	 * This formulation is useful for points close to the wire ("near-field")
	 * at 1/2 <= rhoP <= 2 and |zP| < 1.
	 *
	 * @param rhoP normalized radial coordinate of evaluation location
	 * @param zP normalized axial coordinate of evaluation location
	 * @return normalized radial component of magnetic field
	 */
	static double cwl_B_rho_n(double rhoP, double zP) {
		double rhoP_m_1 = rhoP - 1;
		double rd2 = rhoP_m_1 * rhoP_m_1;

		double n = zP / rhoP_m_1;
		double m = 1 + 2 / rhoP_m_1;

		double sqrt_kCSqNum = Math.hypot(n, 1);
		double sqrt_kCSqDen = Math.hypot(n, m);

		double kCSqNum = sqrt_kCSqNum * sqrt_kCSqNum;
		double kCSqDen = sqrt_kCSqDen * sqrt_kCSqDen;

		double kC = sqrt_kCSqNum / sqrt_kCSqDen;

		double D = CompleteEllipticIntegral.cel(kC, 1, 0, 1);

		double kCp1 = 1 + kC;
		double arg1 = 2 * Math.sqrt(kC) / kCp1;
		double arg2 = 2 / (kCp1 * kCp1 * kCp1);
		double C = arg2 * CompleteEllipticIntegral.cel(arg1, 1, 0, 1);

		// z' / |rho' - 1|^5
		double zP_rd5 = zP / (Math.abs(rhoP_m_1) * rd2 * rd2);

		double prefac = 4 * rhoP / (kCSqDen * sqrt_kCSqDen * kCSqNum);

		return prefac * zP_rd5 * (D - C);
	}

	/**
	 * Compute the normalized radial component of the magnetic field of a circular wire loop.
	 * This formulation is useful for points along rhoP=1 with |zP| < 1.
	 *
	 * @param zP normalized axial coordinate of evaluation location
	 * @return normalized radial component of magnetic field
	 */
	static double cwl_B_rho_v(double zP) {
		double zPSq = zP * zP;

		double kCSq = 1 / (1 + 4 / zPSq);
		double kC = Math.sqrt(kCSq);

		double K = CompleteEllipticIntegral.cel(kC, 1, 1, 1);
		double E = CompleteEllipticIntegral.cel(kC, 1, 1, kCSq);

		return Math.signum(zP) * kC / 2 * ((2 / zPSq + 1) * E - K);
	}

	////// B_z of circular wire loop

	/**
	 * Compute the normalized vertical component of the magnetic field of a circular wire loop.
	 * This formulation is useful for certain points away from the wire ("far-field")
	 * at rhoP < 1/2 or (rhoP <= 2 and |zP| >= 1).
	 *
	 * @param rhoP normalized radial coordinate of evaluation location
	 * @param zP normalized axial coordinate of evaluation location
	 * @return normalized vertical component of magnetic field
	 */
	static double cwl_B_z_f1(double rhoP, double zP) {
		double sqrt_kCSqNum = Math.hypot(zP, 1 - rhoP);
		double sqrt_kCSqDen = Math.hypot(zP, 1 + rhoP);

		double kC = sqrt_kCSqNum / sqrt_kCSqDen;

		double K = CompleteEllipticIntegral.cel(kC, 1, 1, 1);
		double E = CompleteEllipticIntegral.cel(kC, 1, 1, kC * kC);
		double D = CompleteEllipticIntegral.cel(kC, 1, 0, 1);

		double prefac = 1 / (sqrt_kCSqDen * sqrt_kCSqNum * sqrt_kCSqNum);
		double comb = (E - 2 * K + 2 * D);

		return prefac * (E + rhoP * comb);
	}

	/**
	 * Compute the normalized vertical component of the magnetic field of a circular wire loop.
	 * This formulation is useful for certain other points away from the wire ("far-field")
	 * at rhoP > 2.
	 *
	 * @param rhoP normalized radial coordinate of evaluation location
	 * @param zP normalized axial coordinate of evaluation location
	 * @return normalized vertical component of magnetic field
	 */
	static double cwl_B_z_f2(double rhoP, double zP) {
		double sqrt_kCSqNum = Math.hypot(zP, 1 - rhoP);
		double sqrt_kCSqDen = Math.hypot(zP, 1 + rhoP);

		double kC = sqrt_kCSqNum / sqrt_kCSqDen;
		double kCSq = kC * kC;

		double zPSqP1 = zP * zP + 1;
		double rhoPSq = rhoP * rhoP;
		double t1 = zPSqP1 / rhoPSq + 1;
		double t2 = 2 / rhoP;

		// a is sqrt_kCSqDen normalized to rho'^2
		// b is sqrt_kCSqNum normalized to rho'^2
		// a == (z'^2 + (1 + rho')^2) / rho'^2 = (z'^2 + 1)/rho'^2 + 1  +  2/rho'
		// b == (z'^2 + (1 - rho')^2) / rho'^2 = (z'^2 + 1)/rho'^2 + 1  -  2/rho'
		double a = t1 + t2;
		double b = t1 - t2;

		// 1/prefac = sqrt( z'^2 + (1 + rho')^2)           * (z'^2 + (1 - rho')^2)
		//          = sqrt((z'^2 + (1 + rho')^2) / rho'^2) * (z'^2 + (1 - rho')^2) / rho'^2 * rho'^3
		//          = sqrt(a)                              * b                              * rho'^3
		double prefac = 1 / (Math.sqrt(a) * b * rhoPSq * rhoP);

		double cdScale = 1 + (2 + zPSqP1 / rhoP) / rhoP;

		double E = CompleteEllipticIntegral.cel(kC, 1, 1, kCSq);
		double D = CompleteEllipticIntegral.cel(kC, 1, 0, 1);

		double kCP1 = 1 + kC;
		double arg1 = 2 * Math.sqrt(kC) / kCP1;
		double arg2 = 2 / (kCP1 * kCP1 * kCP1);
		double C = arg2 * CompleteEllipticIntegral.cel(arg1, 1, 0, 1);

		// use C - D for (2 * D - E)/kSq
		return prefac * (E + 4 * (C - D) / cdScale);
	}

	/**
	 * Compute the normalized vertical component of the magnetic field of a circular wire loop.
	 * This formulation is useful for points close to the wire ("near-field")
	 * at 1/2 <= rhoP <= 2, but not rhoP=1, and |zP| <= 1.
	 *
	 * @param rhoP normalized radial coordinate of evaluation location
	 * @param zP normalized axial coordinate of evaluation location
	 * @return normalized vertical component of magnetic field
	 */
	static double cwl_B_z_n(double rhoP, double zP) {
		double rp1 = rhoP - 1;

		double n = zP / rp1;
		double m = 1 + 2 / rp1;

		double sqrt_kCSqNum = Math.hypot(n, 1);
		double sqrt_kCSqDen = Math.hypot(n, m);

		double kCSqDen = sqrt_kCSqDen * sqrt_kCSqDen;

		double kC = sqrt_kCSqNum / sqrt_kCSqDen;

		double prefac = 1 / (Math.abs(rp1) * rp1 * rp1 * kCSqDen * sqrt_kCSqDen);

		return prefac * CompleteEllipticIntegral.cel(kC, kC * kC, 1 + rhoP, 1 - rhoP);
	}

	/**
	 * Compute the normalized vertical component of the magnetic field of a circular wire loop.
	 * This formulation is useful for points along rhoP=1 with |zP| <= 1.
	 *
	 * @param zP normalized axial coordinate of evaluation location
	 * @return normalized vertical component of magnetic field
	 */
	static double cwl_B_z_v(double zP) {
		double kCSq = zP * zP / (4 + zP * zP);
		double kC = Math.sqrt(kCSq);

		double f = zP * zP + 4;
		double prefac = 1 / (f * Math.sqrt(f));

		return prefac * CompleteEllipticIntegral.cel(kC, kCSq, 2, 0);
	}
}
