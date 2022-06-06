package de.labathome.abscab;

import java.util.Objects;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;

import org.apache.commons.math3.util.FastMath;

import de.labathome.CompleteEllipticIntegral;

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
			Function<Integer, double[]> vertexSupplier,
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
			Function<Integer, double[]> vertexSupplier,
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
			Function<Integer, double[]> vertexSupplier,
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
			Function<Integer, double[]> vertexSupplier,
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
			Function<Integer, double[]> vertexSupplier,
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
			Function<Integer, double[]> vertexSupplier,
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

			final int nThreads;

			if (numVertices-1 > numEvalPos) {
				// parallelize over nSource-1

				// Note that each thread needs its own copy of the vectorPotential array,
				// so this approach might need quite some memory in case the number of
				// threads and the number of evaluation points is large.

				final int nSourcePerThread;
				if (numVertices-1 < numProcessors) {
					nThreads = numVertices-1;
					nSourcePerThread = 1;
				} else {
					nThreads = numProcessors;

					// It is better that many threads do more
					// than one thread needs to do more.
					nSourcePerThread = (int) Math.ceil( (numVertices-1.0) / nThreads) ;
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
			Function<Integer, double[]> vertexSupplier,
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

			final int nThreads;

			if (numVertices-1 > numEvalPos) {
				// parallelize over nSource-1

				// Note that each thread needs its own copy of the vectorPotential array,
				// so this approach might need quite some memory in case the number of
				// threads and the number of evaluation points is large.

				final int nSourcePerThread;
				if (numVertices-1 < numProcessors) {
					nThreads = numVertices-1;
					nSourcePerThread = 1;
				} else {
					nThreads = numProcessors;

					// It is better that many threads do more
					// than one thread needs to do more.
					nSourcePerThread = (int) Math.ceil( ((double)(numVertices - 1)) / nThreads) ;
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

			final int nThreads;

			if (numVertices-1 > numEvalPos) {
				// parallelize over nSource-1

				// Note that each thread needs its own copy of the vectorPotential array,
				// so this approach might need quite some memory in case the number of
				// threads and the number of evaluation points is large.

				final int nSourcePerThread;
				if (numVertices-1 < numProcessors) {
					nThreads = numVertices-1;
					nSourcePerThread = 1;
				} else {
					nThreads = numProcessors;

					// It is better that many threads do more
					// than one thread needs to do more.
					nSourcePerThread = (int) Math.ceil( (numVertices-1.0) / nThreads) ;
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
			Function<Integer, double[]> vertexSupplier,
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

			final int nThreads;

			if (numVertices-1 > numEvalPos) {
				// parallelize over nSource-1

				// Note that each thread needs its own copy of the vectorPotential array,
				// so this approach might need quite some memory in case the number of
				// threads and the number of evaluation points is large.

				final int nSourcePerThread;
				if (numVertices-1 < numProcessors) {
					nThreads = numVertices-1;
					nSourcePerThread = 1;
				} else {
					nThreads = numProcessors;

					// It is better that many threads do more
					// than one thread needs to do more.
					nSourcePerThread = (int) Math.ceil( ((double)(numVertices - 1)) / nThreads) ;
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

				// r0 projected onto axis of wire segment
				final double rParallelX = alignedZ * eX;
				final double rParallelY = alignedZ * eY;
				final double rParallelZ = alignedZ * eZ;

				// vector perpendicular to axis of wire segment, pointing at evaluation pos
				final double rPerpX = r0x - rParallelX;
				final double rPerpY = r0y - rParallelY;
				final double rPerpZ = r0z - rParallelZ;

				// perpendicular distance between evalPos and axis of wire segment
				final double alignedR = Math.sqrt(rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ);

				// normalized rho component of evaluation location in coordinate system of wire segment
				final double rhoP = alignedR / l;

				// compute parallel component of magnetic vector potential, including current and mu_0
				final double aParallel = aPrefactor * straightWireSegment_A_z(rhoP, zP);

				// add contribution from wire segment to result
				if (useCompensatedSummation) {
					aXSum[idxEval].add(aParallel * eX);
					aYSum[idxEval].add(aParallel * eY);
					aZSum[idxEval].add(aParallel * eZ);
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
			Function<Integer, double[]> vertexSupplier,
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

				// r0 projected onto axis of wire segment
				final double rParallelX = alignedZ * eX;
				final double rParallelY = alignedZ * eY;
				final double rParallelZ = alignedZ * eZ;

				// vector perpendicular to axis of wire segment, pointing at evaluation pos
				final double rPerpX = r0x - rParallelX;
				final double rPerpY = r0y - rParallelY;
				final double rPerpZ = r0z - rParallelZ;

				// perpendicular distance between evalPos and axis of wire segment
				final double alignedR = Math.sqrt(rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ);

				// normalized rho component of evaluation location in coordinate system of wire segment
				final double rhoP = alignedR / l;

				// compute parallel component of magnetic vector potential, including current and mu_0
				final double aParallel = aPrefactor * straightWireSegment_A_z(rhoP, zP);

				// add contribution from wire segment to result
				if (useCompensatedSummation) {
					aXSum[idxEval].add(aParallel * eX);
					aYSum[idxEval].add(aParallel * eY);
					aZSum[idxEval].add(aParallel * eZ);
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

				// r0 projected onto axis of wire segment
				final double rParallelX = alignedZ * eX;
				final double rParallelY = alignedZ * eY;
				final double rParallelZ = alignedZ * eZ;

				// vector perpendicular to axis of wire segment, pointing at evaluation pos
				final double rPerpX = r0x - rParallelX;
				final double rPerpY = r0y - rParallelY;
				final double rPerpZ = r0z - rParallelZ;

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
						bXSum[idxEval].add(bPhi * ePhiX);
						bYSum[idxEval].add(bPhi * ePhiY);
						bZSum[idxEval].add(bPhi * ePhiZ);
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
			Function<Integer, double[]> vertexSupplier,
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

				// r0 projected onto axis of wire segment
				final double rParallelX = alignedZ * eX;
				final double rParallelY = alignedZ * eY;
				final double rParallelZ = alignedZ * eZ;

				// vector perpendicular to axis of wire segment, pointing at evaluation pos
				final double rPerpX = r0x - rParallelX;
				final double rPerpY = r0y - rParallelY;
				final double rPerpZ = r0z - rParallelZ;

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
						bXSum[idxEval].add(bPhi * ePhiX);
						bYSum[idxEval].add(bPhi * ePhiY);
						bZSum[idxEval].add(bPhi * ePhiZ);
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

			// perpendicular distance between evalPos and axis of wire loop
			final double alignedR = Math.sqrt(rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ);

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
			ret[0][idxEval] += aPhi * ePhiX;
			ret[1][idxEval] += aPhi * ePhiY;
			ret[2][idxEval] += aPhi * ePhiZ;
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

			// perpendicular distance between evalPos and axis of wire loop
			final double alignedR = Math.sqrt(rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ);

			final double rhoP;
			if (alignedR > 0.0) {
				// radial unit vector is only defined if evaluation pos is off-axis

				// unit vector in radial direction
				final double eRX = rPerpX / alignedR;
				final double eRY = rPerpY / alignedR;
				final double eRZ = rPerpZ / alignedR;

				// normalized rho component of evaluation location in coordinate system of wire loop
				rhoP = alignedR / radius;

				// compute radial component of normalized magnetic field
				// and scale by current and mu_0
				final double bRho = bPrefactor * circularWireLoop_B_rho(rhoP, zP);

				// add contribution from wire loop to result
				ret[0][idxEval] += bRho * eRX;
				ret[1][idxEval] += bRho * eRY;
				ret[2][idxEval] += bRho * eRZ;
			} else {
				rhoP = 0.0;
			}

			// compute vertical component of normalized magnetic field
			// and scale by current and mu_0
			final double bZ = bPrefactor * circularWireLoop_B_z(rhoP, zP);

			// add contribution from wire loop to result
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


	/**
	 * Full-field magnetic vector potential of straight wire segment; only A_z
	 * component is present.
	 *
	 * @param rhoP
	 * @param zP
	 * @return
	 */
	public static double straightWireSegment_A_z(double rhoP, double zP) {
		if (rhoP == 0.0) {
			return A_z_along_rhoP_0(rhoP, zP);
		} else if (zP == 0.0 || zP == 1.0) {
			return A_z_along_zP_0_or_1(rhoP, zP);
		} else if (-1 < zP && zP <= 2.0 && rhoP < 1.0) {
			// near-field
			if (zP >= 1.0 || rhoP / (1 - zP) >= 1.0) {
				return A_z_6a(rhoP, zP);
			} else if (zP >= 0.0 && rhoP / zP <= 1.0) {
				return A_z_6b(rhoP, zP);
			} else {
				return A_z_6c(rhoP, zP);
			}
		} else {
			return A_z_1(rhoP, zP);
		}
	}

	/**
	 * Full-field magnetic field of straight wire segment; only B_phi component is
	 * present.
	 *
	 * @param rhoP
	 * @param zP
	 * @return
	 */
	public static double straightWireSegment_B_phi(double rhoP, double zP) {
		if (rhoP == 0.0) {
			// line along rho'=0: exactly 0
			return 0.0;
		} else if (zP == 0.0 || zP == 1.0) {
			// line along z'=0 or z'=1
			return B_phi_3(rhoP, zP);
		} else if (zP >= 1.0 || zP <= 0.0 || rhoP >= 1.0 || rhoP / (1 - zP) >= 1.0 || rhoP / zP >= 1.0) {
			// far-field
			return B_phi_4(rhoP, zP);
		} else {
			// near field between z'=0 and z'=1
			return B_phi_5(rhoP, zP);
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
		} else if (zP >= 1.0 || rhoP < 0.5 || rhoP > 2.0) {
			return A_phi_1(rhoP, zP);
		} else if (rhoP != 1.0) {
			return A_phi_6(rhoP, zP);
		} else {
			// rhoP == 1, zP < 1
			return A_phi_5(rhoP, zP);
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
			return 0.0;
		} else if (zP >= 1.0 || rhoP < 0.5 || rhoP > 2.0) {
			return B_rho_3(rhoP, zP);
		} else if (rhoP != 1.0) {
			return B_rho_1(rhoP, zP);
		} else {
			return B_rho_4(rhoP, zP);
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
		if (rhoP < 0.5 || (rhoP <= 2 && zP > 1)) {
			return B_z_1(rhoP, zP);
		} else if (rhoP > 2) {
			return B_z_2(rhoP, zP);
		} else if (rhoP == 1.0) {
			return B_z_4(rhoP, zP);
		} else if (zP != 0.0) {
			return B_z_5(rhoP, zP);
		} else {
			return B_z_6(rhoP, zP);
		}
	}

	/////// A_z of straight wire segment

	/**
	 * Straight-forward implementation of A_z for straight wire segment.
	 * Useful for far-field.
	 * @param rhoP
	 * @param zP
	 * @return
	 */
	static double A_z_1(double rhoP, double zP) {
		double Ri = Math.hypot(rhoP, zP);
		double Rf = Math.hypot(rhoP, 1.0 - zP);
		return FastMath.atanh(1.0 / (Ri + Rf));
	}

	/**
	 * combined solution for rhoP=0
	 *
	 * @param rhoP
	 * @param zP
	 * @return
	 */
	static double A_z_along_rhoP_0(double rhoP, double zP) {
		if (zP < -1 || zP > 2) {
			return A_z_2(rhoP, zP);
		} else {
			return A_z_2b(rhoP, zP);
		}
	}

	/**
	 * special case for rho'=0
	 *
	 * @param rhoP
	 * @param zP
	 * @return
	 */
	static double A_z_2(double rhoP, double zP) {
		return FastMath.atanh(1.0 / (Math.abs(zP) + Math.abs(1 - zP)));
	}

	/**
	 * special case for rho'=0; near-field (excluding wire)
	 *
	 * @param rhoP
	 * @param zP
	 * @return
	 */
	static double A_z_2b(double rhoP, double zP) {
		return Math.signum(zP) * Math.log(Math.abs(zP / (1 - zP))) / 2;
	}

	/**
	 * combined solution for zP=0 or zP = 1
	 *
	 * @param rhoP
	 * @param zP
	 * @return
	 */
	static double A_z_along_zP_0_or_1(double rhoP, double zP) {
		if (rhoP > 1.0) {
			return A_z_3(rhoP, 0.0);
		} else {
			// rhoP <= 1
			return A_z_3b(rhoP, 0.0);
		}
	}

	/**
	 * special case for z'=0
	 *
	 * @param rhoP
	 * @param zP
	 * @return
	 */
	static double A_z_3(double rhoP, double zP) {
		return FastMath.atanh(1.0 / (rhoP + Math.sqrt(rhoP * rhoP + 1)));
	}

	/**
	 * special case for z'=0
	 *
	 * @param rhoP
	 * @param zP
	 * @return
	 */
	static double A_z_3b(double rhoP, double zP) {
		// a little bit more robust --> around rho'=1 +/- one test point we have 15 digits
		double cat = 1 / Math.sqrt(rhoP * rhoP + 1); // cos(atan(...))
		double sat = Math.sin(Math.atan(rhoP) / 2); // sin(atan(...)/2)
		double num = rhoP * cat + 1 + cat;
		double den = rhoP * cat + 2 * sat * sat;
		return Math.log(num / den) / 2;
	}

	// (1): rho' < 1e-15, |z'|>=1
	static double A_z_6a(double rhoP, double zP) {

		double alpha = Math.atan2(rhoP, zP);
		double cosAlpha = Math.cos(alpha);
		double sinAlphaHalf = Math.sin(alpha / 2);

		// nutritious zero: R_i - 1 == (R_i - z') + (z' - 1)
		double Ri_zP = zP * 2 * sinAlphaHalf * sinAlphaHalf / cosAlpha; // R_i - z'

		double zpM1 = zP - 1.0;
		double Rf = Math.sqrt(rhoP * rhoP + zpM1 * zpM1);

		double n = Ri_zP + Rf + zpM1;

		return (Math.log(2 + n) - Math.log(n)) / 2;
	}

	static double A_z_6b(double rhoP, double zP) {

		double alpha = Math.atan2(rhoP, zP);
		double cosAlpha = Math.cos(alpha);
		double sinAlphaHalf = Math.sin(alpha / 2);

		// nutritious zero: R_i - 1 == (R_i - z') + (z' - 1)
		double Ri_zP = 2 * zP * sinAlphaHalf * sinAlphaHalf / cosAlpha; // R_i - z'

		double omz = 1.0 - zP;
		double beta = Math.atan2(rhoP, omz);
		double cosBeta = Math.cos(beta);
		double sinBetaHalf = Math.sin(beta / 2);

		double Rf_p_zM1 = 2 * omz * sinBetaHalf * sinBetaHalf / cosBeta; // R_f - 1 + z'

		double n = Ri_zP + Rf_p_zM1;

		return (Math.log(2 + n) - Math.log(n)) / 2;
	}

	static double A_z_6c(double rhoP, double zP) {
		double omz = 1.0 - zP;

		double beta = Math.atan2(rhoP, omz);
		double sinBetaHalf = Math.sin(beta / 2);

		double R_i = Math.hypot(rhoP, zP);
		double R_f = Math.hypot(rhoP, omz);

		double Rf_m_1 = 2.0 * R_f * sinBetaHalf * sinBetaHalf - zP;

		double n = R_i + Rf_m_1;

		return (Math.log(2 + n) - Math.log(n)) / 2;
	}

	/////// B_phi of straight wire segment

	/**
	 * special case for zP=0 or zP=1
	 *
	 * @param rhoP
	 * @param zP
	 * @return
	 */
	static double B_phi_3(double rhoP, double zP) {
		return 1 / (rhoP * Math.hypot(rhoP, 1));
	}

	// straight-forward implementation; good for far-field
	static double B_phi_4(double rhoP, double zP) {

		double Ri = Math.hypot(rhoP, zP);
		double Rf = Math.hypot(rhoP, 1 - zP);

		double t1 = rhoP * (1/Ri + 1/Rf);

		double den = rhoP * rhoP + zP * (zP - 1) + Ri * Rf;

		return t1 / den;
	}

	// near-field: zP approx 1
	static double B_phi_5(double rhoP, double zP) {

		double Ri = Math.hypot(rhoP, zP);
		double Rf = Math.hypot(rhoP, 1 - zP);

		double t1 = rhoP * (1/Ri + 1/Rf);

		double alpha = Math.atan2(rhoP, zP);
		double cosAlpha = Math.cos(alpha);
		double sinAlphaHalf = Math.sin(alpha / 2.0);

		double beta = Math.atan2(rhoP, 1 - zP);
		double cosBeta = Math.cos(beta);
		double sinBetaHalf = Math.sin(beta / 2.0);

		// (a*b - 1)
		double abm1 = 2.0 / cosAlpha * (sinBetaHalf * sinBetaHalf / cosBeta + sinAlphaHalf * sinAlphaHalf);

		// R_i*R_f - zP*(1-zP) == zP*(1-zP) * (a*b - 1)
		double den = rhoP * rhoP + zP * (1.0 - zP) * abm1;

		return t1 / den;
	}

	///// A_phi of circular wire loop

	static double A_phi_1(double rhoP, double zP) {
		// complementary modulus k_c
		double t = zP * zP + 1 + rhoP * rhoP;
		double kCSqNum = t - 2.0 * rhoP;
		double kCSqDen = t + 2.0 * rhoP;

		double sqrt_kCSqNum = Math.sqrt(kCSqNum);
		double sqrt_kCSqDen = Math.sqrt(kCSqDen);
		double kC = sqrt_kCSqNum / sqrt_kCSqDen;

		double kSq = 4.0 * rhoP / kCSqDen;

		double celPrefactor = 1.0 / sqrt_kCSqDen;

		// Walstrom
		double arg1 = 2 * Math.sqrt(kC) / (1 + kC);
		double a2d = 1 + kC;
		double arg2 = 2 / (a2d * a2d * a2d);
		double C = CompleteEllipticIntegral.cel(arg1, 1, 0, arg2);

		return celPrefactor * kSq * C;
	}

	/**
	 * special case for rho'=1 and z' close to 0
	 *
	 * @param rhoP
	 * @param zP
	 * @return
	 */
	static double A_phi_5(double rhoP, double zP) {
		double kc = Math.sqrt(4.0 + zP * zP) / zP;
		return 1.0 / zP * CompleteEllipticIntegral.cel(kc, 1, 1, -1);
	}

	/**
	 *
	 * @param rhoP
	 * @param zP
	 * @return
	 */
	static double A_phi_6(double rhoP, double zP) {
		double n = zP / (rhoP - 1);
		double m = 1 + 2 / (rhoP - 1);
		double den = n * n + m * m;
		double num = n * n + 1;
		double kcsq = num / den;

		double prefac = 1 / (Math.abs(rhoP - 1) * Math.sqrt(den));
		double celPart = CompleteEllipticIntegral.cel(Math.sqrt(kcsq), 1, -1, 1);
		return prefac * celPart;
	}

	//////// B_rho of circular wire loop

	static double B_rho_1(double rhoP, double zP) {
		double n = zP / (rhoP - 1);
		double m = 1 + 2 / (rhoP - 1);
		double den = n * n + m * m;
		double num = n * n + 1;
		double kCSq = num / den;

		double kC = Math.sqrt(kCSq);

		double D = CompleteEllipticIntegral.cel(kC, 1, 0, 1);

		double arg1 = 2 * Math.sqrt(kC) / (1 + kC);
		double a2d = 1 + kC;
		double arg2 = 2 / (a2d * a2d * a2d);
		double C = CompleteEllipticIntegral.cel(arg1, 1, 0, arg2);

		double rd = rhoP - 1;
		double rd2 = rd * rd;

		double fac1 = 4 * rhoP / (rd2 * rd2) * Math.abs(n);

		return fac1 * (D - C) / (den * Math.sqrt(den) * num);
	}

	static double B_rho_3(double rhoP, double zP) {

		double sqrt_kCSqNum = Math.hypot(zP, 1 - rhoP);
		double sqrt_kCSqDen = Math.hypot(zP, 1 + rhoP);
		double k = 2 * Math.sqrt(rhoP) / sqrt_kCSqDen;
		double kSq = k * k;
		double kCSq = 1.0 - kSq;
		double kC = Math.sqrt(kCSq);

		double D = CompleteEllipticIntegral.cel(kC, 1, 0, 1);

		double arg1 = 2 * Math.sqrt(kC) / (1 + kC);
		double a2d = 1 + kC;
		double arg2 = 2 / (a2d * a2d * a2d);
		double C = CompleteEllipticIntegral.cel(arg1, 1, 0, arg2);

		double prefac = 4 * rhoP * zP / (sqrt_kCSqDen * sqrt_kCSqDen * sqrt_kCSqDen * sqrt_kCSqNum * sqrt_kCSqNum);

		return prefac * (D - C);
	}

	/** special case for rhoP=1, zP --> 0 */
	static double B_rho_4(double rhoP, double zP) {

		double zPSq = zP * zP;
		double pfd = 1 + 4 / zPSq;
		double kCSq = 1 / pfd;
		double kC = Math.sqrt(kCSq);

		double E = CompleteEllipticIntegral.cel(kC, 1, 1, kCSq);
		double K = CompleteEllipticIntegral.cel(kC, 1, 1, 1);

		return 1 / (2 * Math.sqrt(pfd)) * (E / pfd * (1 + (6 + 8 / zPSq) / zPSq) - K);
	}

	////// B_z of circular wire loop

	static double B_z_1(double rhoP, double zP) {

		double sqrt_kCSqNum = Math.hypot(zP, 1 - rhoP);
		double kCSqNum = sqrt_kCSqNum * sqrt_kCSqNum;

		double sqrt_kCSqDen = Math.hypot(zP, 1 + rhoP);

		double k = 2 * Math.sqrt(rhoP) / sqrt_kCSqDen;
		double kSq = k * k;
		double kCSq = 1.0 - kSq;
		double kC = Math.sqrt(kCSq);

		double E = CompleteEllipticIntegral.cel(kC, 1, 1, kCSq);
		double K = CompleteEllipticIntegral.cel(kC, 1, 1, 1);
		double D = CompleteEllipticIntegral.cel(kC, 1, 0, 1);

		double comb = (E - 2 * K + 2 * D);

		double prefac = 1.0 / (sqrt_kCSqDen * kCSqNum);

		return prefac * (E + rhoP * comb);
	}

	static double B_z_2(double rhoP, double zP) {
		// large rho'

		double sqrt_kCSqDen = Math.hypot(zP, 1 + rhoP);
		double k = 2 * Math.sqrt(rhoP) / sqrt_kCSqDen;
		double kSq = k * k;
		double kCSq = 1.0 - kSq;
		double kC = Math.sqrt(kCSq);

		double E = CompleteEllipticIntegral.cel(kC, 1, 1, kCSq);
		double D = CompleteEllipticIntegral.cel(kC, 1, 0, 1);

		double zp2_1 = zP * zP + 1;
		double r6 = zp2_1 * zp2_1 * zp2_1;
		double r5 = r6 / rhoP - 2 * zp2_1 * zp2_1;
		double r4 = r5 / rhoP + 3 * zp2_1 * zp2_1 - 4 * zp2_1;
		double r3 = r4 / rhoP - 4 * zP * zP + 4;
		double r2 = r3 / rhoP + 3 * zP * zP - 1;
		double r1 = r2 / rhoP - 2;
		double r0 = r1 / rhoP + 1;

		// use C-D for (2D-E)/kSq: this finally works without expansion !!!
		double arg1 = 2 * Math.sqrt(kC) / (1 + kC);
		double a2d = 1 + kC;
		double arg2 = 2 / (a2d * a2d * a2d);
		double C = CompleteEllipticIntegral.cel(arg1, 1, 0, arg2);

		double cdScale = 1 + (2 + (zP * zP + 1) / rhoP) / rhoP;

		return 1.0 / (Math.sqrt(r0) * rhoP * rhoP * rhoP) * (E + 4 * (C - D) / cdScale);
	}

	static double B_z_4(double rhoP, double zP) {
		// special case for rhoP=1, zp->0

		double kCSq = zP * zP / (4 + zP * zP);
		double kC = Math.sqrt(kCSq);

		double f = zP * zP + 4;
		double prefac = 1 / (f * Math.sqrt(f));

		return prefac * CompleteEllipticIntegral.cel(kC, kCSq, 2, 0);
	}

	static double B_z_5(double rhoP, double zP) {
		// special case for near-field: rhoP->1, zP->0; but not rhoP=1 or zP=0

		double rp1 = rhoP - 1;

		double n = zP / rp1;
		double m = 1 + 2 / rp1;
		double den = n * n + m * m;
		double num = n * n + 1;
		double kCSq = num / den;

		double prefac = Math.abs(n) / (rp1 * rp1 * den * Math.sqrt(den));

		double ca1 = (1 + rhoP) / zP;
		double ca2 = (1 - rhoP) / zP;
		double cp = CompleteEllipticIntegral.cel(Math.sqrt(kCSq), kCSq, ca1, ca2);

		return prefac * cp;
	}

	static double B_z_6(double rhoP, double zP) {

		double rp = rhoP + 1;

		double kC = (1 - rhoP) / rp;
		double kCSq = kC * kC;

		double ca1 = 1 + rhoP;
		double ca2 = 1 - rhoP;
		double cp = CompleteEllipticIntegral.cel(Math.sqrt(kCSq), kCSq, ca1, ca2);

		return 1 / Math.abs(rp * rp * rp) * cp;
	}
}
