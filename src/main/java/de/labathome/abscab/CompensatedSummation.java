package de.labathome.abscab;

/**
 * 2nd-order Kahan-Babuska summation algorithm as listed in A. Klein, "A
 * Generalized Kahan-Babuska-Summation-Algorithm", Computing 76, 279-293 (2006)
 * doi: 10.1007/s00607-005-0139-x
 */
public class CompensatedSummation {

	/** summation variable */
	private double s;

	/** first-order correction */
	private double cs;

	/** second-order correction */
	private double ccs;

	/** Instantiate a new summation object. */
	public CompensatedSummation() {
		reset();
	}

	/** Reset the accumulated value to zero. */
	public void reset() {
		s = 0.0;
		cs = 0.0;
		ccs = 0.0;
	}

	/**
	 * Add a number of contributions to the sum.
	 *
	 * @param contributions contributions to add to the sum
	 */
	public void add(double[] contributions) {
		for (double contribution : contributions) {
			add(contribution);
		}
	}

	/**
	 * Add a single contribution to the sum.
	 *
	 * @param contribution contribution to add to the sum
	 */
	public void add(double contribution) {
		final double t = s + contribution;
		final double c;
		if (Math.abs(s) >= Math.abs(contribution)) {
			c = (s - t) + contribution;
		} else {
			c = (contribution - t) + s;
		}
		s = t;

		final double t2 = cs + c;
		final double cc;
		if (Math.abs(cs) >= Math.abs(c)) {
			cc = (cs - t2) + c;
		} else {
			cc = (c - t2) + cs;
		}
		cs = t2;
		ccs += cc;
	}

	/**
	 * Compute the second-order corrected sum of all contributions.
	 *
	 * @return compensated sum of all contributions since creation or last reset
	 */
	public double getSum() {
		return s + cs + ccs;
	}
}
