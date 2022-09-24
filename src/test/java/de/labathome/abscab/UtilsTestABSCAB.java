package de.labathome.abscab;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;

public class UtilsTestABSCAB {

	/**
	 * Write a given array to a text file; one line per row.
	 *
	 * @param arr      [numRows] array to write to a file
	 * @param filename file into which to write the given array
	 */
	public static void dumpToFile(double[] arr, String filename) {
		dumpToFile(new double[][] {arr}, filename);
	}

	/**
	 * Write a given array to a text file; one line per row.
	 *
	 * @param arr      [numRows] array to write to a file
	 * @param filename file into which to write the given array
	 */
	public static void dumpToFile(int[] arr, String filename) {
		dumpToFile(new int[][] {arr}, filename);
	}

	/**
	 * Write a given array to a text file; one line per row.
	 *
	 * @param arr      [numColumns][numRows] array to write to a file
	 * @param filename file into which to write the given array
	 */
	public static void dumpToFile(double[][] arr, String filename) {
		File outFile = new File(filename);
		try (PrintWriter pw = new PrintWriter(outFile)) {
			for (int idxRow = 0; idxRow < arr[0].length; ++idxRow) {
				String line = "";
				for (int idxCol=0; idxCol<arr.length; ++idxCol) {
					line += String.format(Locale.ENGLISH, "%+.20e ", arr[idxCol][idxRow]);
				}
				pw.println(line.trim());
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * Write a given array to a text file; one line per row.
	 *
	 * @param arr      [numColumns][numRows] array to write to a file
	 * @param filename file into which to write the given array
	 */
	public static void dumpToFile(int[][] arr, String filename) {
		File outFile = new File(filename);
		try (PrintWriter pw = new PrintWriter(outFile)) {
			for (int idxRow = 0; idxRow < arr[0].length; ++idxRow) {
				String line = "";
				for (int idxCol=0; idxCol<arr.length; ++idxCol) {
					line += String.format("%d ", arr[idxCol][idxRow]);
				}
				pw.println(line.trim());
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * Write the given set of test points, such that the implied values can be
	 * represented exactly in arbitrary-precision software.
	 * The format of an IEEE754 double precision number F is:
	 * <pre>
	 * F = (-1)^s * 2^{E - 1023} * (1 + M/2^{52})
	 * </pre>
	 * where:
	 * <pre>
	 * s = 0, 1            ( 1 bit )
	 * E = 0, 1, ..., 2047 (11 bits)
	 * M = 0, 1, ...       (52 bits)
	 * </pre>
	 * The test point coordinates (rp, zp) are stored in six columns in the output file:
	 * (s, E, M) for rp and then (s, E, M) for zp.
	 * This allows to re-construct the number that is actually implied
	 * by a given double precision variable within arbitrary precision software.
	 *
	 * @param testPointsRp [numCases] set of rp test point coordinates
	 * @param testPointsZp [numCases] set of zp test point coordinates
	 * @param filename     file into which to write the given test points
	 */
	public static void dumpTestPoints(double[] testPointsRp, double[] testPointsZp, String filename) {
		int numCases = testPointsRp.length;

		File outFile = new File(filename);
		try (PrintWriter pw = new PrintWriter(outFile)) {
			pw.println("# rp: sign bit, exponent E, mantiassa M; zp: sign bit, exponent E, mantiassa M");

			for (int i = 0; i < numCases; ++i) {

				long[] rpParts = doubleParts(testPointsRp[i]);
				long[] zpParts = doubleParts(testPointsZp[i]);

				long signRp     = rpParts[0];
				long exponentRp = rpParts[1];
				long mantissaRp = rpParts[2];

				long signZp     = zpParts[0];
				long exponentZp = zpParts[1];
				long mantissaZp = zpParts[2];

				pw.printf(Locale.ENGLISH, "%1d %4d %16d %1d %4d %16d\n",
						signRp, exponentRp, mantissaRp,
						signZp, exponentZp, mantissaZp);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * Decompose a IEEE 754 double precision variable F as:
	 * <pre>
	 * F = (-1)^s * 2^{E - 1023} * (1 + M/2^{52})
	 * </pre>
	 * where:
	 * <pre>
	 * s = 0, 1            ( 1 bit )
	 * E = 0, 1, ..., 2047 (11 bits)
	 * M = 0, 1, ...       (52 bits)
	 * </pre>
	 * @param f variable to decompose
	 * @return {s, E, M}
	 */
	public static long[] doubleParts(double f) {
		//                          64   60   56   52   48   44   40   36   32   28   24   20   16   12   8    4
		final long signMask     = 0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000L;
		final long exponentMask = 0b0111_1111_1111_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000L;
		final long mantissaMask = 0b0000_0000_0000_1111_1111_1111_1111_1111_1111_1111_1111_1111_1111_1111_1111_1111L;

		long bits = Double.doubleToRawLongBits(f);

		long sign     = (bits &     signMask) >>> 63;
		long exponent = (bits & exponentMask) >>> 52;
		long mantissa = (bits & mantissaMask);

		return new long[] {sign, exponent, mantissa};
	}

	/**
	 * Load columns of data from a text file.
	 *
	 * @param filename file to load data from
	 * @return [numColumns][numRows] data loaded from file
	 */
	public static double[][] loadColumnsFromFile(String filename) {
		try {
			List<String> lines = Files.readAllLines(Paths.get(filename));
			return parseColumnsFromLines(lines);
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Load columns of data from a resource.
	 *
	 * @param clazz class to load resource from
	 * @param resourceName name of resource to load data from
	 * @return [numColumns][numRows] data loaded from file
	 */
	public static double[][] loadColumnsFromResource(Class<?> clazz, String resourceName) {
		InputStream is = clazz.getResourceAsStream(resourceName);
		if (is == null) {
			throw new RuntimeException("resource '"+resourceName+"' not found");
		}

		InputStreamReader isr = new InputStreamReader(is);
		try(BufferedReader br = new BufferedReader(isr)) {

			List<String> lines = new LinkedList<>();
			String line = br.readLine();
			while (line != null) {
				lines.add(line);
				line = br.readLine();
			}

			return parseColumnsFromLines(lines);
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Parse columns of data from a given list of lines.
	 * Lines starting with '#' are ignored
	 *
	 * @param lines [numLines] lines of text to parse data from; whitespace-separated items per line
	 * @return [numColumns][numRows] data parsed from given lines
	 */
	private static double[][] parseColumnsFromLines(List<String> lines) {

		// pass 1: count available data lines
		int numRows = 0;
		int numColumns = 0;
		for (String line : lines) {
			// ignore comment lines
			if (!line.startsWith("#")) {

				if (numRows == 0) {
					// first non-comment line defines how many columns there are
					String[] parts = line.trim().split("\\s+");
					numColumns = parts.length;
				}

				// count how many non-comment lines there are
				// --> number of data rows
				numRows++;
			}
		}

		double[][] data = new double[numColumns][numRows];

		// pass 2: actually parse data, line by line
		int idxRow = 0;
		for (String line : lines) {
			// ignore comment lines
			if (!line.startsWith("#")) {
				String[] parts = line.trim().split("\\s+");

				for (int idxCol = 0; idxCol < numColumns; ++idxCol) {
					data[idxCol][idxRow] = Double.valueOf(parts[idxCol]);
				}

				idxRow++;
			}
		}

		return data;
	}

	public static double errorMetric(double ref, double act) {

		double bad = 0.0;
		double good = -16.0;

		if (Math.abs(ref) > 0.0) {
			if (act != ref) {
				return Math.log10(Math.min(Math.pow(10, bad), Math.abs((act-ref)/ref)));
			} else {
				return good;
			}
		} else {
			return ((Math.abs(act) > 0.0) ? bad : good);
		}
	}
}
