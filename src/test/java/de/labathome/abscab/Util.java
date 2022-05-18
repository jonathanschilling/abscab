package de.labathome.abscab;

import java.io.File;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Locale;

public class Util {

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
				pw.println(line.strip());
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

		//                          64   60   56   52   48   44   40   36   32   28   24   20   16   12   8    4
		final long signMask     = 0b1000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000L;
		final long exponentMask = 0b0111_1111_1111_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000_0000L;
		final long mantissaMask = 0b0000_0000_0000_1111_1111_1111_1111_1111_1111_1111_1111_1111_1111_1111_1111_1111L;

		File outFile = new File(filename);
		try (PrintWriter pw = new PrintWriter(outFile)) {
			pw.println("# rp: sign bit, exponent E, mantiassa M; zp: sign bit, exponent E, mantiassa M");

			for (int i = 0; i < numCases; ++i) {

				long rpBits = Double.doubleToRawLongBits(testPointsRp[i]);
				long zpBits = Double.doubleToRawLongBits(testPointsZp[i]);

				long signRp     = (rpBits &     signMask) >>> 63;
				long exponentRp = (rpBits & exponentMask) >>> 52;
				long mantissaRp = (rpBits & mantissaMask);

				long signZp     = (zpBits &     signMask) >>> 63;
				long exponentZp = (zpBits & exponentMask) >>> 52;
				long mantissaZp = (zpBits & mantissaMask);

				pw.printf(Locale.ENGLISH, "%1d %4d %16d %1d %4d %16d\n",
						signRp, exponentRp, mantissaRp,
						signZp, exponentZp, mantissaZp);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
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

			int numRows = 0;
			int numColumns = 0;
			for (String line : lines) {
				// ignore comment lines
				if (!line.startsWith("#")) {

					if (numRows == 0) {
						// first non-comment line defines how many columns there are
						String[] parts = line.strip().split("\\s+");
						numColumns = parts.length;
					}

					// count how many non-comment lines there are
					// --> number of data rows
					numRows++;
				}
			}

			double[][] data = new double[numColumns][numRows];

			int idxRow = 0;
			for (String line : lines) {
				// ignore comment lines
				if (!line.startsWith("#")) {
					String[] parts = line.strip().split("\\s+");

					for (int idxCol = 0; idxCol < numColumns; ++idxCol) {
						data[idxCol][idxRow] = Double.valueOf(parts[idxCol]);
					}

					idxRow++;
				}
			}

			return data;
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}
}
