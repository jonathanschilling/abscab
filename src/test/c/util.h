#ifndef UTIL_H
#define UTIL_H

#include <ctype.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/**
 * Check if two values are approximately equal within a prescribed tolerance.
 * For values much smaller than 1, this is similar to a comparison of the
 * absolute values. For values much greater than 1, this is similar to a
 * comparison of the relative values.
 *
 * This method is described in Gill, Murray & Wright, "Practical Optimization" (1984).
 *
 * @param expected  expected result
 * @param actual    actual result
 * @param tolerance relative or absolute tolerance on the mismatch between the
 *                  expected and the actual values
 * @return 0 if the values match within the prescribed tolerance; 1 otherwise
 */
int assertRelAbsEquals(double expected, double actual, double tolerance) {
	double relAbsError = fabs(actual - expected) / (1.0 + fabs(expected));
	if (relAbsError > tolerance) {
		printf("expected %g, actual %g (rel/abs error %g, tolerance %g)\n",
				expected, actual, relAbsError, tolerance);
		return 1;
	}
	return 0;
}

/**
 * Error metric for computing the "number of matching digits" between two numbers.
 */
double errorMetric(double ref, double act) {
	double bad = 0.0;
	double good = -16.0;

	double tenToBad = pow(10, bad);

	if (fabs(ref) > 0.0) {
		if (act != ref) {
			double relErr = fabs((act-ref)/ref);
			return log10(tenToBad < relErr ? tenToBad : relErr);
		} else {
			return good;
		}
	} else {
		return ((fabs(act) > 0.0) ? bad : good);
	}
}

// https://stackoverflow.com/a/122721
// Note: This function returns a pointer to a substring of the original string.
// If the given string was allocated dynamically, the caller must not overwrite
// that pointer with the returned value, since the original pointer must be
// deallocated using the same allocator with which it was allocated.  The return
// value must NOT be deallocated using free() etc.
char* trim_whitespace(char *str) {

	// Trim leading space
	while (isspace((unsigned char) *str)) {
		str++;
	}

	// All spaces?
	if (*str == 0) {
		return str;
	}

	// Trim trailing space
	char *end = str + strlen(str) - 1;
	while (end > str && isspace((unsigned char) *end)) {
		end--;
	}

	// Write new null terminator character
	end[1] = '\0';

	return str;
}

/**
 * Load columns of data from a text file.
 * The number of rows found will be written into numRows.
 * The number of columns found will be written into numColumns.
 * The return value is an array of arrays (pointer to an array of pointers),
 * where one entry corresponds to one column of the data in the file.
 * Rows in the file starting with "#" are not counted and not parsed.
 */
double** loadColumnsFromFile(char *filename, int *numRows, int *numColumns) {

	// try to open given file
	FILE *fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("failed to open file '%s'\n", filename);
		return NULL;
	}

	int bufSize = 256;
	char *line = (char *) malloc(bufSize);
	if (line == NULL) {
		printf("could not allocate read buffer\n");
		return NULL;
	}

	int rows = 0;
	int cols = 0;
	char *delims = " \t";

	// first pass: read file and count lines and columns
	while (fgets(line, bufSize, fp) != NULL) {
		// skip empty lines and comment lines
		if (line[0] != '#') {
			rows++;

			// now count how many columns there are in the current line
			char* trimmed_line = trim_whitespace(line);

			// https://en.cppreference.com/w/c/string/byte/strtok
			int cols_in_this_line = 0;
			char *token = strtok(trimmed_line, delims);
			while (token) {
				cols_in_this_line++;
				token = strtok(NULL, delims);
			}

			// cols <-- max(cols, cols_in_this_line)
			cols = (cols_in_this_line > cols ? cols_in_this_line : cols);
		}
	}

//	printf("found %d rows of up to %d column(s)\n", rows, cols);

	// allocate target array
	double** data = (double**) malloc(cols * sizeof(double*));
	for (int i = 0; i < cols; ++i) {
		data[i] = (double*) malloc(rows * sizeof(double));
	}

	// seek back to start of file
	int status = fseek(fp, 0, SEEK_SET);
	if (status) {
		printf("failed to seek to start of file '%s': status = %d\n", filename, status);
		return NULL;
	}

	// second pass: read data from file and parse into allocated array
	int row = 0;
	while (fgets(line, bufSize, fp) != NULL) {
		// skip empty lines and comment lines
		if (line[0] != '#') {

			// now count how many columns there are in the current line
			char* trimmed_line = trim_whitespace(line);

			int col = 0;
			char *token = strtok(trimmed_line, delims);
			while (token) {

//				printf("parse '%s' into data[%d][%d]", token, col, row);

				// parse data
				data[col][row] = atof(token);

//				printf(" => %g\n", data[col][row]);

				token = strtok(NULL, delims);
				col++;
			}
			row++;
		}
	}

	// close file and check that this was successful
	status = fclose(fp);
	if (status) {
		printf("failed to close file '%s': status = %d\n", filename, status);
	}

	free(line);

	// return row and column counts
	*(numRows) = rows;
	*(numColumns) = cols;

	return data;
}

void dumpToFile(int numCols, int numRows, double *data, char *filename) {

	FILE *fp = fopen(filename, "w+");
	if (fp == NULL) {
		printf("failed to open file '%s'\n", filename);
		return;
	}

	for (int row = 0; row < numRows; ++row) {
		for (int col = 0; col < (numCols-1); ++col) {
			fprintf(fp, "%.20e ", data[row * numCols + col]);
		}
		fprintf(fp, "%.20e\n", data[row * numCols + (numCols-1)]);
	}

	int status = fclose(fp);
	if (status) {
		printf("failed to close file '%s': status = %d\n", filename, status);
	}
}

#endif // UTIL_H
