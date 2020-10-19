#include <stdio.h>
#include <stdlib.h>
#include "cubic2quad.c"

#define assertTrue(a) do { \
	if (!(a)) { \
		fprintf(stderr, "assertion failed (value: %d). line %d\n", (bool)(a), __LINE__); \
	} \
} while(0)

#define assertEqual(a, b) do { \
	if ((a) != (b)) { \
		fprintf(stderr, "assertion failed: %llu != %llu. line %d\n", (unsigned long long)(a), (unsigned long long)(b), __LINE__); \
	} \
} while(0)

#define assertCloseRes(a, b, res) do { \
	if (fabs((a) - (b)) > (res)) { \
		fprintf(stderr, "assertion failed: %f not close to %f. line %d\n", (a), (b), __LINE__); \
	} \
} while(0)
#define assertClose(a, b) assertCloseRes(a, b, 1e-15)

#define assertArraysCloseRes(a, b, n, res) do { \
	for (int i = 0; i < n; i++) { \
		if (fabs((a)[i] - (b)[i]) > (res)) { \
			fprintf(stderr, "assertion failed: %f not close to %f (index %d). line %d\n", (a)[i], (b)[i], (i), __LINE__); \
		} \
	} \
} while(0)
#define assertArraysClose(a, b, n) assertArraysCloseRes(a, b, n, 1e-15)

static int cmp_doubles(const void *pa, const void *pb)
{
	double a = *((double*)pa), b = *((double*)pb);
	if (a > b) return 1;
	if (a < b) return -1;
	return 0;
}

static void sort_doubles(double *arr, int n)
{
	qsort(arr, n, sizeof(double), cmp_doubles);
}

static void test_cubic_equation_solver()
{
	double roots[3];

	// equation with zero coefficients
	{
		int n = cubic_solve(0, 0, 0, 0, roots);
		assertEqual(n, 0);
	}

	// linear equation
	{
		int n = cubic_solve(0, 0, 1, -1, roots);
		assertEqual(n, 1);
		assertEqual(roots[0], 1);
	}

	// quadratic equation with no real roots
	{
		int n = cubic_solve(0, 1, 2, 2, roots);
		assertEqual(n, 0);
	}

	// quadratic equation with one real root
	{
		int n = cubic_solve(0, 1, 2, 1, roots);
		assertEqual(n, 1);
		assertClose(roots[0], -1.0);
	}

	// quadratic equation with two real roots
	{
		int n = cubic_solve(0, 1, 1, 0, roots);
		assertEqual(n, 2);
		sort_doubles(roots, 2);
		assertClose(roots[0], -1.0);
		assertClose(roots[1], 0.0);
	}

	// cubic equation with one real root
	{
		int n = cubic_solve(1, 0, 0, 1, roots);
		assertEqual(n, 1);
		assertClose(roots[0], -1.0);
	}

	// cubic equation with two real roots
	{
		int n = cubic_solve(1, 1, 0, 0, roots);
		assertEqual(n, 2);
		sort_doubles(roots, 2);
		assertClose(roots[0], -1.0);
		assertClose(roots[1], 0.0);
	}

	// cubic equation with three real roots
	{
		int n = cubic_solve(1, 0, -1, 0, roots);
		assertEqual(n, 3);
		sort_doubles(roots, 3);
		assertClose(roots[0], -1.0);
		assertClose(roots[1], 0.0);
		assertClose(roots[2], 1.0);
	}
}

static bool is_approximation_close(
	// cubic
	double p1x, double p1y, double c1x, double c1y, double c2x, double c2y, double p2x, double p2y,
	// quadratics
	const QBezier * const quadCurves, const int quadCurvesLen,
	double errorBound)
{
	Point pc[4];
	calc_power_coefficients(
		p_new(p1x, p1y),
		p_new(c1x, c1y),
		p_new(c2x, c2y),
		p_new(p2x, p2y),
		pc);
	return _is_approximation_close(pc[0], pc[1], pc[2], pc[3], quadCurves, quadCurvesLen, errorBound);
}

static void test__is_approximation_close()
{
	// straight line is the same as a straight line (error ~ 1e-8)
	{
		QBezier expect = { { 0, 0 }, { 15, 15 }, { 30, 30 } };
		assertTrue(is_approximation_close(0, 0, 10, 10, 20, 20, 30, 30,
			&expect, 1, 1e-8));
	}

	// two parallel lines (distance between is 0.01) are not same (error ~ 1e-8)
	{
		QBezier expect = { { 0, 0.01 }, { 15, 0.01 }, { 30, 0.01 } };
		assertTrue(!is_approximation_close(0, 0, 10, 0, 20, 0, 30, 0,
			&expect, 1, 1e-8));
	}

	// two parallel lines (distance between is 0.01) are same (error 0.1)
	{
		QBezier expect = { { 0, 0.01 }, { 15, 0.01 }, { 30, 0.01 } };
		assertTrue(is_approximation_close(0, 0, 10, 0, 20, 0, 30, 0,
			&expect, 1, 0.1));
	}

	// actually quadratic curve is the same as converted (error ~ 1e-8)
	{
		QBezier expect = { { 0, 0 }, { 15, 15 }, { 30, 0 } };
		assertTrue(is_approximation_close(0, 0, 10, 10, 20, 10, 30, 0,
			&expect, 1, 1e-8));
	}

	// two quadratic curves (distance between CP is 0.01) are not same (error ~ 1e-8)
	{
		QBezier expect = { { 0, 0 }, { 15, 15.01 }, { 30, 0 } };
		assertTrue(!is_approximation_close(0, 0, 10, 10, 20, 10, 30, 0,
			&expect, 1, 1e-8));
	}

	// two quadratic curves (distance between CP is 0.01) are same (error ~ 0.1)
	{
		QBezier expect = { { 0, 0 }, { 15, 15.01 }, { 30, 0 } };
		assertTrue(is_approximation_close(0, 0, 10, 10, 20, 10, 30, 0,
			&expect, 1, 0.1));
	}
}

static void test_cubic2quad()
{
	double out[MAX_DOUBLES_OUT];

	// straight line to the same straight line (error ~ 1e-8)
	{
		double in[] = { 0, 0, 10, 10, 20, 20, 30, 30 };
		double expect[] = { 0, 0, 15, 15, 30, 30 };
		int n = cubic2quad(in, 1e-8, out);
		assertEqual(n, 1);
		assertArraysClose(out, expect, 6);
	}

	// quadratic curve to the same quadratic curve (error ~ 1e-8)
	{
		double in[] = { 0, 0, 10, 10, 20, 10, 30, 0 };
		double expect[] = { 0, 0, 15, 15, 30, 0 };
		int n = cubic2quad(in, 10, out);
		assertEqual(n, 1);
		assertArraysClose(out, expect, 6);
	}

	// cubic curve that is close to quadratic to a quadratic one (error ~ 0.1)
	{
		double in[] = { 0, 0, 10, 9, 20, 11, 30, 0 };
		int n = cubic2quad(in, 0.1, out);
		assertEqual(n, 1);
	}

	// should be able to handle inflections (error ~ 0.5)
	{
		double in[] = { 858, -113, 739, -68, 624, -31, 533, 0 };
		int n = cubic2quad(in, 0.5, out);
		// All points should be in the bounding box of the source curve
		for (int i = 0; i < (n * 6); i += 2) {
			assertTrue(out[i] <= in[0]);
			assertTrue(out[i] >= in[6]);
			assertTrue(out[i + 1] >= in[1]);
			assertTrue(out[i + 1] <= in[7]);
		}
	}

	// should be able to handle inflections (error ~ 0.005)
	{
		double in[] = { 858, -113, 739, -68, 624, -31, 533, 0 };
		int n = cubic2quad(in, 0.005, out);
		// All points should be in the bounding box of the source curve
		for (int i = 0; i < (n * 6); i += 2) {
			assertTrue(out[i] <= in[0]);
			assertTrue(out[i] >= in[6]);
			assertTrue(out[i + 1] >= in[1]);
			assertTrue(out[i + 1] <= in[7]);
		}
	}

	// should split curve at inflection point
	{
		double in[] = { 0, 100, 70, 0, 30, 0, 100, 100 };
		cubic2quad(in, 1000, out);
		assertCloseRes(out[4], 34.33, 0.01);
		assertCloseRes(out[5], 45.45, 0.01);
		assertCloseRes(out[10], 65.67, 0.01);
		assertCloseRes(out[11], 45.45, 0.01);
	}

	// cubic curve have to be converted to two or more quads for large errorBound (error ~ 100)
	{
		double in[] = { 0, 0, -5, 10, 35, 10, 30, 0 };
		int n = cubic2quad(in, 100, out);
		assertTrue(n > 1);
	}
}

static void test_compare_to_original()
{
	/*
	This test compares the C cubic2quad with the original JS version by
	having the original version generate lots of random input cubics and
	"expected output" quadratic sequences. Use the following JS (Node)
	code to generate the test data, then use the C code it generates as
	input (see 'Replace' comment below) for this test.

const count = 100; // number of beziers to generate, can adjust as needed
const res = 0.1; // resolution / errorBound, can adjust too
const cubic2quad = require('cubic2quad');
function s(n) { return Number(n.toFixed(5)); }
function decompress(p) {
	let result = [];
	for (let i = 0; i < p.length - 2; i += 4) {
		result.push(s(p[i]), s(p[i+1]), s(p[i+2]), s(p[i+3]), s(p[i+4]), s(p[i+5]));
	}
	return result;
}
let inp = [];
for (let i = 0; i < (count * 8); i++) {
	inp.push(s((Math.random() * 100) - 50));
}
let out = [];
for (let i = 0; i < (count * 8); i += 8) {
	const arr = decompress(cubic2quad(...inp.slice(i, i+8), res));
	out = out.concat([arr.length, ...arr]);
}
console.log(`const double res = ${res}, count = ${count};`);
console.log(`const double in[] = { ${inp.join()} };`);
console.log(`const double expect[] = { ${out.join()} };`);

	*/

	/// Replace these with the output generated by the script above:
	const double res = 0, count = 0;
	const double in[] = {};
	const double expect[] = {};
	///

	double out[MAX_DOUBLES_OUT];
	for (int i = 0, j = 0; i < count; i++, j += 1+((int)expect[j]) ) {
		int nexpect = ((int)expect[j])/6;
		int n = cubic2quad(&in[i*8], res, out);
		if (n != nexpect) {
			fprintf(stderr, "bezier %d (%d != %d)\n", i, n, nexpect);
			// Minor differences in floating point measurements
			// can result in different numbers of output quads.
			// So just check that both the expected and output
			// beziers are close to the original cubic.
			assertTrue(is_approximation_close(
				in[i*8], in[i*8+1], in[i*8+2], in[i*8+3],
				in[i*8+4], in[i*8+5], in[i*8+6], in[i*8+7],
				(QBezier *)&expect[j+1], nexpect, res));
			assertTrue(is_approximation_close(
				in[i*8], in[i*8+1], in[i*8+2], in[i*8+3],
				in[i*8+4], in[i*8+5], in[i*8+6], in[i*8+7],
				(QBezier *)out, n, res));
			continue;
		}
		assertArraysCloseRes(out, &expect[j+1], (int)expect[j], 0.00001);
	}
}

int main() {
	test_cubic_equation_solver();
	test__is_approximation_close();
	test_cubic2quad();
	test_compare_to_original();
	return 0;
}
