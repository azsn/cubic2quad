#include <stdio.h>
#include <stdlib.h>
#include "cubic2quad.c"

#define assertTrue(a) do { if (!(a)) { fprintf(stderr, "assertion failed (value: %llu). line %d\n", (unsigned long long)a, __LINE__); } } while(0)

#define assertEqual(a, b) do { if ((a) != (b)) { fprintf(stderr, "assertion failed: %llu != %llu. line %d\n", (unsigned long long)(a), (unsigned long long)(b), __LINE__); } } while(0)

#define assertClose(a, b) do { if (fabs((a) - (b)) > 1e-15) { fprintf(stderr, "assertion failed: %f not close to %f. line %d\n", (a), (b), __LINE__); } } while(0)

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
	// equation with zero coefficients
	{
		double roots[3];
		int n = cubic_solve(0, 0, 0, 0, roots);
		assertEqual(n, 0);
	}

	// linear equation
	{
		double roots[3];
		int n = cubic_solve(0, 0, 1, -1, roots);
		assertEqual(n, 1);
		assertEqual(roots[0], 1);
	}

	// quadratic equation with no real roots
	{
		double roots[3];
		int n = cubic_solve(0, 1, 2, 2, roots);
		assertEqual(n, 0);
	}

	// quadratic equation with one real root
	{
		double roots[3];
		int n = cubic_solve(0, 1, 2, 1, roots);
		assertEqual(n, 1);
		assertClose(roots[0], -1.0);
	}

	// quadratic equation with two real roots
	{
		double roots[3];
		int n = cubic_solve(0, 1, 1, 0, roots);
		assertEqual(n, 2);
		sort_doubles(roots, 2);
		assertClose(roots[0], -1.0);
		assertClose(roots[1], 0.0);
	}

	// cubic equation with one real root
	{
		double roots[3];
		int n = cubic_solve(1, 0, 0, 1, roots);
		assertEqual(n, 1);
		assertClose(roots[0], -1.0);
	}

	// cubic equation with two real roots
	{
		double roots[3];
		int n = cubic_solve(1, 1, 0, 0, roots);
		assertEqual(n, 2);
		sort_doubles(roots, 2);
		assertClose(roots[0], -1.0);
		assertClose(roots[1], 0.0);
	}

	// cubic equation with three real roots
	{
		double roots[3];
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

int main() {
	test_cubic_equation_solver();
	test__is_approximation_close();
	return 0;
}
