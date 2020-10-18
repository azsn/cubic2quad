#include <stdio.h>
#include <stdlib.h>
#include "cubic2quad.c"

#define assertEqual(a, b) do { if (a != b) { fprintf(stderr, "assertion failed: %llu != %llu. line %d\n", (unsigned long long)a, (unsigned long long)b, __LINE__); } } while(0)

#define assertClose(a, b) do { if (fabs(a - b) > 1e-15) { fprintf(stderr, "assertion failed: %f not close to %f. line %d\n", a, b, __LINE__); } } while(0)

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

int main() {
	test_cubic_equation_solver();
	return 0;
}
