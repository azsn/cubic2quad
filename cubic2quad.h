// cubic2quad generates a spline of quadratic beziers to approximate a single
// cubic bezier.
//
// Parameters:
// in: The input cubic bezier in the form
//     p1x, p1y, c1x, c1y, c2x, c2y, p2x, p2y
//
// precision: How close the output spline should be to the original cubic.
//     Smaller values for precision will result in a more accurate spline but
//     will require more quadratic beziers to form it.
//
// out: The output quadratic beziers, a repetition of 6 doubles
//     p1x, p1y, cx, cy, p2x, p2y
//     Note that (p2x,p2y) of one quadratic will always equal the (p1x,p1y)
//     of the next quadratic because they are placed end-to-end.
//
// return value: The number of output quadratics written to `out`.
//     `out` is filled with [return value]*6 doubles.
int cubic2quad(const double in[8], const double precision, double out[144])
