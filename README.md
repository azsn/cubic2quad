# cubic2quad

A C library to convert cubic bezier curves to quadratic bezier curves by
approximation.

**This library is a near-direct translation of
[github.com/fontello/cubic2quad](https://github.com/fontello/cubic2quad) to C.**
This version tries to closely match the logic of the original, with only some
minor changes to avoid requiring heap memory allocation.

## Usage

This library has only one function, `cubic2quad()`. See
[`cubic2quad.h`](cubic2quad.h) for usage details. The simplest way to use this
code is to directly copy `cubic2quad.c`/`.h` into your project.

## Tests

The [`tests.c`](tests.c) file ports all of the tests provided from the original JS version.
It also has a test that tries to compare the original version to this one,
which requires some manual configuration -- see `test_compare_to_original()`
in `tests.c` for details.

To run tests, run `make`. No output means all tests passed with no problems.

## License

[MIT](LICENSE)
