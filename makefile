run_tests: clean tests
	./tests

clean:
	-rm tests

tests: tests.c
