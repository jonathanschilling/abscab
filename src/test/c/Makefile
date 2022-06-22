
ABSCAB_DIR=../../main/c

.PHONY: all test clean

test_cel: ${ABSCAB_DIR}/cel.h test_cel.c util.h
	gcc -I${ABSCAB_DIR} -lm -o test_cel test_cel.c 

all: test_cel

test:
	@ ./test_cel

clean:
	rm -f test_cel
