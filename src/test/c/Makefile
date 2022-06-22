

ABSCAB_DIR=../../main/c
ABSCAB_H=${ABSCAB_DIR}/abscab.h

.PHONY: all test

test_cel: ${ABSCAB_H} test_cel.c util.h
	gcc -I${ABSCAB_DIR} -lm -o test_cel test_cel.c 

all: test_cel

test:
	@ ./test_cel

clean:
	rm -f test_cel
