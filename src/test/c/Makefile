
CC?=gcc
CFLAGS=--std=c99

ABSCAB_DIR=../../main/c

.PHONY: all test clean demo_programs test_programs

all: test_programs demo_programs

test_programs: test_cel test_abscab

demo_programs: demo_abscab

test: test_programs
	@ ./test_cel
	@ ./test_abscab
	
demo: demo_programs
	@ ./demo_abscab

clean:
	rm -f test_cel test_abscab demo_abscab

test_cel: $(ABSCAB_DIR)/cel.h util.h test_cel.c
	$(CC) $(CFLAGS) test_cel.c -I$(ABSCAB_DIR) -o test_cel -lm

test_abscab: $(ABSCAB_DIR)/abscab.h util.h test_abscab.c
	$(CC) $(CFLAGS) test_abscab.c -I$(ABSCAB_DIR) -o test_abscab -lm

demo_abscab: $(ABSCAB_DIR)/abscab.h util.h demo_abscab.c
	$(CC) $(CFLAGS) demo_abscab.c -I$(ABSCAB_DIR) -o demo_abscab -lm -O3 -fopenmp
