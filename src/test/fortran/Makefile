
F08?=gfortran
CFLAGS=-std=f2008 -fopenmp -O3

ABSCAB_DIR=../../main/fortran

.PHONY: all test clean test_programs demo_programs

all: test_programs demo_programs

test_programs: test_cel test_abscab

demo_programs: demo_abscab

test: test_programs
	@ ./test_cel
	@ ./test_abscab

demo: demo_programs
	@ ./demo_abscab

clean:
	rm -f *.mod *.o
	rm -f test_cel test_abscab demo_abscab

mod_kinds.o: $(ABSCAB_DIR)/mod_kinds.f08
	$(F08) $(CFLAGS) $(ABSCAB_DIR)/mod_kinds.f08 -c

mod_cel.o: mod_kinds.o $(ABSCAB_DIR)/mod_cel.f08
	$(F08) $(CFLAGS) $(ABSCAB_DIR)/mod_cel.f08 -c

mod_compsum.o: mod_kinds.o $(ABSCAB_DIR)/mod_compsum.f08
	$(F08) $(CFLAGS) $(ABSCAB_DIR)/mod_compsum.f08 -c

mod_testutil.o: mod_kinds.o mod_testutil.f08
	$(F08) $(CFLAGS) mod_testutil.f08 -c

test_cel: mod_cel.o test_cel.f08
	$(F08) $(CFLAGS) test_cel.f08 -c
	$(F08) $(CFLAGS) mod_kinds.o mod_cel.o test_cel.o -o test_cel

abscab.o: mod_cel.o mod_compsum.o $(ABSCAB_DIR)/abscab.f08
	$(F08) $(CFLAGS) $(ABSCAB_DIR)/abscab.f08 -c

test_abscab: abscab.o mod_testutil.o test_abscab.f08
	$(F08) $(CFLAGS) test_abscab.f08 -c
	$(F08) $(CFLAGS) mod_kinds.o mod_cel.o mod_compsum.o abscab.o mod_testutil.o test_abscab.o -o test_abscab

demo_abscab: abscab.o mod_testutil.o demo_abscab.f08
	$(F08) $(CFLAGS) demo_abscab.f08 -c
	$(F08) $(CFLAGS) mod_kinds.o mod_cel.o mod_compsum.o abscab.o mod_testutil.o demo_abscab.o -o demo_abscab
