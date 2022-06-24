
F08?=gfortran
CFLAGS=-std=f2008

ABSCAB_DIR=../../main/fortran

.PHONY: all test clean test_programs # demo_programs

all: test_programs # demo_programs

test_programs: test_cel # test_abscab

#demo_programs: demo_abscab

test: test_programs
	@ ./test_cel
#	@ ./test_abscab

#demo: demo_programs
#	@ ./demo_abscab

clean:
	rm -f test_cel #test_abscab demo_abscab

mod_kinds: $(ABSCAB_DIR)/mod_kinds.f08
	$(F08) $(CFLAGS) $(ABSCAB_DIR)/mod_kinds.f08 -c

mod_cel: mod_kinds $(ABSCAB_DIR)/mod_cel.f08
	$(F08) $(CFLAGS) $(ABSCAB_DIR)/mod_cel.f08 -c

test_cel: mod_cel mod_kinds test_cel.f08
	$(F08) $(CFLAGS) test_cel.f08 -c
	$(F08) mod_kinds.o mod_cel.o test_cel.o -o test_cel
#
#test_abscab: $(ABSCAB_DIR)/abscab.f08 test_abscab.f08
#	$(FC) $(CFLAGS) test_abscab.f08 -I$(ABSCAB_DIR) -o test_abscab
#
#demo_abscab: $(ABSCAB_DIR)/abscab.f08 demo_abscab.f08
#	$(FC) $(CFLAGS) demo_abscab.f08 -I$(ABSCAB_DIR) -o demo_abscab
