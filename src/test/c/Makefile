
CC?=gcc

ABSCAB_DIR=../../main/c

.PHONY: all test clean

all: test_cel test_abscab
	
test: all
	@ ./test_cel
	@ ./test_abscab

clean:
	rm -f test_cel test_abscab

test_cel: $(ABSCAB_DIR)/cel.h util.h test_cel.c
	$(CC) test_cel.c -I$(ABSCAB_DIR) -o test_cel -lm 

test_abscab: $(ABSCAB_DIR)/abscab.h util.h test_abscab.c
	$(CC) test_abscab.c -I$(ABSCAB_DIR) -o test_abscab -lm
