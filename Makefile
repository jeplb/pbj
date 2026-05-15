CC      ?= cc
CFLAGS  ?= -O3 -Wall -Wextra -Wno-unused-parameter -std=c11 -pthread
LDFLAGS ?=
LDLIBS  ?= -lz -lpthread -lm

SRCS = src/main.c src/fastq.c src/revcomp.c src/overlap.c src/score.c src/stats.c src/merge.c src/threads.c src/adapter.c
OBJS = $(SRCS:.c=.o)
HDRS = src/pbj.h src/fastq.h src/revcomp.h src/overlap.h src/score.h src/stats.h src/merge.h src/threads.h src/adapter.h src/kseq.h

all: pbj

pbj: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)

%.o: %.c $(HDRS)
	$(CC) $(CFLAGS) -c -o $@ $<

native: CFLAGS += -march=native
native: pbj

debug: CFLAGS = -O0 -g -Wall -Wextra -std=c11 -pthread -fsanitize=address,undefined
debug: LDFLAGS += -fsanitize=address,undefined
debug: pbj tests/unit/test_unit tests/unit/test_edge_cases

UNIT_TEST_OBJS = src/fastq.o src/revcomp.o src/overlap.o src/score.o \
                 src/stats.o src/merge.o src/adapter.o
UNIT_TEST_HDRS = $(HDRS) tests/unit/test_helpers.h

tests/unit/test_unit: tests/unit/test_unit.c $(UNIT_TEST_OBJS) $(UNIT_TEST_HDRS)
	$(CC) $(CFLAGS) -o $@ tests/unit/test_unit.c $(UNIT_TEST_OBJS) \
		$(LDFLAGS) $(LDLIBS)

tests/unit/test_edge_cases: tests/unit/test_edge_cases.c $(UNIT_TEST_OBJS) $(UNIT_TEST_HDRS)
	$(CC) $(CFLAGS) -o $@ tests/unit/test_edge_cases.c $(UNIT_TEST_OBJS) \
		$(LDFLAGS) $(LDLIBS)

clean:
	rm -f src/*.o pbj tests/unit/test_unit tests/unit/test_edge_cases
	rm -rf tests/out

test: pbj unit-test
	./tests/run_tests.sh

unit-test: tests/unit/test_unit tests/unit/test_edge_cases
	./tests/unit/test_unit
	./tests/unit/test_edge_cases

bench: pbj
	./bench/throughput.sh

.PHONY: all clean test unit-test bench native debug
