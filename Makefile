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
debug: pbj

clean:
	rm -f src/*.o pbj
	rm -rf tests/out

test: pbj
	./tests/run_tests.sh

bench: pbj
	./bench/throughput.sh

.PHONY: all clean test bench native debug
