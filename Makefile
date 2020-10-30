CC = g++
CFLAGS = -O3 -Wall -Wextra -Wno-c++11-extensions
SRCDIR = ./src
BINDIR = ./bin

default: mpi

mpi: $(SRCDIR)/n-body.cpp bin-dir
	mpicxx $(CFLAGS) -o $(BINDIR)/n-body $(SRCDIR)/n-body.cpp

debug: CFLAGS += -DDEBUG
debug: default

bin-dir:
	mkdir -p $(BINDIR)

.PHONY: clean
clean:
	rm -rf $(BINDIR)
