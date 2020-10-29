CC = g++
CFLAGS = -O3 -Wall -Wextra
SRCDIR = ./src
BINDIR = ./bin

default: mpi

mpi: $(SRCDIR)/n-body.cpp BinDir
	mpicxx $(CFLAGS) -o $(BINDIR)/n-body $(SRCDIR)/n-body.cpp

debug: CFLAGS += -DDEBUG
debug: default

BinDir:
	mkdir -p $(BINDIR)

clean:
	rm -rf $(BINDIR)
