###############################################################################
# Makefile for Project, High Performance Programming
###############################################################################

CFLAGS = -g -O3 -Wall -fopenmp
LIBS = -lm

BIN = MC

CC = gcc

all: $(BIN)

$(BIN): MC_parallel.c MC.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
	
clean:
	$(RM) $(BIN)

cleanall:
	rm -f *.csv