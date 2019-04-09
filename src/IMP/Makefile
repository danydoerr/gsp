CC = g++
CFLAGS = -std=c++14 -Wall -g -O2

all: bin debug

bin:
	$(CC) $(CFLAGS) *.cpp -o atomizer -fopenmp -O3

debug:
	$(CC) $(CFLAGS) *.cpp -g -o atomizer_debug -fopenmp

