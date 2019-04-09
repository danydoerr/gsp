CC = g++
CFLAGS = -std=c++14 -Wall -O3 -fopenmp
DEBUG_FLAGS = -std=c++14 -Wall -g -fopenmp

all: bin debug

bin:
	$(CC) $(CFLAGS) *.cpp -o atomizer

debug:
	$(CC) $(DEBUG_FLAGS) *.cpp -o atomizer_debug

clean:
	rm -f atomizer atomizer_debug
