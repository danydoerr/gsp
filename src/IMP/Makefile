CC = g++
CFLAGS = -std=c++14 -Wall -fopenmp
RM_CLEAN = *.o atomizer atomizer_debug

.PHONY: debug atomizer

all: atomizer

AlignmentRecord.o: AlignmentRecord.h AlignmentRecord.cpp
	@echo "**Compiling AlignmentRecord.cpp**"
	$(CC) $(CFLAGS) -c AlignmentRecord.cpp
	@echo

Breakpoints.o: Breakpoints.h Breakpoints.cpp
	@echo "**Compiling Breakpoints.cpp**"
	$(CC) $(CFLAGS) -c Breakpoints.cpp
	@echo

Classify.o: Classify.h IMP.h Classify.cpp
	@echo "**Compiling Classify.cpp**"
	$(CC) $(CFLAGS) -c Classify.cpp
	@echo

IMP.o: Util.h IMP.h IMP.cpp
	@echo "**Compiling IMP.cpp**"
	$(CC) $(CFLAGS) -c IMP.cpp
	@echo

InputParser.o: InputParser.h AlignmentRecord.h Util.h InputParser.cpp
	@echo "**Compiling InputParser.cpp**"
	$(CC) $(CFLAGS) -c InputParser.cpp
	@echo

Util.o: Util.h Util.cpp
	@echo "**Compiling Util.cpp**"
	$(CC) $(CFLAGS) -c Util.cpp
	@echo

Atomizer.o: AlignmentRecord.h InputParser.h Breakpoints.h IMP.h Classify.h Util.h Atomizer.cpp
	@echo "**Compiling Atomizer.cpp**"
	$(CC) $(CFLAGS) -c Atomizer.cpp
	@echo

rm_obj:
	@rm -f *.o

debug: CFLAGS+=-g

debug: debug_bin

# when building debug, must remove all .o, use them, and remove them again (otherwise the not-debug bin may use them)
debug_bin: rm_obj AlignmentRecord.o Breakpoints.o Classify.o IMP.o InputParser.o Util.o Atomizer.o
	@echo "**Linking files**"
	$(CC) $(CFLAGS) AlignmentRecord.o Breakpoints.o Classify.o IMP.o InputParser.o Util.o Atomizer.o -o atomizer_debug
	@rm -f *.o
	@echo

atomizer: CFLAGS+=-O3

atomizer: atomizer_bin

atomizer_bin: AlignmentRecord.o Breakpoints.o Classify.o IMP.o InputParser.o Util.o Atomizer.o
	@echo "**Linking files**"
	$(CC) $(CFLAGS) AlignmentRecord.o Breakpoints.o Classify.o IMP.o InputParser.o Util.o Atomizer.o -o atomizer
	@echo

clean: ;
	@echo "**Cleaning files**"
	rm -f $(RM_CLEAN)
	@echo

