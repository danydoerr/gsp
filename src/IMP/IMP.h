#pragma once

#include "AlignmentRecord.h"

/* Returns index of the last element in xList that is <= x.
If all elements in xList are > x, result is 0. Expects xList to be sorted ascending. */
unsigned int binSearch(unsigned long x, const std::vector<unsigned long>& xList);

/* Returns index of the last element in bpList whose starting position is <= x.
If there are none, result is 0. Expects bpList to be sorted ascending. */
unsigned int binSearchRegion(unsigned long x, const std::vector<WasteRegion>& bpList);

/* Maps input breakpoint from alignment query to alignment target. */
unsigned int mapBreakpoint(unsigned long bpPosition, const AlignmentRecord& aln);

/* Maps an atom to the target of an alignment covering that atom.
This means it returns a region that is aligned to the input atom. */
Region mapAtomThroughAln(const Region& atom, const AlignmentRecord& aln);

/* Paritions input regions in two set, one containing the regions that cover other regions,
the other one containing the ones that don't.
Input must be sorted according to operator < in Region. */
void partitionCoveringRegion(const std::vector<Region>& input, unsigned int minLength,
	std::vector<Region>& covering, std::vector<Region>& notCovering);

/* Creates a new optimal set of waste region from notCovering and covering
via dynamic programming and stores it in result. */
void createNewWasteRegions(const std::vector<Region>& notCovering, const std::vector<Region>& covering,
	double epsilon, unsigned int minLength, unsigned long atomStart, std::vector<Region>& result);

/* Joins newly added waste regions with older ones. */
void consolidateRegions(std::vector<WasteRegion> &regions, unsigned int minLength);

/* Checks if both vectors contain the same elements.
Expects both input vectors to be sorted in the same way, e.g. by atom length. */
bool areDifferent(std::vector<Region> &first, std::vector<Region> &second);