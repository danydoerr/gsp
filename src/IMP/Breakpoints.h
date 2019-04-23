#pragma once
#include <vector>
#include <memory>
#include "AlignmentRecord.h"

/* Creates initial breakpoints from alignment and species boundaries and stores them in result. */
void initBreakpoints(const std::vector<AlignmentRecord *>& alns,
	const std::vector<unsigned long>& speciesBounds,
	std::vector<Breakpoint>& result);

/* Stores a list of Regions in result, created from input breakpoints.
The result will be sorted by position. Expects input breakpoints to be sorted by position as well. */
void createWaste(const std::vector<Breakpoint>& breakpoints, unsigned int minLength,
	std::vector<WasteRegion>& result);

/* Creates atoms as regions in between waste regions and stores them in result.
The result will be sorted by length, ascending. */
void atomsFromWaste(std::vector<WasteRegion>& wasteRegions, std::vector<Region>& result);