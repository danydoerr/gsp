#include <iostream>
#include <algorithm>
#include "Breakpoints.h"

void initBreakpoints(const std::deque<AlignmentRecord *>& alns,
	const std::vector<unsigned long>& speciesBounds,
	std::vector<Breakpoint>& result) {
	for (auto bp : speciesBounds)
		result.push_back(Breakpoint(bp));
	for (auto aln : alns) {
		result.push_back(Breakpoint(aln->tStart));
		result.push_back(Breakpoint(aln->tEnd));
	}
	std::sort(result.begin(), result.end()); // sort breakpoints by position
	auto last = std::unique(result.begin(), result.end()); // remove duplicate breakpoints
	result.erase(last, result.end());
}

void createWaste(const std::vector<Breakpoint>& breakpoints, unsigned int minLength,
	std::vector<WasteRegion>& result) {
	if (breakpoints.empty()) {
		std::cerr << "ERROR: Got empty breakpoint list when trying to create regions.";
		exit(EXIT_FAILURE);
	}
	result.push_back(WasteRegion(breakpoints[0].position));
	for (size_t i = 1; i < breakpoints.size(); i++) {
		auto prev = &(result.back());
		unsigned long distance = breakpoints[i].position - prev->last;
		if (distance <= minLength) // too close for atom to be in between
			prev->last = breakpoints[i].position;
		else// distance > minLength, create new region
			result.push_back(WasteRegion(breakpoints[i].position));
	}
}

void atomsFromWaste(std::vector<WasteRegion>& wasteRegions, std::vector<Region>& result) {
	for (size_t i = 0; i < wasteRegions.size() - 1; i++)
		result.push_back(Region(wasteRegions[i].last,wasteRegions[i+1].first));
}