#include <algorithm>
#include <map>
#include <set>
#include <iostream>
#include <utility>

#include "IMP.h"

unsigned int binSearch(unsigned long x, const std::vector<unsigned long>& xList) {
	unsigned int result = std::distance(xList.begin(), std::upper_bound(xList.begin(), xList.end(), x));
	if (result == 0) return result;
	else return result - 1;
}

unsigned int binSearchRegion(unsigned long x, const std::vector<WasteRegion>& bpList) {
	const WasteRegion tmp(x);
	unsigned int result = std::distance(bpList.begin(), std::upper_bound(bpList.begin(), bpList.end(), tmp,
		[](WasteRegion lbp, WasteRegion rbp) {return lbp.first < rbp.first; }));
	if (result == 0) return result;
	else return result - 1;
}

unsigned int mapBreakpoint(unsigned long bpPosition, const AlignmentRecord& aln) {
	auto idx = binSearch(bpPosition, aln.tStarts);
	unsigned int result;
	unsigned long dist = (bpPosition >= aln.tStarts[idx]) ? bpPosition - aln.tStarts[idx] : 0;
	if (dist > aln.blockSizes[idx]) dist = aln.blockSizes[idx];
	if (aln.strand == '+')
		result = aln.qStarts[idx] + dist;
	else
		result = aln.qStarts[idx] - dist;
	return result;
}

Region mapAtomThroughAln(const Region& atom, const AlignmentRecord& aln) {
	auto firstMapped = mapBreakpoint(atom.first, aln);
	auto lastMapped = mapBreakpoint(atom.last, aln);
	if (firstMapped <= lastMapped) return Region(firstMapped, lastMapped);
	else return Region(lastMapped, firstMapped);
}

void partitionCoveringRegion(const std::vector<Region>& input, unsigned int minL,
	std::vector<Region>& covering, std::vector<Region>& notCovering) {
	int minLength = static_cast<signed int>(minL);
	for (size_t i = 0; i < input.size(); i++) {
		const Region* currentRegion = &input[i];
		long lastShortStart, lastLongStart;

		if (notCovering.empty()) lastShortStart = 0 - minLength - 2;
		else lastShortStart = notCovering.back().first;
		if (covering.empty()) lastLongStart = 0 - minLength - 2;
		else lastLongStart = covering.back().first;

		if (lastLongStart >= 0 && static_cast<unsigned long>(lastLongStart) >= currentRegion->first) {
			while (lastLongStart >= 0 && static_cast<unsigned long>(lastLongStart) >= currentRegion->first) {
				if (!covering.empty()) covering.pop_back();
				if (covering.empty()) lastLongStart = 0 - minLength - 2;
				else lastLongStart = covering.back().first;
			}
			covering.push_back(*currentRegion);
		} else {
			if (lastShortStart >= 0 && static_cast<unsigned long>(lastShortStart) >= currentRegion->first)
				covering.push_back(*currentRegion);
			else notCovering.push_back(*currentRegion);
		}
	}
}

/* Calculates the optimal cost for a new waste region set with waste regions at pos and in region.
Previous best results for positions are stored for dynamic programming. */
void dpFindOptimal(Region region, std::map<unsigned long, dpPosition>& allPositions,
	unsigned long pos, double epsilon, unsigned int minLength) {
	std::vector<dpStats> tmp;
	for (auto l = region.first; l <= region.last; l++) { // iterate over P(j,k)
		if ((pos - l) < minLength) // join waste regions
			tmp.push_back(dpStats(allPositions.find(l)->second.cost + pos - l, true, l));
		else {
			bool alignedToWaste = false;
			for (auto i : allPositions.find(l)->second.notCoveringIds) {
				for (auto j : allPositions.find(pos)->second.notCoveringIds) {
					if (i == j) alignedToWaste = true;
					goto endloop;
				}
			}
		endloop:
			if (alignedToWaste) tmp.push_back(dpStats(allPositions.find(l)->second.cost + pos - l, true, l));
			else tmp.push_back(dpStats(allPositions.find(l)->second.cost + epsilon, false, l));

			// now check covering regions
			alignedToWaste = false;
			for (auto i : allPositions.find(l)->second.coveringIds) {
				for (auto j : allPositions.find(pos)->second.coveringIds) {
					if (i == j) alignedToWaste = true;
					goto endloop2;
				}
			}
		endloop2:
			if (alignedToWaste) tmp.push_back(dpStats(allPositions.find(l)->second.cost + pos - l, true, l));
			else tmp.push_back(dpStats(allPositions.find(l)->second.cost + epsilon, false, l));
		}
	}
	dpStats optimal = *std::min_element(tmp.begin(), tmp.end(),
		[](dpStats a, dpStats b) {return a.cost < b.cost; });
	dpPosition* toChange = &allPositions.find(pos)->second;
	toChange->cost = optimal.cost;
	toChange->dist = optimal.dist;
	toChange->prev = optimal.prev;
}

/* After the cost of an optimal solution is computed, the optimal set for that solution
is created by tracing back the stored positions. The optimal set will be stored in result. */
void dpTraceBack(std::map<unsigned long, dpPosition>& allPositions,
	unsigned long lastPos, unsigned long atomFirst, std::vector<Region>& result) {
	unsigned long currentPos = allPositions.find(lastPos)->second.prev;
	dpPosition *posData = &(allPositions.find(currentPos)->second);
	std::vector<std::pair<bool, Region>> tmpRegions;
	while (currentPos >= atomFirst) {
		if (tmpRegions.empty())
			tmpRegions.push_back(std::make_pair(posData->dist, Region(currentPos, currentPos)));
		else {
			std::pair<bool, Region> *tmp = &tmpRegions.back();
			if (tmp->first) {
				tmp->second.first = currentPos;
				tmp->first = posData->dist;
			}
			else {
				tmpRegions.push_back(std::make_pair(posData->dist, Region(currentPos, currentPos)));
			}
		}
		if (!currentPos) break; // atom starts at 0
		currentPos = posData->prev;
		posData = &(allPositions.find(currentPos)->second);
	}
	for (auto i : tmpRegions)
		result.push_back(i.second);
}

void createNewWasteRegions(const std::vector<Region>& notCovering, const std::vector<Region>& covering,
	double epsilon, unsigned int minLength, unsigned long atomStart, std::vector<Region>& result) {
	std::set<unsigned long> allPos;
	// stores each position contained in both region vectors
	std::map<unsigned long, dpPosition> allPositions;
	for (size_t i = 0; i < notCovering.size(); i++) { // collect positions of nonCovering intervals
		const Region* curRegion = &notCovering[i];
		for (auto pos = curRegion->first; pos <= curRegion->last; pos++) {
			allPos.insert(pos);
			allPositions.insert(std::pair<int, dpPosition>(pos, dpPosition(i))).
				first->second.notCoveringIds.push_back(i);
		}
	}
	for (size_t i = 0; i < covering.size(); i++) { // add positions also contained in covering regions
		const Region* curRegion = &covering[i];
		for (auto pos = curRegion->first; pos <= curRegion->last; pos++)
			if (allPositions.count(pos)) // TODO there may be a mistake here
				allPositions.find(pos)->second.coveringIds.push_back(i);
	}

	std::set<unsigned int> currentShortIntervals, lastShortIntervals;
	unsigned int lastFinishedIdx = 0;
	for (auto pos : allPos) {
		auto position = allPositions.find(pos);
		for (auto i : position->second.notCoveringIds)
			currentShortIntervals.insert(i);

		if (pos == *(allPos.begin())) continue; // only do init for first element
		// get ID of rightmost region not containing pos but left of pos
		for (auto previous : lastShortIntervals)
			if (currentShortIntervals.find(previous) == currentShortIntervals.end()) // not in set
				lastFinishedIdx = previous;
		dpFindOptimal(notCovering[lastFinishedIdx], allPositions, pos, epsilon, minLength);
		lastShortIntervals = currentShortIntervals;
		currentShortIntervals.clear();
	}
	dpTraceBack(allPositions, notCovering.back().last, atomStart, result);
}

void consolidateRegions(std::vector<WasteRegion> &regions, unsigned int minLength) {
	std::vector<WasteRegion> tmp(regions);
	std::sort(tmp.begin(), tmp.end());
	regions.clear();
	size_t i = 0;
	WasteRegion currentRegion = *tmp.begin();
	while (i < tmp.size()) {
		if (i + 1 < tmp.size()) {
			WasteRegion nextRegion = tmp[i + 1];
			i++;
			if (nextRegion.first <= currentRegion.last + minLength) {
				// join regions
				currentRegion.last = std::max(currentRegion.last, nextRegion.last);
				for (auto breakpoint : nextRegion.bpPositions)
					currentRegion.bpPositions.insert(breakpoint);
			} else { // not close enough to join
				regions.push_back(currentRegion);
				currentRegion = nextRegion;
			}
		} else { // push last element region
			regions.push_back(currentRegion);
			i++;
		}
	}
}

bool areDifferent(std::vector<Region> &first, std::vector<Region> &second) {
	if (first.size() != second.size()) return true;
	// if size is equal, compare each elements positions
	for (size_t i = 0; i < first.size(); i++)
		if (first[i].first != second[i].first || first[i].last != second[i].last)
			return true;
	return false;
}