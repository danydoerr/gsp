#include <algorithm>
#include <map>
#include <set>
#include <iostream>
#include <utility>

#include "Util.h"
#include "IMP.h"


void IMP(std::vector<Region>& protoAtoms,
	std::vector<WasteRegion>& wasteRegions,
	const std::vector<std::vector<AlignmentRecord *>>& buckets,
	unsigned int bucketSize, unsigned int minLength, double epsilon,
	const std::chrono::time_point<std::chrono::high_resolution_clock> start,
	unsigned int numThreads) {
	
	auto startIMP = std::chrono::high_resolution_clock::now();
	
	int iterationCount = 0;
	while (true) {
		#pragma omp declare reduction (merge : std::vector<Region> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
		std::vector<Region> newRegions;
		#pragma omp parallel for num_threads(numThreads) reduction(merge: newRegions)
		for (size_t i = 0; i < protoAtoms.size(); i++) { // iterate over all current atoms
			Region* atom = &protoAtoms[i];
			unsigned long bucketIdx = atom->getMiddlePos() / bucketSize;
			auto alns = &buckets[bucketIdx]; // get all alignments that contain middlePos
			std::vector<Region> intervals; // waste region set W
			for (auto aln : *alns) { // iterate over all alignments covering the atom
				if (aln->tStart > atom->first || aln->tEnd < atom->last) continue; // skip alns that don't cover atom
				Region mappedRegion = mapAtomThroughAln(*atom, *aln);
				auto regionFirst = binSearchRegion(mappedRegion.first, wasteRegions);
				auto regionLast = binSearchRegion(mappedRegion.last, wasteRegions);
				for (auto j = regionFirst; j <= regionLast; j++) { // iterate over waste regions in mappedRegion
					WasteRegion* currentRegion = &wasteRegions[j];
					if (mappedRegion.first > currentRegion->last || currentRegion-> first > mappedRegion.last) continue;
					// map waste region back to atom
					auto inverseRegionFirst = mapBreakpoint(currentRegion->first, *(aln->sym));
					auto inverseRegionLast = mapBreakpoint(currentRegion->last, *(aln->sym));
					// skip if inversely mapped region does not overlap atom
					if ((inverseRegionFirst < atom->first && inverseRegionLast < atom->first)
						|| (inverseRegionFirst > atom->last && inverseRegionLast > atom->last))
						continue;
					// else push region to interval list
					if (inverseRegionFirst > inverseRegionLast)
						std::swap(inverseRegionFirst, inverseRegionLast);
					Region inverselyMappedRegion(inverseRegionFirst, inverseRegionLast);
					// clip ends to atom
					if (inverselyMappedRegion.first < atom->first) inverselyMappedRegion.first = atom->first;
					if (inverselyMappedRegion.last > atom->last) inverselyMappedRegion.last = atom->last;
					intervals.push_back(inverselyMappedRegion);
				} // end of iteration over waste in mappedRegion
			} // end of iteration over alignments containing middlepos
			// add waste regions at ends of atom
			intervals.push_back(Region(atom->first, atom->first));
			intervals.push_back(Region(atom->last, atom->last));
			std::sort(intervals.begin(), intervals.end()); // sorting before removing duplicates
			intervals.erase(std::unique(intervals.begin(), intervals.end()), intervals.end()); // remove duplicates

			// create waste region set set W_new from W
			std::vector<Region> covering, notCovering, newWasteRegions;
			partitionCoveringRegion(intervals, minLength, covering, notCovering);
			createNewWasteRegions(notCovering, covering, epsilon, minLength, atom->first, newWasteRegions);
			// add W_new to all new regions
			newRegions.insert(newRegions.end(), newWasteRegions.begin(), newWasteRegions.end());
		}
		wasteRegions.insert(wasteRegions.end(), newRegions.begin(), newRegions.end());
		consolidateRegions(wasteRegions, minLength); // join new and old waste regions
		std::vector<Region> newAtoms;
		atomsFromWaste(wasteRegions, newAtoms);
		if (!areDifferent(protoAtoms, newAtoms)) break; // stop if there is no improvement
		protoAtoms = newAtoms;
		std::cerr << "INFO: " << wasteRegions.size() << " waste regions after IMP iteration "
			<< ++iterationCount << ".";
		shoutTime(start);
	}
	std::cerr << "IMP algorithm done.";

	auto endIMP = std::chrono::high_resolution_clock::now();
	auto timeIMP = std::chrono::duration_cast<std::chrono::milliseconds>(endIMP - startIMP).count();
	std::cerr << " Algorithm time: " << timeIMP << " milliseconds.";
	
	shoutTime(start);
}

void fillBuckets(std::deque<AlignmentRecord *>& alns, unsigned int bucketSize,
	std::vector<std::vector<AlignmentRecord *>>& result) {
	unsigned int firstBucket, lastBucket;
        for (auto bucket : result)
            bucket.reserve(bucketSize); // preallocate vector of the necessary size
	for (auto alnPtr : alns) {
		firstBucket = alnPtr->tStart / bucketSize;
		lastBucket = alnPtr->tEnd / bucketSize;
		for (auto i = firstBucket; i <= lastBucket; i++)
			result[i].push_back(alnPtr);
	}
}

unsigned int binSearch_tStarts(unsigned long x, const AlignmentRecord& aln) {
	unsigned int result = std::distance(aln.begin_tStarts(), std::upper_bound(aln.begin_tStarts(), aln.end_tStarts(), x));
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
	auto idx = binSearch_tStarts(bpPosition, aln);
	unsigned int result;
	unsigned long dist = (bpPosition >= aln.get_tStarts(idx)) ? bpPosition - aln.get_tStarts(idx) : 0;
	if (dist > aln.blockSizes[idx]) dist = aln.blockSizes[idx];
	if (aln.strand == '+')
		result = aln.get_qStarts(idx) + dist;
	else
		result = aln.get_qStarts(idx) - dist;
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

/* Calculates the optimal cost for a new waste region set with waste regions at pos and in closestLeftRegion.
Best results for each position are stored for dynamic programming. */
void dpFindOptimal(Region closestLeftRegion, std::map<unsigned long, dpPosition>& allPositions,
	unsigned long pos, double epsilon, unsigned int minLength) {
	std::vector<dpStats> positionCost; // contains min cost for each pos & position which achieved it
	for (auto l = closestLeftRegion.first; l <= closestLeftRegion.last; l++) { // iterate over P(j,k)
		if ((pos - l) < minLength) // join waste regions
			positionCost.push_back(dpStats(allPositions.find(l)->second.cost + pos - l, true, l));
		else {
			bool alignedToWaste = false;
			for (auto i : allPositions.find(l)->second.notCoveringIds)
				for (auto j : allPositions.find(pos)->second.notCoveringIds)
					if (i == j) alignedToWaste = true;
			if (alignedToWaste)  // join waste regions
				positionCost.push_back(dpStats(allPositions.find(l)->second.cost + pos - l, true, l));
			else // don't join - create new atom in between
				positionCost.push_back(dpStats(allPositions.find(l)->second.cost + epsilon, false, l));

			// now check covering regions
			alignedToWaste = false;
			for (auto i : allPositions.find(l)->second.coveringIds)
				for (auto j : allPositions.find(pos)->second.coveringIds)
					if (i == j) alignedToWaste = true;
			if (alignedToWaste)
				positionCost.push_back(dpStats(allPositions.find(l)->second.cost + pos - l, true, l));
			else
				positionCost.push_back(dpStats(allPositions.find(l)->second.cost + epsilon, false, l));
		}
	}
	dpStats optimal = *std::min_element(positionCost.rbegin(), positionCost.rend(),
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
	std::vector<bool> tmpRegionBools;
        bool is_first = true;
	while (currentPos >= atomFirst) {
		if (is_first) {
                        result.push_back(Region(currentPos, currentPos));
			tmpRegionBools.push_back(posData->dist);
                        is_first = false;
                }
		else {
                        bool tmpBool = tmpRegionBools.back();
                        Region *tmpRegion = &result.back();
			if (tmpBool) {
				tmpRegion->first = currentPos;
                                tmpRegionBools[tmpRegionBools.size()-1] = posData->dist;
			}
			else {
                                tmpRegionBools.push_back(posData->dist);
				result.push_back(Region(currentPos, currentPos));
			}
		}
		if (!currentPos) break; // atom starts at 0
		currentPos = posData->prev;
		posData = &(allPositions.find(currentPos)->second);
	}
}

void createNewWasteRegions(const std::vector<Region>& notCovering, const std::vector<Region>& covering,
	double epsilon, unsigned int minLength, unsigned long atomStart, std::vector<Region>& result) {
	std::set<unsigned long> nonCovPos; // positions in noncovering
	std::map<unsigned long, dpPosition> allPositions; // positions in either vector

	// collect positions of nonCovering intervals
	for (size_t i = 0; i < notCovering.size(); i++) {
		const Region* curRegion = &notCovering[i];
		for (auto pos = curRegion->first; pos <= curRegion->last; pos++) {
			nonCovPos.insert(pos);
			allPositions.insert(std::pair<unsigned long, dpPosition>(pos, dpPosition(i))).
				first->second.notCoveringIds.push_back(i);
		}
	}
	// add positions also contained in covering regions
	for (size_t i = 0; i < covering.size(); i++) {
		const Region* curRegion = &covering[i];
		for (auto pos = curRegion->first; pos <= curRegion->last; pos++)
			if (allPositions.count(pos))
				allPositions.find(pos)->second.coveringIds.push_back(i);
	}

	std::set<unsigned int> *currentShortIntervals = new std::set<unsigned int>(), *lastShortIntervals = new std::set<unsigned int>();
	unsigned int lastFinishedIdx = 0;
	for (auto pos : nonCovPos) { // iterate over all viable positions i from left to right
		auto position = allPositions.find(pos);
		for (auto i : position->second.notCoveringIds)
			currentShortIntervals->insert(i);
		if (pos == *(nonCovPos.begin())) continue; // only init for first (leftmost) position
		// get ID of rightmost region not containing pos but left of pos
		for (auto previous : *lastShortIntervals)
			if (currentShortIntervals->find(previous) == currentShortIntervals->end()) // not in set
				lastFinishedIdx = previous;
		dpFindOptimal(notCovering[lastFinishedIdx], allPositions, pos, epsilon, minLength);
		delete lastShortIntervals;
                lastShortIntervals = currentShortIntervals; //lastShortIntervals = currentShortIntervals;
		currentShortIntervals = new std::set<unsigned int>(); //currentShortIntervals.clear();
	}
	dpTraceBack(allPositions, notCovering.back().last, atomStart, result);
        delete lastShortIntervals; delete currentShortIntervals;
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
