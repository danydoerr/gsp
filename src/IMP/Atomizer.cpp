#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <vector>
#include <algorithm>

#include "AlignmentRecord.h"
#include "InputParser.h"
#include "Breakpoints.h"
#include "IMP.h"
#include "Classify.h"

void IMP(std::vector<Region>& , std::vector<WasteRegion>&,
	const std::vector<std::vector<std::shared_ptr<AlignmentRecord>>>&,
	unsigned int, unsigned int, double,
	const std::chrono::time_point<std::chrono::high_resolution_clock>,
	unsigned int);
void printResult(const std::vector<WasteRegion>&,
	const std::vector<int>&, 
	const std::map<std::string, unsigned long>&);
void shoutTime(const std::chrono::time_point<std::chrono::high_resolution_clock>);

int main(int argc, char** argv) {
	// only reason the following vars are not const is for cmd arg parsing
	char* pslPath;
	unsigned int maxGapLength = 13, minAlnLength = 13, minLength = 250, bucketSize = 1000, numThreads = 1;
	float minAlnIdentity = 0.8f;
	parseCmdArgs(argc, argv, pslPath, minLength, maxGapLength, minAlnLength, minAlnIdentity, bucketSize, numThreads);

	// init maps and vectors
	std::map<std::string, unsigned long> speciesStarts; // maps species name to their starting position in concatenated string
	std::vector<unsigned long> speciesBoundaries; // contains starting positions in concatenated sequence
	std::vector<std::shared_ptr<AlignmentRecord>> alignments;
	std::vector<Breakpoint> breakPoints;
	std::vector<WasteRegion> wasteRegions;
	std::vector<Region> protoAtoms;

	std::cerr << "Starting: " << pslPath << " is used for segmentation with parameters:\n"
		<< "minLength: " << minLength << ", minIdent: " << minAlnIdentity * 100 << ", maxGap: "
		<< maxGapLength << ", minAlnLength: " << minAlnLength
		<<  ", bucketSize: " << bucketSize << std::endl;
	auto start = std::chrono::high_resolution_clock::now();
	speciesStarts = { {"$", 0} };
	parsePsl(pslPath, speciesStarts, maxGapLength, minAlnLength, minAlnIdentity, alignments);
	for (auto i : speciesStarts) speciesBoundaries.push_back(i.second);
	std::cerr << "INFO: PSL parsing done, considering " << alignments.size() << " alignments between "
		<< speciesStarts.size() - 1 << " sequences.";
	shoutTime(start);
	std::vector<std::vector<std::shared_ptr<AlignmentRecord>>>
		buckets((speciesStarts.find("$")->second / bucketSize) + 1); // reserve with appropiate size
	fillBuckets(alignments, bucketSize, buckets);
	std::cerr << "INFO: Filled " << buckets.size() << " buckets.";
	shoutTime(start);
	const double epsilon = 1 / (static_cast<double>(bucketSize)*buckets.size());
	initBreakpoints(alignments, speciesBoundaries, breakPoints);
	createWaste(breakPoints, minLength, wasteRegions);
	atomsFromWaste(wasteRegions, protoAtoms);
	std::cerr << "INFO: Created " << wasteRegions.size() << " initial waste regions from initial breakpoints.";
	shoutTime(start);
	IMP(protoAtoms, wasteRegions, buckets, bucketSize, minLength, epsilon, start, numThreads);
	std::vector<int> classes;
	int nrClasses = 0;
	classify(wasteRegions, buckets, bucketSize, minAlnIdentity, classes, nrClasses);
	std::cerr << "Put " << wasteRegions.size() - 1 << " atoms in " << nrClasses << " classes. "
		<< "Printing result." << std::endl;
	shoutTime(start);
	printResult(wasteRegions, classes, speciesStarts);
	return EXIT_SUCCESS;
}

void IMP(std::vector<Region>& protoAtoms,
	std::vector<WasteRegion>& wasteRegions,
	const std::vector<std::vector<std::shared_ptr<AlignmentRecord>>>& buckets,
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

void printResult(const std::vector<WasteRegion> &regions, const std::vector<int> &classes,
	const std::map<std::string, unsigned long> &speciesStarts) {
	// flipped is speciesStarts sorted by position
	std::map<unsigned long, std::string> flipped;
	for (auto specStart : speciesStarts)
		flipped[specStart.second] = specStart.first;
	std::vector<std::string> names;
	std::vector<unsigned long> starts;
	for (auto specStart : flipped) {
		starts.push_back(specStart.first);
		names.push_back(specStart.second);
	}

	std::cout << "#name\tatom_nr\tclass\tstrand\tstart\tend" << "\n"; // header line
	for (size_t i = 0; i + 1 < regions.size(); i++) {
		auto j = binSearch(regions[i].last, starts);
		auto move = starts[j];
		auto start = (regions[i].last > move) ? regions[i].last - move : 0;
		auto end = regions[i + 1].first - move;
		if (end > starts[j + 1]) end = starts[j + 1];
		auto classNr = abs(classes[i]);
		auto strand = (classes[i] > 0) ? '+' : '-';
		std::cout << names[j] << "\t" << i+1 << "\t" << classNr << "\t" << strand << "\t"
			<< start << "\t" << end << "\n";
	}
}

void shoutTime(const std::chrono::time_point<std::chrono::high_resolution_clock> start) {
	auto end = std::chrono::high_resolution_clock::now();
	auto diff = end - start;
	auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();
	std::cerr << " Time passed since start: " << ms << " milliseconds." << std::endl;
}
