#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <vector>
#include <deque>
#include <algorithm>

#include "AlignmentRecord.h"
#include "InputParser.h"
#include "Breakpoints.h"
#include "IMP.h"
#include "Classify.h"
#include "Util.h"

int main(int argc, char** argv) {
	// only reason the following vars are not const is for cmd arg parsing
	unsigned int maxGapLength, minAlnLength, minLength, bucketSize, numThreads;
	float minAlnIdentity;
        
        InputParser parser;
        parser.parseCmdArgs(argc, argv);
        parser.getCmdLineArgs(minLength, maxGapLength, minAlnLength, minAlnIdentity, bucketSize, numThreads);

	// init maps and vectors
	std::map<std::string, unsigned long> speciesStarts; // maps species name to their starting position in concatenated string
	std::vector<unsigned long> speciesBoundaries; // contains starting positions in concatenated sequence
	std::deque<AlignmentRecord *> alignments;
	std::vector<Breakpoint> breakPoints;
	std::vector<WasteRegion> wasteRegions;
	std::vector<Region> protoAtoms;

	std::cerr << "Starting with parameters:\n"
		<< "minLength: " << minLength << ", minIdent: " << minAlnIdentity * 100 << ", maxGap: "
		<< maxGapLength << ", minAlnLength: " << minAlnLength
		<<  ", bucketSize: " << bucketSize
                <<  ", numThreads: " << numThreads << std::endl;
	auto start = std::chrono::high_resolution_clock::now();
	speciesStarts = { {"$", 0} };
        
        try{
            parser.parsePsl(speciesStarts, alignments);
        }catch(const std::exception &e){
            std::cerr << e.what() << std::endl;
            throw;
        }
	
	for (auto i : speciesStarts) speciesBoundaries.push_back(i.second);
	std::cerr << "INFO: PSL parsing done, considering " << alignments.size() << " alignments between "
		<< speciesStarts.size() - 1 << " sequences.";
	shoutTime(start);
	std::vector<std::vector<AlignmentRecord *>>
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
        for (auto aln : alignments)
            delete aln;
	return EXIT_SUCCESS;
}
