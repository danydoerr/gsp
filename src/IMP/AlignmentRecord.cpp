#include <iostream>
#include "AlignmentRecord.h"

AlignmentRecord::AlignmentRecord(char strand,
	unsigned long qStart, unsigned long qEnd,
	unsigned long tStart, unsigned long tEnd,
	unsigned int blockCount, std::vector<unsigned int> blockSizes,
	std::vector<unsigned long> qStarts, std::vector<unsigned long> tStarts)
	: strand(strand), qStart(qStart), qEnd(qEnd), tStart(tStart), tEnd(tEnd),
	blockCount(blockCount), blockSizes(blockSizes), qStarts(qStarts), tStarts(tStarts), sym(nullptr) {}

void AlignmentRecord::printRecord() const {
	std::cout << "Strand: " << strand << "\n";
	std::cout << "qStart: " << qStart << "\n";
	std::cout << "qEnd: " << qEnd << "\n";
	std::cout << "tStart: " << tStart << "\n";
	std::cout << "tEnd: " << tEnd << "\n";
	std::cout << "blockCount: " << blockCount << "\n";
	std::cout << "blockSizes: ";
	for (auto i : blockSizes) std::cout << i << ",";
	std::cout << std::endl;
	std::cout << "qStarts: ";
	for (auto i : qStarts) std::cout << i << ",";
	std::cout << "\n";
	std::cout << "tStarts: ";
	for (auto i : tStarts) std::cout << i << ",";
	std::cout << "\n";
	if (sym != nullptr)
		std::cout << "sym hast tStart " << sym->tStart << " and tEnd " << sym->tEnd <<
		". This ones tStart according to sym is " << sym->sym->tStart << "\n";
}

AlignmentRecord *AlignmentRecord::revert() const {
	// all it really does is swap query and target
	std::vector<unsigned int> newBlockSizes;
        std::vector<unsigned long> newQStarts, newTStarts;
	if (strand == '+') {
		newBlockSizes = blockSizes;
		newQStarts = tStarts;
		newTStarts = qStarts;
	}
	else { // reverse strand - revert order, and make endpoints startpoints
		for (int i = qStarts.size() - 1; i >= 0; i--) {
			newQStarts.push_back(tStarts[i] + blockSizes[i]);
			newTStarts.push_back(qStarts[i] - blockSizes[i]);
			newBlockSizes.push_back(blockSizes[i]);
		}
	}
	return new AlignmentRecord(strand, tStart, tEnd, qStart, qEnd,
		blockCount, newBlockSizes, newQStarts, newTStarts);
}

Breakpoint::Breakpoint(unsigned long position)
: position(position) {}

WasteRegion::WasteRegion(unsigned long pos)
: Region(pos,pos) {}

WasteRegion::WasteRegion(Region atom)
: Region(atom.first, atom.last) {}

Region::Region(unsigned long first, unsigned long last)
: first(first), last(last) {}

dpPosition::dpPosition(unsigned int idx)
: idx(idx), cost(0.0), dist(false), prev(0) {}

dpStats::dpStats(double cost, bool dist, unsigned long prev)
: cost(cost), dist(dist), prev(prev) {}