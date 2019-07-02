#include <iostream>
#include <cstring>
#include "AlignmentRecord.h"

// ulong2ushort, ulong2uint or nothing, depending on BLOCKS_SIZE (see AlignmentRecord.h)
/* NO NOT CHANGE THE NEXT DEFINES */
#if BLOCKS_SIZE == BLOCKS_USHORT
#define ulong2block_local_t ulong2ushort
#elif BLOCKS_SIZE == BLOCKS_UINT
#define ulong2block_local_t ulong2uint
#else
#define ulong2block_local_t 
#endif


AlignmentRecord::AlignmentRecord(char strand,
	unsigned long qStart, unsigned long qEnd,
	unsigned long tStart, unsigned long tEnd,
	unsigned int blockCount, std::vector<unsigned int> blockSizes,
	std::vector<unsigned long> qStarts, std::vector<unsigned long> tStarts)
	: strand(strand), qStart(qStart), qEnd(qEnd), tStart(tStart), tEnd(tEnd),
          blockCount(blockCount), sym(nullptr) {
        
        int i = 0;
        this->blockSizes = new block_local_t[blockCount];
        for (unsigned long x : blockSizes)
            this->blockSizes[i++] = ulong2block_local_t(x);
        
        i = 0;
        this->qStarts = new block_local_t[blockCount];
        for (unsigned long x : qStarts)
            this->qStarts[i++] = ulong2block_local_t(x - qStart); // converting to local coordinate
        
        i = 0;
        this->tStarts = new block_local_t[blockCount];
        for (unsigned long x : tStarts)
            this->tStarts[i++] = ulong2block_local_t(x - tStart); // converting to local coordinate
}

AlignmentRecord::AlignmentRecord(char strand,
	unsigned long qStart, unsigned long qEnd,
	unsigned long tStart, unsigned long tEnd,
	unsigned int blockCount, std::vector<unsigned int> blockSizes,
	std::vector<unsigned long> qStarts, std::vector<unsigned long> tStarts,
        unsigned int start_pos)
	: strand(strand), blockCount(blockCount), sym(nullptr) {
        
        unsigned int end_pos = start_pos + blockCount - 1;
        
        this->tStart = tStarts[start_pos];
	this->tEnd = tStarts[end_pos] + blockSizes[end_pos];
        if (strand == '+') {
            this->qStart = qStarts[start_pos];
            this->qEnd = qStarts[end_pos] + blockSizes[end_pos];
        } else { // strand == '-'
            this->qStart = qStarts[end_pos] - blockSizes[end_pos];
            this->qEnd = qStarts[start_pos];
        }
        
        unsigned int i;
        this->blockSizes = new block_local_t[blockCount];
        for (i = 0; i < blockCount; ++i)
            this->blockSizes[i] = ulong2block_local_t(blockSizes[start_pos + i]);
        
        this->qStarts = new block_local_t[blockCount];
        for (i = 0; i < blockCount; ++i)
            this->qStarts[i] = ulong2block_local_t(qStarts[start_pos + i] - this->qStart); // converting to local coordinate
        
        this->tStarts = new block_local_t[blockCount];
        for (i = 0; i < blockCount; ++i)
            this->tStarts[i] = ulong2block_local_t(tStarts[start_pos + i] - this->tStart); // converting to local coordinate
}

AlignmentRecord::AlignmentRecord(const AlignmentRecord &other)
        : strand(other.strand), qStart(other.qStart), qEnd(other.qEnd),
          tStart(other.tStart), tEnd(other.tEnd),
          blockCount(other.blockCount), sym(other.sym) {
        unsigned long membytes = sizeof(block_local_t) * other.blockCount;
        
        this->blockSizes = new block_local_t[other.blockCount];
        memcpy(this->blockSizes,  other.blockSizes, membytes);
        
        this->qStarts = new block_local_t[other.blockCount];
        memcpy(this->qStarts, other.qStarts, membytes);
        
        this->tStarts = new block_local_t[other.blockCount];
        memcpy(this->tStarts, other.tStarts, membytes);
}

AlignmentRecord::~AlignmentRecord() {
        delete[] blockSizes;
        delete[] qStarts;
        delete[] tStarts;
}

void AlignmentRecord::printRecord() const {
	std::cout << "Strand: " << strand << "\n";
	std::cout << "qStart: " << qStart << "\n";
	std::cout << "qEnd: " << qEnd << "\n";
	std::cout << "tStart: " << tStart << "\n";
	std::cout << "tEnd: " << tEnd << "\n";
	std::cout << "blockCount: " << blockCount << "\n";
	std::cout << "blockSizes: ";
	for (block_local_t i = 0; i < blockCount; i++) std::cout << blockSizes[i] << ",";
	std::cout << std::endl;
	std::cout << "qStarts: ";
	for (block_local_t i = 0; i < blockCount; i++) std::cout << qStarts[i] << ",";
	std::cout << "\n";
	std::cout << "tStarts: ";
	for (block_local_t i = 0; i < blockCount; i++) std::cout << tStarts[i] << ",";
	std::cout << "\n";
	if (sym != nullptr)
		std::cout << "sym hast tStart " << sym->tStart << " and tEnd " << sym->tEnd <<
		". This ones tStart according to sym is " << sym->sym->tStart << "\n";
}

AlignmentRecord *AlignmentRecord::revert() const {
	// all it really does is swap query and target
	std::vector<unsigned int> newBlockSizes;
        std::vector<unsigned long> newQStarts, newTStarts;
        newBlockSizes.reserve(blockCount);
        newQStarts.reserve(blockCount);
        newTStarts.reserve(blockCount);
	if (strand == '+') {
		for (unsigned int i = 0; i < blockCount; i++) {
			newQStarts.push_back(get_tStarts(i)); // must transform to global coordinates to pass to the constructor
			newTStarts.push_back(get_qStarts(i)); // must transform to global coordinates to pass to the constructor
			newBlockSizes.push_back(blockSizes[i]);
		}
	}
	else { // reverse strand - revert order, and make endpoints startpoints
		for (long i = blockCount - 1; i >= 0; i--) {
			newQStarts.push_back(get_tStarts(i) + blockSizes[i]); // must transform to global coordinates to pass to the constructor
			newTStarts.push_back(get_qStarts(i) - blockSizes[i]); // must transform to global coordinates to pass to the constructor
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