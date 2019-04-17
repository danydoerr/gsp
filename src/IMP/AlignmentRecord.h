#pragma once
#include <vector>
#include <set>
#include <memory>

// Type used for storing block sizes and starts in alignment records
#define align_rec_block_t unsigned long

/* Represantation of all needed information of a single psl line.
Additionally, contains a pointer to sym, the AlignmentRecord of its inverse alignment. */
struct AlignmentRecord {
	char strand; // + (forward) or - (reverse)
	unsigned long qStart; // alignment start position in query
	unsigned long qEnd; // alignment end position in query
	unsigned long tStart; // alignment start position in target
	unsigned long tEnd; // alignment end position in target
	//unsigned int blockCount; // number of blocks in aln (not needed, we can just check blockSizes.size()
	std::vector<unsigned int> blockSizes; // size of each block
	std::vector<align_rec_block_t> qStarts; // start position of each block in query
	std::vector<align_rec_block_t> tStarts; // start position of each block in target
	AlignmentRecord *sym; // pointer to inverse alignment

	/* Constructor. */
	AlignmentRecord(char strand,
		unsigned long qStart, unsigned long qEnd,
		unsigned long tStart, unsigned long tEnd,
		std::vector<unsigned int> blockSizes,
		std::vector<align_rec_block_t> qStarts, std::vector<align_rec_block_t> tStarts);
	bool operator < (const AlignmentRecord other) { return tEnd < other.tEnd; }; // for sorting

	/* Prints all attributes of an AlignmentRecord to STDOUT. */
	void printRecord() const;
	unsigned long getLength() const { return tEnd - tStart; };

	/* Calculates AlignmentRecord of the inverse alignment and returns a pointer to it. */
	AlignmentRecord *revert() const;
        
        /* Returns the number of blocks (old blockCount member variable) */
        inline unsigned int blockCount() const { return blockSizes.size(); }; // this would be implicitly inline in a classes, not sure in structs
};

struct Breakpoint {
	unsigned long position;

	Breakpoint(unsigned long position);
	bool operator < (const Breakpoint other) { return position < other.position; }; // for sorting
	bool operator == (const Breakpoint other) { return position == other.position; };
};

/* A Regions in a sequence, defined by two position (start and end). */
struct Region {
	unsigned long first;
	unsigned long last;

	Region(unsigned long first, unsigned long last);
	unsigned long getLength() const { return last - first + 1; };
	unsigned long getMiddlePos() const { return (first + last) / 2; }
	bool operator == (const Region other) { return last == other.last && first == other.first; };

	/* Region with last position further to the right is greater.
	If last of both is equal, region with first position further to the right is greater */
	bool operator < (const Region other) {
		if (last == other.last) return first > other.first;
		else return last < other.last;
	}
};

/* A waste region. Basically qual to region, but sorted differently. */
struct WasteRegion : public Region {
	WasteRegion(unsigned long pos);
	WasteRegion(Region atom);

	bool operator < (const WasteRegion other) {
		if (first == other.first) return last < other.last;
		else return first < other.first;
	}
};

struct dpPosition {
	unsigned int idx;
	double cost;
	bool dist;
	unsigned long prev;
	std::vector<unsigned int> coveringIds;
	std::vector<unsigned int> notCoveringIds;
	dpPosition(unsigned int idx);
};

struct dpStats {
	double cost;
	bool dist;
	unsigned long prev;
	dpStats(double cost, bool dist, unsigned long prev);
};