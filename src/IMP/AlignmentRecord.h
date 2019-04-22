#pragma once
#include <vector>
#include <set>
#include <memory>
#include <string>

/* ADJUSTABLE MEMORY OPTIMIZATION */
/* HERE WE CAN SET THE TYPE USED FOR STORING BLOCKS SIZES AND STARTS AS LOCAL COORDINATES */
/* NO NOT CHANGE THE NEXT DEFINES */
#define BLOCKS_ULONG 0  // unsigned long (no optimization)
#define BLOCKS_UINT 1   // unsigned int (should fit blocks and their corresponding data)
#define BLOCKS_USHORT 2 // unsigned short (need to be careful)
/* SET BLOCKS_SIZE TO ONE OF ABOVE TO DEFINE THE VARIABLE SIZE USED FOR BLOCKS */
#define BLOCKS_SIZE BLOCKS_USHORT // <--- set here the variable size
/* NO NOT CHANGE THE NEXT DEFINES */
#if BLOCKS_SIZE == BLOCKS_USHORT
#define block_local_t unsigned short
#elif BLOCKS_SIZE == BLOCKS_UINT
#define block_local_t unsigned int
#else
#define block_local_t unsigned long
#endif


/* Representation of all needed information of a single psl line.
Additionally, contains a pointer to sym, the AlignmentRecord of its inverse alignment. */
class AlignmentRecord {
public:
	char strand; // + (forward) or - (reverse)
	unsigned long qStart; // alignment start position in query
	unsigned long qEnd; // alignment end position in query
	unsigned long tStart; // alignment start position in target
	unsigned long tEnd; // alignment end position in target
	block_local_t blockCount; // number of blocks in aln
        block_local_t *blockSizes; // size of each block

private:
        // store local coordinates, public accessible by global_qStarts/global_tStarts methods
        // unlike the psl file, when the strand is "-", the starts are relative to the beginning instead of to the end of sequence
	block_local_t *qStarts; // start position of each block in query
	block_local_t *tStarts; // start position of each block in target
        
        /* Converts unsigned long to unsigned int, throwing an exception if doesn't fit 
         * (even that this adds some overhead, we have to do this to prevent
         * malfunctioning since we use smaller variables to try to save some
         * memory, and we cannot allow the program to continue if some value
         * can't fit these variables
         */
        inline unsigned int ulong2uint(const long &ul) const;
        
        /* Converts unsigned long to unsigned short, throwing an exception if doesn't fit
         * (same as above)
         */
        inline unsigned short ulong2ushort(const long &ul) const;
        
public:
	AlignmentRecord *sym; // pointer to inverse alignment

	/* Constructor (qStarts and tStarts are global coordinates) */
	AlignmentRecord(char strand,
		unsigned long qStart, unsigned long qEnd,
		unsigned long tStart, unsigned long tEnd,
		unsigned int blockCount, std::vector<unsigned int> blockSizes,
		std::vector<unsigned long> qStarts, std::vector<unsigned long> tStarts);
        
        /* Same as before, but considers only a subinterval of the alignment,
         * consisting of "blockCount" blocks >= 1 starting from start_pos >= 0*/
        AlignmentRecord(char strand,
		unsigned long qStart, unsigned long qEnd,
		unsigned long tStart, unsigned long tEnd,
		unsigned int blockCount, std::vector<unsigned int> blockSizes,
		std::vector<unsigned long> qStarts, std::vector<unsigned long> tStarts,
                unsigned int start_pos);
        
        /* Copy constructor */
        AlignmentRecord(const AlignmentRecord &other);
        
        /* Destructor */
        ~AlignmentRecord();
        
	bool operator < (const AlignmentRecord other) { return tEnd < other.tEnd; }; // for sorting

	/* Prints all attributes of an AlignmentRecord to STDOUT. */
	void printRecord() const;
	unsigned long getLength() const { return tEnd - tStart; };

	/* Calculates AlignmentRecord of the inverse alignment and returns a pointer to it. */
	AlignmentRecord *revert() const;
        
        /* Returns one index of qStarts in global coordinates. */
        inline unsigned long get_qStarts(unsigned int idx) const { return qStarts[idx] + qStart; };
        
        /* Returns one index of qStarts in global coordinates */
        inline unsigned long get_tStarts(unsigned int idx) const { return tStarts[idx] + tStart; };
        
        /* Iterator over qStarts, tStarts and blockSizes (global coordinates) implementation */
        class iterator : public std::iterator<std::forward_iterator_tag, unsigned long>
        {           
        public:
            enum Type : char { QUERY = 'Q', TARGET = 'T', BLOCK_SIZES = 'B'};
            inline iterator(const AlignmentRecord *record, Type t, unsigned int idx = 0);
            inline iterator(const iterator& i);
            inline iterator& operator=(const iterator& i);
            inline iterator& operator++();
            inline iterator operator++(int);
            inline iterator operator+(const int& rhs);
            inline unsigned long operator*() const;
            inline unsigned long operator->() const;
            inline bool operator==(const iterator& i) const;
            inline bool operator!=(const iterator& i) const;
            
        private:
            const AlignmentRecord *record;
            Type type;
            unsigned int cur_idx;
        };

        inline iterator begin_qStarts() const;
        inline iterator end_qStarts() const;
        inline iterator begin_tStarts() const;
        inline iterator end_tStarts() const;
        inline iterator begin_blockSizes() const;
        inline iterator end_blockSizes() const;
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



/* Alignment Record inline methods */

inline unsigned int AlignmentRecord::ulong2uint(const long &ul) const {
        unsigned int ui = ul;
        if (ui != ul) throw std::range_error("Cannot fit this number in an unsigned int: " + std::to_string(ul) + " (" + __FILE__ + ":" + std::to_string(__LINE__) + ")");
        return ui;   
}

inline unsigned short AlignmentRecord::ulong2ushort(const long &ul) const {
        unsigned short uh = ul;
        if (uh != ul) throw std::range_error("Cannot fit this number in an unsigned short: " + std::to_string(ul) + " (" + __FILE__ + ":" + std::to_string(__LINE__) + ")");
        return uh; 
}

inline AlignmentRecord::iterator AlignmentRecord::begin_qStarts() const
{
  return iterator(this, iterator::Type::QUERY);
}

inline AlignmentRecord::iterator AlignmentRecord::end_qStarts() const
{
  return iterator(this, iterator::Type::QUERY, blockCount);
}

inline AlignmentRecord::iterator AlignmentRecord::begin_tStarts() const
{
  return iterator(this, iterator::Type::TARGET);
}

inline AlignmentRecord::iterator AlignmentRecord::end_tStarts() const
{
  return iterator(this, iterator::Type::TARGET, blockCount);
}

inline AlignmentRecord::iterator AlignmentRecord::begin_blockSizes() const
{
  return iterator(this, iterator::Type::BLOCK_SIZES);
}

inline AlignmentRecord::iterator AlignmentRecord::end_blockSizes() const
{
  return iterator(this, iterator::Type::BLOCK_SIZES, blockCount);
}

inline AlignmentRecord::iterator::iterator(const AlignmentRecord *record, Type type, unsigned int idx) :
    record(record),
    type(type),
    cur_idx(idx)
{}

inline AlignmentRecord::iterator::iterator(const iterator& i) :
    record(i.record),
    type(i.type),
    cur_idx(i.cur_idx)
{}

inline AlignmentRecord::iterator& AlignmentRecord::iterator::operator=(const iterator& i)
{ 
    record = i.record;
    type = i.type;
    cur_idx = i.cur_idx;
    return *this; 
}

inline AlignmentRecord::iterator& AlignmentRecord::iterator::operator++()
{
    ++cur_idx;
    return *this; 
}

inline AlignmentRecord::iterator AlignmentRecord::iterator::operator++(int)
{ 
    iterator tmp(*this);
    ++cur_idx;
    return tmp; 
}

inline AlignmentRecord::iterator AlignmentRecord::iterator::operator+(const int& rhs)
{
    return iterator(record, type, cur_idx + rhs);
}

inline unsigned long AlignmentRecord::iterator::operator*() const
{
    if (type == QUERY)
        return record->get_qStarts(cur_idx);
    else if (type == TARGET)
        return record->get_tStarts(cur_idx);
    else
        return record->blockSizes[cur_idx];
}

inline unsigned long AlignmentRecord::iterator::operator->() const
{
    if (type == QUERY)
        return record->get_qStarts(cur_idx);
    else if (type == TARGET)
        return record->get_tStarts(cur_idx);
    else
        return record->blockSizes[cur_idx];
}

inline bool AlignmentRecord::iterator::operator==(const iterator& i) const
{
  return record == i.record && type == i.type && cur_idx == i.cur_idx; 
}

inline bool AlignmentRecord::iterator::operator!=(const iterator& i) const
{
  return record != i.record || type != i.type || cur_idx != i.cur_idx;
}
