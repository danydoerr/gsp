#pragma once
#include <string>
#include <map>
#include <vector>
#include <deque>
#include <memory>
#include "AlignmentRecord.h"

class InputParser {
    
public:
    /* Constructor */
    InputParser();
    
    /* Parses command line arguments or prints help if none are given. */
    void parseCmdArgs(int argc, char** &argv);
    
    /* Places in variables command line arguments parsed */
    void getCmdLineArgs(std::vector<std::string> &pslPaths,
            unsigned int &minLength, unsigned int &maxGap, unsigned int &minAlnLength,
            float &minAlnIdentity, unsigned int &bucketSize, unsigned int &numThreads);
    
    /* Places in variables command line arguments parsed, except for pslPaths */
    void getCmdLineArgs(unsigned int &minLength, unsigned int &maxGap,
            unsigned int &minAlnLength, float &minAlnIdentity,
            unsigned int &bucketSize, unsigned int &numThreads);

    /* Reads a psl file. 
    Each line is parsed to an AlignmentRecord. Pointers to all records are stored in result.
    Result is sorted by the alignment's starting position in the target sequence. */
    void parsePsl(std::map<std::string, unsigned long>& speciesStart,
            std::deque<AlignmentRecord *>& result);
    
    /* Reads psl files and finds the maximum block size and maximum block start
     * (using local coordinates) inside an alignment */
    void getMaxBlockSizeAndLocalStart(unsigned long &max_bsize, unsigned long &max_start);

private:
    std::vector<std::string> pslPaths;
    unsigned int minLength;
    unsigned int maxGapLength;
    unsigned int minAlnLength;
    float minAlnIdentity;
    unsigned int bucketSize;
    unsigned int numThreads;
    bool printZeroLines;
    bool inputNotPsl;
    
    // We suppose psl lines won't be longer than that
    static const unsigned int MAX_LINE = 32768;
    
    // Used during parse
    char *line; // current line
    unsigned int pos; // position in current line
    unsigned long line_num; // current line number
    std::vector<unsigned long> zeroBlockLines; // lines containing blocks of size 0
    

    /* Parses a single psl line to alignment records (original and reverse,
     * sometimes split) and add them to records vector, returns the number of
     * records added */
    unsigned long recordsFromPsl(std::deque<AlignmentRecord *>& records,             
            std::map<std::string, unsigned long>& speciesStart);
    
    /* Reads a string field  */
    inline std::string getStringField();

    /* Reads and returns a long field value (we assume no sign, just digits) */
    inline unsigned long getLongField();

    /* Reads and returns an int field value (we assume no sign, just digits) */
    inline unsigned int getIntField();

    /* Reads and returns a long subfield value (we assume no sign, just digits, ends with comma) */
    inline unsigned long getLongSubField();

    /* Reads and returns an int subfield value (we assume no sign, just digits, ends with comma) */
    inline unsigned int getIntSubField();

    /* Reads and returns an integer vector from a field composed by a set of int subfields separated and ending by comma + \t */
    inline std::vector<unsigned int> getIntArrayField(unsigned int numberOfSubfields);

    /* Reads and returns an integer vector from a field composed by a set of long subfields separated and ending by comma + \t */
    inline std::vector<unsigned long> getLongArrayField(unsigned int numberOfSubfields);

    /* Advances in line skipping a number of fields */
    inline void skipFields(unsigned int numberOfFields);

    /* Check if sequences in current line were already read. If not, add with its related offset */
    inline void updateSpeciesStart(std::map<std::string, unsigned long>& speciesStart,
            std::string name, unsigned long size);

    /* Adds record and reverse to vector and setup sym pointers */
    inline void setupSymAndAdd(std::deque<AlignmentRecord *>& records, AlignmentRecord *rec);
    
    /* Removes blocks of size 0 and updates related data */
    inline void removeZeroBlocks(unsigned int &blockCount, std::vector<unsigned int> &blockSizes,
            std::vector<unsigned long> &qStarts, std::vector<unsigned long> &tStarts);
    
    /* Prints to stderr message about lines that have removed blocks of size 0 */
    void printZeroBlockInfo(void);
    
    /* For each file in pslPaths read that file and get the psl paths inside it,
     * then replace the strings in pslPaths with the actual psl paths */
    void readPslPaths(void);
};


/* InputParser inline methods */

inline std::string InputParser::getStringField() {
        std::string str;
        str.reserve(16); // should be enough in most cases
        while (line[pos] != '\t')
            str.push_back(line[pos++]);
        ++pos; // move to after \t
        return str;
}

inline unsigned long InputParser::getLongField() {
        unsigned long v = 0;
        while (line[pos] != '\t') {
            v *= 10;
            v += line[pos++] - '0';
        }
        ++pos; // move to after \t
        return v;
}

inline unsigned int InputParser::getIntField() {
        unsigned int v = 0;
        while (line[pos] != '\t') {
            v *= 10;
            v += line[pos++] - '0';
        }
        ++pos; // move to after \t
        return v;
}

inline unsigned long InputParser::getLongSubField() {
// we could join this function with getNumericField function, but an extra || comparison would make it slower
        unsigned long v = 0;
        while (line[pos] != ',') {
            v *= 10;
            v += line[pos++] - '0';
        }
        ++pos; // move to after ,
        return v;
}

inline unsigned int InputParser::getIntSubField() {
// we could join this function with getNumericField function, but an extra || comparison would make it slower
        unsigned int v = 0;
        while (line[pos] != ',') {
            v *= 10;
            v += line[pos++] - '0';
        }
        ++pos; // move to after ,
        return v;
}

inline std::vector<unsigned int> InputParser::getIntArrayField(unsigned int numberOfSubfields) {
        std::vector<unsigned int> values;
        values.reserve(numberOfSubfields);
        for (unsigned int i = 0; i < numberOfSubfields; i++)
            values.push_back(getIntSubField());
        ++pos; // move to after \t (or \n if this is the last field)
        return values;
}

inline std::vector<unsigned long> InputParser::getLongArrayField(unsigned int numberOfSubfields) {
        std::vector<unsigned long> values;
        values.reserve(numberOfSubfields);
        for (unsigned int i = 0; i < numberOfSubfields; i++)
            values.push_back(getLongSubField());
        ++pos; // move to after \t (or \n if this is the last field)
        return values;
}

inline void InputParser::skipFields(unsigned int numberOfFields) {
        for (unsigned int skipped = 0; skipped < numberOfFields; ++pos)
            if (line[pos] == '\t')
                ++skipped;
}

inline void InputParser::updateSpeciesStart(std::map<std::string, unsigned long>& speciesStart,
        std::string name, unsigned long size) {
        if (!speciesStart.count(name)) {
            auto last = speciesStart.find("$");
            auto curLen = last->second;
            speciesStart.insert(std::pair<std::string, unsigned long>(name, curLen));
            last->second = curLen + size;
        }
}

inline void InputParser::setupSymAndAdd(std::deque<AlignmentRecord *>& records, AlignmentRecord *rec) {
        AlignmentRecord *rev = rec->revert();
        rec->sym = rev;
        rev->sym = rec;
        records.push_back(rec);
        records.push_back(rev);
}

inline void InputParser::removeZeroBlocks(unsigned int &blockCount, std::vector<unsigned int> &blockSizes,
        std::vector<unsigned long> &qStarts, std::vector<unsigned long> &tStarts) {
    unsigned int i;
    for (i = 0; i < blockCount && blockSizes[i] != 0; ++i)
        ;
    
    if (i == blockCount) // no zero blocks
        return; // this is the most usual case and we want to know fast
    
    ++i; // i is in the 1st occurrence of 0, we go to the next position
    unsigned int zeros = 1;
    for ( ; i < blockCount; ++i)
        if (blockSizes[i] == 0)
            ++zeros;
        else {
            blockSizes[i - zeros] = blockSizes[i];
            qStarts[i - zeros] = qStarts[i];
            tStarts[i - zeros] = tStarts[i];
        }
    
    blockCount -= zeros;
    blockSizes.resize(blockCount);
    qStarts.resize(blockCount);
    tStarts.resize(blockCount);
    
    zeroBlockLines.push_back(line_num);
}
