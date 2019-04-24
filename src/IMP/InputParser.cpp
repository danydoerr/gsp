#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm>

#include "InputParser.h"
#include "AlignmentRecord.h"
#include "Util.h"

InputParser::InputParser() {
    // Default values
    maxGapLength = 13;
    minAlnLength = 13;
    minLength = 250;
    bucketSize = 1000;
    numThreads = 1;
    minAlnIdentity = 0.8f;
    printZeroLines = false;
}

void InputParser::parseCmdArgs(int argc, char** &argv) {
	if (argc <= 1) {
		std::cerr << "Usage: atomizer <psl file(s)> [options]\n\n"
                        << "Multiple input psl files may be given, but always as first arguments.\n"
			<< "Optional arguments are given after their descriptor. The descriptor is NOT case-sensitive. \n"
			<< "If an optional argument is not given, the default value will be used.\n"
			<< "--minLength <minLength>: The minimum length an atom must have (defualt: 250).\n"
			<< "--minIdent <minIdent>: Minimum identity an alignment must have to be considered, "
			<< "alignments with lower identity will be skipped (default: 80).\n"
			<< "--maxGap <maxGap>: The maximum gap length inside of an alignment. If this length is "
			<< "exceeded, the alignment will be split in two (default: 13).\n"
			<< "--minAlnLength <minAlnLength>: The minimal length an alignment must have to be considered. "
			<< "Shorter alignments are ignored (default: 13).\n"
			<< "--bucketSize <size>: Size of buckets used to find covering alignments, "
			<< "increase if you run out of memory (default: 1000).\n"
			<< "--numThreads <num>: Number of threads to run IMP algorithm (default: 1).\n"
                        << "--printZeroLines: Print line numbers with blocks of size 0 (default: no)."
			<< std::endl;
		exit(EXIT_SUCCESS);
	}
	int mandatoryArgs = 0;
	for (int i = 1; i < argc; i++) {
		std::string arg = argv[i];
		if (arg.substr(0, 2) != "--") mandatoryArgs++;
		else break;
	}
	if (mandatoryArgs < 1) {
		std::cerr << "Insufficient arguments, call without arguments to get instructions." << std::endl;
		exit(EXIT_FAILURE);
	}
        for (int i = 1; i <= mandatoryArgs; i++)
            pslPaths.push_back(argv[i]);
	for (int i = mandatoryArgs + 1; i < argc; i++) {
		std::string arg = argv[i];
		std::transform(arg.begin(), arg.end(), arg.begin(), tolower);
		try {
			if (arg == "--minlength") minLength = std::stoul(argv[++i]);
			else if (arg == "--minident") minAlnIdentity = std::stoul(argv[++i]) / 100.0f;
			else if (arg == "--maxgap") maxGapLength = std::stoul(argv[++i]);
			else if (arg == "--minalnlength") minAlnLength = std::stoul(argv[++i]);
			else if (arg == "--bucketsize") bucketSize = std::stoul(argv[++i]);
			else if (arg == "--numthreads") numThreads = std::stoul(argv[++i]);
                        else if (arg == "--printzerolines") printZeroLines = true;
			else {
				std::cerr << "Unknown argument " << arg << ". Call without arguments for instructions." << std::endl;
				exit(EXIT_FAILURE);
			}
		}
		catch (std::invalid_argument) {
			std::cerr << "The value for argument " << arg << " could not be parsed." << std::endl;
			exit(EXIT_FAILURE);
		}
	}
}

void InputParser::getCmdLineArgs(std::vector<char*> &pslPaths,
        unsigned int &minLength, unsigned int &maxGapLength, unsigned int &minAlnLength,
        float &minAlnIdentity, unsigned int &bucketSize, unsigned int &numThreads) {
    
    pslPaths = this->pslPaths;
    minLength = this->minLength;
    maxGapLength = this->maxGapLength;
    minAlnLength = this->minAlnLength;
    minAlnIdentity = this->minAlnIdentity;
    bucketSize = this->bucketSize;
    numThreads = this->numThreads;
}

/* Parses a single psl line to alignment records (original and reverse,
 * sometimes split) and add them to records vector, returns the number of
 * records added */
unsigned int InputParser::recordsFromPsl(std::vector<AlignmentRecord *>& records,
        std::map<std::string, unsigned long>& speciesStart) {
    
        unsigned int orig_size = 0; // records size before adding new records
        pos = 0; // position in line
        
        { // skip low quality alignments
            unsigned int matches = getIntField();
            unsigned int mismatches = getIntField();
            unsigned int repmatches = getIntField();
            matches += repmatches;
            if (matches == 0) return 0;
            unsigned int all = matches + mismatches;
            if (static_cast<float>(matches) / static_cast<float>(all) < minAlnIdentity) return 0;    
            // Comments from original parser:
            /* removed these filters for now - filter input psl by hand instead when needed
             * if (curRec.tStart > curRec.qStart) continue; // only one version of symmetric alignments 
             * if (curRec.tStart == curRec.qStart && curRec.tEnd == curRec.qEnd)
             *	continue; // skip alignments that align a region to itself*/
        }
        
        skipFields(5);
        
        // fields variables, in the order they appear
        const char strand = line[pos++];
        ++pos; // we should be at \t now, move past it
        
        const std::string qName = getStringField();
        const unsigned long qSize = getLongField();
        updateSpeciesStart(speciesStart, qName, qSize); // check if sequence is in the map, if not, add it
        const unsigned long qOffset = speciesStart.find(qName)->second; // offset positions for concatenated sequence
        const unsigned long qStart = getLongField() + qOffset;
        const unsigned long qEnd = getLongField() + qOffset;
        
        const std::string tName = getStringField();
        const unsigned long tSize = getLongField();
        updateSpeciesStart(speciesStart, tName, tSize);
        const unsigned long tOffset = speciesStart.find(tName)->second;
        const unsigned long tStart = getLongField() + tOffset;
        const unsigned long tEnd = getLongField() + tOffset;
        
        unsigned int blockCount = getIntField();
        
        std::vector<unsigned int> blockSizes = getIntArrayField(blockCount);
        std::vector<unsigned long> qStarts = getLongArrayField(blockCount);
        std::vector<unsigned long> tStarts = getLongArrayField(blockCount);
        
        removeZeroBlocks(blockCount, blockSizes, qStarts, tStarts);
        
	// shift start positions by offset, if on reverse strand (only query) recompute w.r.t. starts to start of sequence
	if (strand == '+')
            for (auto i = qStarts.begin(); i != qStarts.end(); i++)
                *i += qOffset;
	else
            for (auto i = qStarts.begin(); i != qStarts.end(); i++)
                *i = qSize - *i + qOffset;
	for (auto i = tStarts.begin(); i != tStarts.end(); i++)
            *i += tOffset;
        
        // Splits alignment in parts if it contains gaps longer than maxGapLength,
        // adds to results only if split parts are longer than minAlnLength
	unsigned int start = 0, length, end;
	for (end = 0; end < blockCount - 1; end++)	
            if (tStarts[end + 1] - (tStarts[end] + blockSizes[end]) > maxGapLength
                    || (strand == '+' && qStarts[end + 1] - (qStarts[end] + blockSizes[end]) > maxGapLength)
                    || (strand == '-' && qStarts[end] - (qStarts[end + 1] + blockSizes[end]) > maxGapLength)) {
                length = (tStarts[end] + blockSizes[end]) - tStarts[start];
                if (length > minAlnLength)
                    setupSymAndAdd(records, new AlignmentRecord(strand, qStart, qEnd, tStart, tEnd, end-start+1, blockSizes, qStarts, tStarts, start));
                start = end + 1;
            }
	length = (tStarts[end] + blockSizes[end]) - tStarts[start];
	if (length > minAlnLength)
            setupSymAndAdd(records, new AlignmentRecord(strand, qStart, qEnd, tStart, tEnd, end-start+1, blockSizes, qStarts, tStarts, start));
        
        return records.size() - orig_size;
}

void InputParser::parsePsl(std::map<std::string, unsigned long>& speciesStart,
	std::vector<AlignmentRecord *>& result) {
    
        line = new char[MAX_LINE]; // I'm not sure if it is a good idea to allocate this big block in the stack
        zeroBlockLines.reserve(1024);
                
	std::ifstream pslFile;
        for (auto psl : pslPaths) {
            std::cerr << "Reading " << psl << "... ";
            pslFile.open(psl);
            if (pslFile.is_open()) {
                    line_num = 0;
                    while (!pslFile.getline(line, MAX_LINE).eof()) {
                            ++line_num;
                            if (line[0] == '#' || line[0] == '\0') continue; // skip comments and empty lines
                            
                            recordsFromPsl(result, speciesStart);
                    }
                    pslFile.close();
                    std::cerr << "Done." << std::endl;
                    printZeroBlockInfo();
                    zeroBlockLines.clear();
            }
            else {
                    std::cerr << "ERROR: psl file could not be opened!" << std::endl;
                    exit(EXIT_FAILURE);
            }
        }
        delete[] line;
}

void InputParser::getMaxBlockSizeAndLocalStart(unsigned long &max_bsize, unsigned long &max_start) {
        max_bsize = 0;
        max_start = 0;
        line = new char[MAX_LINE]; // I'm not sure if it is a good idea to allocate this big block in the stack
        
        std::map<std::string, unsigned long> speciesStarts; // maps species name to their starting position in concatenated string
        speciesStarts = { {"$", 0} };
        
        std::vector<AlignmentRecord *> records;
        records.reserve(128);
        
        std::vector<std::string> dontFit;
        dontFit.reserve(1024);
        unsigned long last_size = 0;
        
	std::ifstream pslFile;
        for (auto psl : pslPaths) {
            std::cout << "Reading " << psl << "... ";
            pslFile.open(psl);
            if (pslFile.is_open()) {
                    line_num = 0;
                    while (!pslFile.getline(line, MAX_LINE).eof()) {
                            ++line_num;
                            if (line[0] == '#' || line[0] == '\0') continue; // skip comments and empty lines
                            
                            try {
                                recordsFromPsl(records, speciesStarts);
                            } catch (std::range_error& e) { 
                                dontFit.push_back(std::string("Input line ") + std::to_string(line_num) + ": " + std::string(e.what()));
                            }
                            
                            for (auto curRec : records) {
                                for (auto size = curRec->begin_blockSizes(); size != curRec->end_blockSizes(); size++)
                                    if (*size > max_bsize)
                                        max_bsize = *size;
                                for (auto start = curRec->begin_qStarts(); start != curRec->end_qStarts(); start++)
                                    if (*start - curRec->qStart > max_bsize)
                                        max_start = *start - curRec->qStart;
                                for (auto start = curRec->begin_tStarts(); start != curRec->end_tStarts(); start++)
                                    if (*start - curRec->tStart > max_bsize)
                                        max_start = *start - curRec->tStart;
                                delete curRec;
                            }
                            
                            records.clear();
                    }
                    pslFile.close();
                    std::cout << "Done." << std::endl;
                    unsigned long newEntries = dontFit.size() - last_size;
                    if (newEntries > 0) {
                        std::cout << "\t" << newEntries << " lines have values that don't fit in <" << sizeof(block_local_t) << " bytes>" << ", the first ones for each line are:" << std::endl;
                        for (auto msg = dontFit.end() - newEntries; msg != dontFit.end(); ++msg)
                            std::cout << "\t" << *msg << std::endl;
                    }
                    last_size = dontFit.size();
            }
            else {
                    std::cerr << "ERROR: psl file \"" << psl <<"\" could not be opened!" << std::endl;
                    exit(EXIT_FAILURE);
            }
        }
        std::cerr << dontFit.size() << " lines with values that don't fit in <" << sizeof(block_local_t) << " bytes>" << std::endl;
        delete[] line;
}

void InputParser::printZeroBlockInfo(void) {
    if (zeroBlockLines.size() == 0)
        return;
    
    std::cerr << "\tBlocks of size 0 found and removed in " << zeroBlockLines.size() << " alignments";
    if (!printZeroLines) {
        std::cerr << std::endl;
        return;
    }
    
    std::cerr << ", lines ";
    bool first = true;
    for (auto l : zeroBlockLines)
        if (first) {
            std::cerr << l;
            first = false;
        }
        else
            std::cerr << ", " << l;
    std::cerr << std::endl;
}