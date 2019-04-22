#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <string>

#include "Util.h"
#include "InputParser.h"

// We suppose lines won't be longer than that
#define MAX_LINE 32768


void analyzePsl(std::vector<char*> &pslPath,
	unsigned int maxGapLength, unsigned int minAlnLength, float minAlnIdentity) {
    
        unsigned long max_bsize = 0, max_start = 0;
        char *line = new char[MAX_LINE]; // I'm not sure if it is a good idea to allocate this big block in the stack
        
        std::map<std::string, unsigned long> speciesStarts; // maps species name to their starting position in concatenated string
        speciesStarts = { {"$", 0} };
        
        std::vector<AlignmentRecord *> records;
        records.reserve(128);
        
	std::ifstream pslFile;
        for (auto psl : pslPath) {
            std::cerr << "Reading " << psl << "... ";
            pslFile.open(psl);
            if (pslFile.is_open()) {
                    while (!pslFile.getline(line, MAX_LINE).eof()) {
                            if (line[0] == '#') continue; // skip comments
                            
                            recordsFromPsl(records, line, maxGapLength, minAlnLength, minAlnIdentity, speciesStarts);
                            
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
                            }
                            
                            records.clear();
                    }
                    pslFile.close();
                    std::cerr << "Done." << std::endl;
            }
            else {
                    std::cerr << "ERROR: psl file could not be opened!" << std::endl;
                    exit(EXIT_FAILURE);
            }
        }
        std::cerr << "MAX SIZE: " << max_bsize << ", MAX START: " << max_start << std::endl;
        delete[] line;
}


int main(int argc, char *argv[]) {
    std::vector<char *> infiles;
    unsigned int maxGapLength = 13, minAlnLength = 13, minLength = 250, bucketSize = 1000, numThreads = 1;
    float minAlnIdentity = 0.8f;
    
    parseCmdArgs(argc, argv, infiles, minLength, maxGapLength, minAlnLength, minAlnIdentity, bucketSize, numThreads);
    
    analyzePsl(infiles, maxGapLength, minAlnLength, minAlnIdentity);
    
    return 0;
}