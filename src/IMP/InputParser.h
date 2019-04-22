#pragma once
#include <string>
#include <map>
#include <vector>
#include <memory>
#include "AlignmentRecord.h"

/* Parses command line arguments or prints help if none are given. */
void parseCmdArgs(int argc, char** &argv, std::vector<char*> &pslPath,
	unsigned int &minLengh, unsigned int &maxGap, unsigned int &minAlnLength,
	float &minAlnIdentity, unsigned int &bucketSize, unsigned int &numThreads);

/* Reads a psl file. 
Each line is parsed to an AlignmentRecord. Pointers to all records are stored in result.
Result is sorted by the alignment's starting position in the target sequence. */
void parsePsl(std::vector<char*> &pslPath, std::map<std::string, unsigned long>& speciesStart,
	unsigned int maxGapLength, unsigned int minAlnLength, float minAlnIdentity,
	std::vector<AlignmentRecord *>& result);

/* Parses a single psl line to alignment records (original and reverse,
 * sometimes split) and add them to records vector, returns the number of
 * records added */
unsigned int recordsFromPsl(std::vector<AlignmentRecord *>& records, const char *line,
        unsigned int maxGapLength, unsigned int minAlnLength, float minAlnIdentity,
        std::map<std::string, unsigned long>& speciesStart);