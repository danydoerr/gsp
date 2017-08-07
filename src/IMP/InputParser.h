#pragma once
#include <string>
#include <map>
#include <vector>
#include <memory>
#include "AlignmentRecord.h"

/* Parses command line arguments or prints help if none are given. */
void parseCmdArgs(int argc, char** &argv, char* &pslPath,
	unsigned int &minLengh, unsigned int &maxGap, unsigned int &minAlnLength,
	float &minAlnIdentity);

/* Reads a psl file. 
Each line is parsed to an AlignmentRecord. Pointers to all records are stored in result.
Result is sorted by the alignment's starting position in the target sequence. */
void parsePsl(const char* pslPath, std::map<std::string, unsigned long>& speciesStart,
	unsigned int maxGapLength, unsigned int minAlnLength, float minAlnIdentity,
	std::vector<std::shared_ptr<AlignmentRecord>>& result);

/* Organizes AlignmentRecords into buckets with regards to their target positions.
A bucket represents a number of sequence positions, said number being equal to bucketSize.
This makes finding alignments covering a certain positions much faster. */
void fillBuckets(std::vector<std::shared_ptr<AlignmentRecord>>& alns, unsigned int bucketSize,
	std::vector<std::vector<std::shared_ptr<AlignmentRecord>>>& result);