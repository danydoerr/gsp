#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm>

#include "InputParser.h"
#include "AlignmentRecord.h"
#include "Util.h"

void parseCmdArgs(int argc, char** &argv, std::vector<char*> &pslPath,
	unsigned int &minLengh, unsigned int &maxGap, unsigned int &minAlnLength,
	float &minAlnIdentity, unsigned int &bucketSize, unsigned int &numThreads) {
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
			<< "--bucketSize: Size of buckets used to find covering alignments, "
			<< "increase if you run out of memory (default: 1000).\n"
			<< "--numThreads: Number of threads to run IMP algorithm (default: 1)."
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
            pslPath.push_back(argv[i]);
	for (int i = mandatoryArgs + 1; i < argc; i++) {
		std::string arg = argv[i];
		std::transform(arg.begin(), arg.end(), arg.begin(), tolower);
		try {
			if (arg == "--minlength") minLengh = std::stoul(argv[++i]);
			else if (arg == "--minident") minAlnIdentity = std::stoul(argv[++i]) / 100.0f;
			else if (arg == "--maxgap") maxGap = std::stoul(argv[++i]);
			else if (arg == "--minalnlength") minAlnLength = std::stoul(argv[++i]);
			else if (arg == "--bucketsize") bucketSize = std::stoul(argv[++i]);
			else if (arg == "--numthreads") numThreads = std::stoul(argv[++i]);
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

/* Splits input string at each occurence of char delim. */
const std::vector<std::string> splitString(const std::string& input, const char delim,
	unsigned int expectedSize = 1) {
	std::string buffer{ "" };
	std::vector<std::string> result;
	result.reserve(expectedSize);
	for (auto c : input) {
		if (c != delim) buffer+=c;
		else {
			result.push_back(buffer);
			buffer.clear();
		}
	}
	if (!buffer.empty()) result.push_back(buffer); // don't forget last part of string
	return result;
}

/* Turns input vector of strings into a vector of numbers.
Will end the program if one of the strings does not represent an integer number. */
const std::vector<unsigned long> stringVecToLongVec(const std::vector<std::string>& input) {
	std::vector<unsigned long> result;
        result.reserve(input.size());
	for (auto i : input) {
		try {
			result.push_back(std::stoul(i));
		}
		catch (std::invalid_argument) {
			std::cerr << "ERROR: Input parser failed to parse string to long, input was " << i << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	return result;
}

/* Turns input vector of strings into a vector of numbers.
Will end the program if one of the strings does not represent an integer number. */
const std::vector<unsigned int> stringVecToIntVec(const std::vector<std::string>& input) {
	std::vector<unsigned int> result;
        result.reserve(input.size());
	for (auto i : input) {
		try {
			result.push_back(stoui(i));
		}
		catch (std::invalid_argument) {
			std::cerr << "ERROR: Input parser failed to parse string to unsigned int, input was " << i << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	return result;
}

/* Turns input vector of strings into a vector of numbers.
Will end the program if one of the strings does not represent an integer number. */
const std::vector<unsigned short> stringVecToShortVec(const std::vector<std::string>& input) {
	std::vector<unsigned short> result;
        result.reserve(input.size());
	for (auto i : input) {
		try {
			result.push_back(stouh(i));
		}
		catch (std::invalid_argument) {
			std::cerr << "ERROR: Input parser failed to parse string to unsigned short, input was " << i << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	return result;
}

/* Parses a single psl line to an AlignmentRecord. */
const AlignmentRecord recordFromPsl(const std::vector<std::string>& pslLine,
	std::map<std::string, unsigned long>& speciesStart) {
	const std::string qName = pslLine[9], tName = pslLine[13];
	const char strand = pslLine[8].front();
	const unsigned long qSize = std::stoul(pslLine[10]), tSize = std::stoul(pslLine[14]);
	// check if sequences in current line were already read. If not, add them.
	if (!speciesStart.count(qName)) {
		auto last = speciesStart.find("$");
		auto curLen = last->second;
		speciesStart.insert(std::pair<std::string, unsigned long>(qName, curLen));
		last->second = curLen + qSize;
	}
	if (!speciesStart.count(tName)) {
		auto last = speciesStart.find("$");
		auto curLen = last->second;
		speciesStart.insert(std::pair<std::string, unsigned long>(tName, curLen));
		last->second = curLen + tSize;
	}

	// offset positions for concatenated sequence
	const unsigned long qOffset = speciesStart.find(qName)->second, tOffset = speciesStart.find(tName)->second;
	unsigned long qStart = std::stoul(pslLine[11]) + qOffset, qEnd = std::stoul(pslLine[12]) + qOffset,
		tStart = std::stoul(pslLine[15]) + tOffset, tEnd = std::stoul(pslLine[16]) + tOffset;
        unsigned int blockCount = stoui(pslLine[17]);
	std::vector<unsigned long> qStarts = stringVecToLongVec(splitString(pslLine[19], ',', blockCount)),
		tStarts = stringVecToLongVec(splitString(pslLine[20], ',', blockCount));
	// shift start positions by offest. If on reverse strand, recompute query w.r.t. starts to start of sequence
	if (strand == '+')
		for (auto i = qStarts.begin(); i != qStarts.end(); i++)
			*i += qOffset;
	else {
		//auto qSize = std::stoul(pslLine[10]); // already computed
		for (auto i = qStarts.begin(); i != qStarts.end(); i++)
			*i = qSize - *i + qOffset;
	}
	for (auto i = tStarts.begin(); i != tStarts.end(); i++)
		*i += tOffset;

        std::vector<unsigned int> blockSizes = stringVecToIntVec(splitString(pslLine[18], ',', blockCount));        
	return AlignmentRecord(strand, qStart, qEnd, tStart, tEnd, blockCount,
		blockSizes, qStarts, tStarts);
}

/* Returns a part of the input AlignmentRecord, from startBlock to endBlock, including both. */
const AlignmentRecord cutRecord(const AlignmentRecord &aln, unsigned int startBlock, unsigned int endBlock) {
	if (startBlock == 0 && endBlock == ((unsigned int) aln.blockCount) - 1) return aln;
	unsigned long qStart, qEnd, tStart = aln.get_tStarts(startBlock),
		tEnd = aln.get_tStarts(endBlock) + aln.blockSizes[endBlock];
	if (aln.strand == '+') {
		qStart = aln.get_qStarts(startBlock);
		qEnd = aln.get_qStarts(endBlock) + aln.blockSizes[endBlock];
	} else { // strand == '-'
		qStart = aln.get_qStarts(endBlock) - aln.blockSizes[endBlock];
		qEnd = aln.get_qStarts(startBlock);
	}
	unsigned int blockCount = endBlock - startBlock + 1;
	std::vector<unsigned int> blockSizes(aln.begin_blockSizes() + startBlock, aln.begin_blockSizes() + endBlock + 1);
        std::vector<unsigned long> qStarts(aln.begin_qStarts() + startBlock, aln.begin_qStarts() + endBlock + 1),
		tStarts(aln.begin_tStarts() + startBlock, aln.begin_tStarts() + endBlock + 1);
	return AlignmentRecord(aln.strand, qStart, qEnd, tStart, tEnd, blockCount, blockSizes, qStarts, tStarts);
}

/* Splits input AlignmentRecord in parts if it contains gaps longer than maxGapLength .
Stores splitted parts in result only if they are longer than minAlnLength. */
void splitRecord(const AlignmentRecord &aln,
	unsigned int maxGapLength, unsigned int minAlnLength, std::vector<AlignmentRecord> &result) {
	int start = 0;
	for (size_t i = 0; i < ((unsigned int) aln.blockCount) - 1; i++) {
		// check for long gaps
		if (aln.get_tStarts(i + 1) - (aln.get_tStarts(i) + aln.blockSizes[i]) > maxGapLength
			|| (aln.strand == '+' && aln.get_qStarts(i + 1) - (aln.get_qStarts(i) + aln.blockSizes[i]) > maxGapLength)
			|| (aln.strand == '-' && aln.get_qStarts(i) - (aln.get_qStarts(i + 1) + aln.blockSizes[i]) > maxGapLength)) {
			AlignmentRecord splitRec = cutRecord(aln, start, i);
			start = i + 1;
			if (splitRec.getLength() > minAlnLength)
				result.push_back(splitRec);
		}
	}
	AlignmentRecord splitRec = cutRecord(aln, start, aln.blockCount-1);
	if (splitRec.getLength() > minAlnLength)
		result.push_back(splitRec);
}

void parsePsl(std::vector<char*> &pslPath, std::map<std::string, unsigned long>& speciesStart,
	unsigned int maxGapLength, unsigned int minAlnLength, float minAlnIdentity,
	std::vector<AlignmentRecord *>& result) {           
	std::ifstream pslFile;
        for (auto psl : pslPath) {
            std::cerr << "Reading " << psl << "... ";
            pslFile.open(psl);
            if (pslFile.is_open()) {
                    std::string line;
                    while (std::getline(pslFile, line)) {
                            if (line.front() == '#') continue; // skip comments
                            auto splitLine = splitString(line, '\t', 21);
                            if (splitLine.size() != 21) {
                                    std::cerr << "ERROR: read illegal line while reading psl file. Line was:" << std::endl;
                                    std::cerr << line << std::endl;
                                    exit(EXIT_FAILURE);
                            }
                            try { // skip low quality alignments
                                    int matches = std::stoi(splitLine[0]) + std::stoi(splitLine[2]); // matches + repMatches
                                    if (matches == 0) continue;
                                    int all = matches + std::stoi(splitLine[1]); // + misMatches
                                    if (static_cast<float>(matches) / static_cast<float>(all) < minAlnIdentity) continue;
                            }
                            catch (std::invalid_argument) {
                                    std::cerr << "ERROR: invalid argument when parsing matches while reading psl." << std::endl;
                                    std::cerr << "These three should be numbers: " << splitLine[0] << ", " << splitLine[2]
                                            << ", " << splitLine[1] << std::endl;
                                    exit(EXIT_FAILURE);
                            }
                            AlignmentRecord curRec = recordFromPsl(splitLine, speciesStart);
                            /* removed these filters for now - filter input psl by hand instead when needed
                             * if (curRec.tStart > curRec.qStart) continue; // only one version of symmetric alignments 
                             * if (curRec.tStart == curRec.qStart && curRec.tEnd == curRec.qEnd)
                             *	continue; // skip alignments that align a region to itself*/
                            std::vector<AlignmentRecord> splitAlns;
                            splitRecord(curRec, maxGapLength, minAlnLength, splitAlns);
                            for (auto aln : splitAlns) {
                                    AlignmentRecord *a = new AlignmentRecord(aln);
                                    AlignmentRecord *rev = aln.revert();
                                    a->sym = rev;
                                    rev->sym = a;
                                    result.push_back(a);
                                    result.push_back(rev);
                            }
                    }
                    pslFile.close();
                    std::cerr << "Done." << std::endl;
            }
            else {
                    std::cerr << "ERROR: psl file could not be opened!" << std::endl;
                    exit(EXIT_FAILURE);
            }
        }
}

void fillBuckets(std::vector<AlignmentRecord *>& alns, unsigned int bucketSize,
	std::vector<std::vector<AlignmentRecord *>>& result) {
	unsigned int firstBucket, lastBucket;
	for (auto alnPtr : alns) {
		firstBucket = alnPtr->tStart / bucketSize;
		lastBucket = alnPtr->tEnd / bucketSize;
		for (auto i = firstBucket; i <= lastBucket; i++)
			result[i].push_back(alnPtr);
	}
}

