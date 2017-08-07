#include "InputParser.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm>

void parseCmdArgs(int argc, char** &argv, char* &pslPath,
	unsigned int &minLengh, unsigned int &maxGap, unsigned int &minAlnLength,
	float &minAlnIdentity) {
	if (argc <= 1) {
		std::cerr << "Usage: Call with at least one arguments:\n"
			<< "atomizer <psl file> [options]\n\n"
			<< "Optional arguments are given after their descriptor. The descriptor is NOT case-sensitive. \n"
			<< "If an optional argument is not given, the default value will be used.\n"
			<< "--minLength <minLength>: The minimum length an atom must have (defualt: 250).\n"
			<< "--minIdent <minIdent>: Minimum identity an alignment must have to be considered, "
			<< "alignments with lower identity will be skipped. Give as percent, e.g. --minIdent 100 if "
			<< "you only want perfect alignments (default: 80).\n"
			<< "--maxGap <maxGap>: The maximum gap length inside of an alignment. If this length is "
			<< "exceeded, the alignment will be split into two (default: 13).\n"
			<< "--minAlnLength <minAlnLength>: The minimal length an alignment must have to be considered. "
			<< "Shorter alignments are droppe. Note that this is applied after the splitting at long gaps, "
			<< "so a sufficiently long alignment might be split into two that are two short each and "
			<< "be completely omitted (default:13)."
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
	pslPath = argv[1];
	for (int i = mandatoryArgs + 1; i < argc; i++) {
		std::string arg = argv[i];
		std::transform(arg.begin(), arg.end(), arg.begin(), tolower);
		try {
			if (arg == "--minlength") minLengh = std::stoul(argv[++i]);
			else if (arg == "--minident") minAlnIdentity = std::stoul(argv[++i]) / 100.0f;
			else if (arg == "--maxgap") maxGap = std::stoul(argv[++i]);
			else if (arg == "--minalnlength") minAlnLength = std::stoul(argv[++i]);
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

/* Parses a single psl line to an AlignmentRecord. */
const AlignmentRecord recordFromPsl(const std::vector<std::string>& pslLine,
	std::map<std::string, unsigned long>& speciesStart) {
	const std::string strand = pslLine[8], qName = pslLine[9], tName = pslLine[13];
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
	auto qStart = std::stoul(pslLine[11]) + qOffset, qEnd = std::stoul(pslLine[12]) + qOffset,
		tStart = std::stoul(pslLine[15]) + tOffset, tEnd = std::stoul(pslLine[16]) + tOffset,
		blockCount = std::stoul(pslLine[17]);
	auto qStarts = stringVecToLongVec(splitString(pslLine[19], ',', blockCount)),
		tStarts = stringVecToLongVec(splitString(pslLine[20], ',', blockCount));
	// shift start positions by offest. If on reverse strand, recompute query w.r.t. starts to start of sequence
	if (strand == "+")
		for (auto i = qStarts.begin(); i != qStarts.end(); i++)
			*i += qOffset;
	else {
		auto qSize = std::stoul(pslLine[10]);
		for (auto i = qStarts.begin(); i != qStarts.end(); i++)
			*i = qSize - *i + qOffset;
	}
	for (auto i = tStarts.begin(); i != tStarts.end(); i++)
		*i += tOffset;

	return AlignmentRecord(strand.front(), qStart, qEnd, tStart, tEnd, blockCount,
		stringVecToLongVec(splitString(pslLine[18], ',', blockCount)), // blockSizes
		qStarts, tStarts);
}

/* Returns a part of the input AlignmentRecord, from startBlock to endBlock, including both. */
const AlignmentRecord cutRecord(const AlignmentRecord &aln, unsigned int startBlock, unsigned int endBlock) {
	if (startBlock == 0 && endBlock == aln.blockCount - 1) return aln;
	unsigned long qStart, qEnd, tStart = aln.tStarts[startBlock],
		tEnd = aln.tStarts[endBlock] + aln.blockSizes[endBlock];
	if (aln.strand == '+') {
		qStart = aln.qStarts[startBlock];
		qEnd = aln.qStarts[endBlock] + aln.blockSizes[endBlock];
	} else { // strand == '-'
		qStart = aln.qStarts[endBlock] - aln.blockSizes[endBlock];
		qEnd = aln.qStarts[startBlock];
	}
	unsigned int blockCount = endBlock - startBlock + 1;
	std::vector<unsigned long> blockSizes(aln.blockSizes.begin() + startBlock, aln.blockSizes.begin() + endBlock + 1),
		qStarts(aln.qStarts.begin() + startBlock, aln.qStarts.begin() + endBlock + 1),
		tStarts(aln.tStarts.begin() + startBlock, aln.tStarts.begin() + endBlock + 1);
	return AlignmentRecord(aln.strand, qStart, qEnd, tStart, tEnd, blockCount, blockSizes, qStarts, tStarts);
}

/* Splits input AlignmentRecord in parts if it contains gaps longer than maxGapLength .
Stores splitted parts in result only if they are longer than minAlnLength. */
void splitRecord(const AlignmentRecord &aln,
	unsigned int maxGapLength, unsigned int minAlnLength, std::vector<AlignmentRecord> &result) {
	int start = 0;
	for (size_t i = 0; i < aln.blockCount - 1; i++) {
		// check for long gaps
		if (aln.tStarts[i + 1] - (aln.tStarts[i] + aln.blockSizes[i]) > maxGapLength
			|| (aln.strand == '+' && aln.qStarts[i + 1] - (aln.qStarts[i] + aln.blockSizes[i]) > maxGapLength)
			|| (aln.strand == '-' && aln.qStarts[i] - (aln.qStarts[i + 1] + aln.blockSizes[i]) > maxGapLength)) {
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

void parsePsl(const char * pslPath, std::map<std::string, unsigned long>& speciesStart,
	unsigned int maxGapLength, unsigned int minAlnLength, float minAlnIdentity,
	std::vector<std::shared_ptr<AlignmentRecord>>& result) {
	std::ifstream pslFile;
	pslFile.open(pslPath);
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
			if (curRec.tStart > curRec.qStart) continue; // only one version of symmetric alignments
			if (curRec.tStart == curRec.qStart && curRec.tEnd == curRec.qEnd)
				continue; // skip regions that are aligned to themselves
			std::vector<AlignmentRecord> splitAlns;
			splitRecord(curRec, maxGapLength, minAlnLength, splitAlns);
			for (auto aln : splitAlns) {
				std::shared_ptr<AlignmentRecord> a(new AlignmentRecord(aln));
				std::shared_ptr<AlignmentRecord> rev(aln.revert());
				a->sym = rev;
				rev->sym = a;
				result.push_back(a);
				result.push_back(rev);
			}
		}
		pslFile.close();
	}
	else {
		std::cerr << "ERROR: psl file could not be opened!" << std::endl;
		exit(EXIT_FAILURE);
	}
}

void fillBuckets(std::vector<std::shared_ptr<AlignmentRecord>>& alns, unsigned int bucketSize,
	std::vector<std::vector<std::shared_ptr<AlignmentRecord>>>& result) {
	unsigned int firstBucket, lastBucket;
	for (auto alnPtr : alns) {
		firstBucket = alnPtr->tStart / bucketSize;
		lastBucket = alnPtr->tEnd / bucketSize;
		for (auto i = firstBucket; i <= lastBucket; i++)
			result[i].push_back(alnPtr);
	}
}
