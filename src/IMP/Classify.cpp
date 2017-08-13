#include <iostream>
#include <algorithm>
#include "Classify.h"
#include "IMP.h"

void chooseAtom(const std::vector<WasteRegion>& regions,
	const Region mappedAtom, unsigned int regionFirst, unsigned int regionLast,
	Region &atomResult, unsigned int &jResult) {
	unsigned int maxJ = 0;
	int maxLength = 0;
	for (auto j = regionFirst; j < regionLast; j++) {
		int newLength;
		if (j == regionFirst)
			newLength = regions[j+1].first - mappedAtom.first;
		else if (j + 1 != regionLast)
			newLength = regions[j + 1].first - regions[j].last + 1;
		else { // j == regionLast - 1
			if (regions[j + 1].last < mappedAtom.last)
				newLength = mappedAtom.last - regions[j + 1].last;
			else newLength = 0;
		}
		if (newLength > maxLength) {
			maxLength = newLength;
			maxJ = j;
		}
	}
	if (!maxLength) {
		std::cerr << "ERROR: Graph construction failed! maxLength is still 0 at the end of chooseAtom.";
		exit(EXIT_FAILURE);
	}
	jResult = maxJ;
	atomResult = Region(regions[maxJ].last, regions[maxJ+1].first);
}

/* Returns the portion of the atom that is covered by the interval [sndStart,sndEnd]. */
float coverage(const Region atom, unsigned long sndStart, unsigned long sndEnd) {
	unsigned long length = atom.getLength() - 1;
	if (length == 0) length = 1;
	long lastStart = std::max(atom.first, sndStart);
	long firstEnd = std::min(atom.last, sndEnd);
	float result = static_cast<float>(firstEnd - lastStart) / static_cast<float>(length);
	return result;
}

/* Connects atoms only if they are aligned to each other and exceed minAlnCoverage. */
void constructAtomGraph(const std::vector<WasteRegion>& regions,
	const std::vector<std::vector<std::shared_ptr<AlignmentRecord>>>& buckets,
	unsigned int bucketSize, float minAlnCoverage,
	std::vector<std::map<unsigned int, int>> &graph) {
	for (size_t i = 0; i < regions.size() - 1; i++) {
		Region atom(regions[i].last, regions[i+1].first);
		unsigned long bucketIdx = atom.getMiddlePos() / bucketSize;
		for (auto aln : buckets[bucketIdx]) { // iterate over alignments that could cover atom
			if (aln->tStart > atom.first || aln->tEnd < atom.last) continue; // alignment doesn't cover atom
			Region mappedAtom = mapAtomThroughAln(atom, *aln);
			auto regionFirst = binSearchRegion(mappedAtom.first, regions);
			auto regionLast = binSearchRegion(mappedAtom.last, regions);
			unsigned int jfinal;
			Region newAtom(0,0);
			if (regionFirst == regionLast) {
				newAtom = Region(regions[regionFirst].last, regions[regionFirst+1].first);
				jfinal = regionFirst;
			} else if (regionFirst == regionLast - 1) {
				if (mappedAtom.last <= regions[regionLast].last) {
					newAtom = Region(regions[regionFirst].last, regions[regionFirst + 1].first);
					jfinal = regionFirst;
				} else
					chooseAtom(regions, mappedAtom, regionFirst, regionLast, newAtom, jfinal);
			} else
				chooseAtom(regions, mappedAtom, regionFirst, regionLast, newAtom, jfinal);
			if (coverage(newAtom, aln->qStart, aln->qEnd) < minAlnCoverage) continue; // coverage too low
			if (coverage(newAtom, mappedAtom.first, mappedAtom.last) <= 0.0f) continue;
			if (coverage(mappedAtom, newAtom.first, newAtom.last) <= 0.0f) continue;
			if (coverage(newAtom, aln->tStart, aln->tEnd) >= minAlnCoverage
				&& coverage(atom, aln->qStart, aln->qEnd) >= minAlnCoverage) continue; // text and query cover both atoms
			signed char strand = (aln->strand == '+') ? 1 : -1;
			if (graph[i].count(jfinal))
				graph[i].find(jfinal)->second += strand;
			else graph[i].insert(std::make_pair(jfinal, strand));
			if (graph[jfinal].count(i))
				graph[jfinal].find(i)->second += strand;
			else graph[jfinal].insert(std::make_pair(i, strand));
		}
	}
	// remove unnecessary (i.e. empty) maps at back
	while (!graph.empty() && graph.back().empty())
		graph.pop_back();
}

/* Sets class of each atom in the connected component to classNr. */
void fillComponent(const std::vector<std::map<unsigned int, int>> &graph,
	std::vector<int> &classes, unsigned int i, int classNr) {
	if (classes[i]) {
		if (classes[i] != classNr && classes[i] != -classNr) {
			std::cerr << "ERROR: Bad class when filling component. Class is "
				<< classes[i] << ", but should be " << classNr;
			exit(EXIT_FAILURE);
		}
		return;
	}
	classes[i] = classNr;
	for (auto m : graph[i]) {
		auto other = classNr;
		if (m.second < 0) other = -classNr;
		fillComponent(graph, classes, m.first, other);
	}
}

void classify(const std::vector<WasteRegion>& regions,
	const std::vector<std::vector<std::shared_ptr<AlignmentRecord>>>& buckets,
	unsigned int bucketSize, float minAlnCoverage,
	std::vector<int> &classes, int &classNr) {
	if (regions.size() < 2) {
		std::cerr << "ERROR: Too small input region vector for atom classification.";
		return;
	}
	std::vector<std::map<unsigned int, int>> graph(regions.size() - 1);
	constructAtomGraph(regions, buckets, bucketSize, minAlnCoverage, graph);
	classes.resize(regions.size() - 1, 0);
	classNr = 0;
	for (size_t i = 0; i < regions.size()-1; i++) {
		if (!classes[i]) {
			classNr++;
			fillComponent(graph, classes, i, classNr);
		}
	}
}