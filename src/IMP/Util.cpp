#include <iostream>
#include "Util.h"
#include "IMP.h"

void printResult(const std::vector<WasteRegion> &regions, const std::vector<int> &classes,
	const std::map<std::string, unsigned long> &speciesStarts) {
	// flipped is speciesStarts sorted by position
	std::map<unsigned long, std::string> flipped;
	for (auto specStart : speciesStarts)
		flipped[specStart.second] = specStart.first;
	std::vector<std::string> names;
	std::vector<unsigned long> starts;
	for (auto specStart : flipped) {
		starts.push_back(specStart.first);
		names.push_back(specStart.second);
	}

	std::cout << "#name\tatom_nr\tclass\tstrand\tstart\tend" << "\n"; // header line
	for (size_t i = 0; i + 1 < regions.size(); i++) {
		auto j = binSearch(regions[i].last, starts);
		auto move = starts[j];
		auto start = (regions[i].last > move) ? regions[i].last - move : 0;
		auto end = regions[i + 1].first - move;
		if (end > starts[j + 1]) end = starts[j + 1];
		auto classNr = abs(classes[i]);
		auto strand = (classes[i] > 0) ? '+' : '-';
		std::cout << names[j] << "\t" << i+1 << "\t" << classNr << "\t" << strand << "\t"
			<< start << "\t" << end << "\n";
	}
}

void shoutTime(const std::chrono::time_point<std::chrono::high_resolution_clock> start) {
	auto end = std::chrono::high_resolution_clock::now();
	auto diff = end - start;
	auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();
	std::cerr << " Time passed since start: " << ms << " milliseconds." << std::endl;
}
