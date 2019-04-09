#pragma once

#include <vector>
#include <map>
#include <chrono>
#include "AlignmentRecord.h"

/* Prints result */
void printResult(const std::vector<WasteRegion>&,
	const std::vector<int>&, 
	const std::map<std::string, unsigned long>&);

/* Prints the time elapsed since beginning */
void shoutTime(const std::chrono::time_point<std::chrono::high_resolution_clock>);
