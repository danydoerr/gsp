#pragma once

#include <map>
#include "AlignmentRecord.h"

/* Finds connected components. */
void classify(const std::vector<WasteRegion> &regions,
	const std::vector<std::vector<AlignmentRecord *>> &buckets,
	unsigned int bucketSize, float minAlnCoverage,
	std::vector<int> &classes, int &nrClasses);