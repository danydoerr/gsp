#pragma once

#include <map>
#include "AlignmentRecord.h"

/* Finds connected components. */
void classify(const std::vector<WasteRegion> &regions,
	const std::vector<std::vector<std::shared_ptr<AlignmentRecord>>> &buckets,
	unsigned int bucketSize, float minAlnCoverage,
	std::vector<int> &classes);