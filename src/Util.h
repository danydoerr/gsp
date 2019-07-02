#pragma once

#include <vector>
#include <map>
#include <chrono>
#include <string>
#include "AlignmentRecord.h"

/* Returns index of the last element in xList that is <= x.
If all elements in xList are > x, result is 0. Expects xList to be sorted ascending. */
unsigned int binSearch(unsigned long x, const std::vector<unsigned long>& xList);

/* Prints result */
void printResult(const std::vector<WasteRegion>&,
	const std::vector<int>&, 
	const std::map<std::string, unsigned long>&);

/* Prints the time elapsed since beginning */
void shoutTime(const std::chrono::time_point<std::chrono::high_resolution_clock>);

/* Converts string to unsigned int, throwing an exception if the number doesn't fit. */
inline unsigned int stoui(const std::string& s)
{
        unsigned long lresult = stoul(s, 0, 10);
        unsigned int result = lresult;
        if (result != lresult) throw std::range_error("Cannot fit this number in an unsigned int: " + s + " (" + __FILE__ + ":" + std::to_string(__LINE__) + ")");
        return result;
}
/* Converts string to unsigned short, throwing an exception if the number doesn't fit. */
inline unsigned short stouh(const std::string& s)
{
        unsigned long lresult = stoul(s, 0, 10);
        unsigned short result = lresult;
        if (result != lresult) throw std::range_error("Cannot fit this number in an unsigned short: " + s + " (" + __FILE__ + ":" + std::to_string(__LINE__) + ")");
        return result;
}
