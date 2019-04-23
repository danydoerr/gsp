#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <string>

#include "Util.h"
#include "InputParser.h"

int main(int argc, char *argv[]) {
    unsigned long max_bsize, max_start;
    
    InputParser parser;
    parser.parseCmdArgs(argc, argv);
    parser.getMaxBlockSizeAndLocalStart(max_bsize, max_start);
    std::cerr << "For values that fit, MAX SIZE: " << max_bsize << ", MAX START: " << max_start << std::endl;
    
    return 0;
}