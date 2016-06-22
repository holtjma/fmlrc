#ifndef RLE_BWT_H
#define RLE_BWT_H

//C headers
#include <stdint.h>

//C++ headers
#include <string>
#include <vector>

//Custom headers
#include "base_bwt.h"

using namespace std;

class RLE_BWT : public BaseBWT {
private:
    //loaded from disk
    uint8_t bitPower;
    uint64_t binSize;
    
    //constructFMIndex() - IMPORTANT: THIS IS TRANSPOSED COMPARED TO PYTHON IMPL
    //aka row = symbol; column = index
    uint64_t** fmIndex;
    vector<uint64_t> refFM;
    uint64_t offsetSum;
    
    //these functions build all auxiliary structures required for the FM-index lookups
    void constructFMIndex();
    
public:
    //constructor and destructor
    RLE_BWT(string inFN, uint8_t bitPower=8);
    ~RLE_BWT();
    
    //query sub-routines
    bwtRange constrainRange(uint8_t sym, bwtRange inRange);
};

#endif