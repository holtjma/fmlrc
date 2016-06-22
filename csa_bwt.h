#ifndef CSA_BWT_H
#define CSA_BWT_H

//C headers
#include <stdint.h>

//C++ headers

//custom headers
#include "bit_array.h"
#include "base_bwt.h"

using namespace std;

class CSA_BWT : public BaseBWT{
private:
    //constructFMIndex()
    vector<BitArray*> csa;
    
    //these functions build all auxiliary structures required for the FM-index lookups
    void constructFMIndex(bool storeDN);
    
public:
    //constructor and destructor
    CSA_BWT(string inFN, bool storeD=true);
    ~CSA_BWT();
    
    //query sub-routines
    bwtRange constrainRange(uint8_t sym, bwtRange inRange);
};

#endif