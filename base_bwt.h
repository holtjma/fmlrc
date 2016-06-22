#ifndef BASE_BWT_H
#define BASE_BWT_H

//C headers
#include <stdint.h>

//C++ headers
#include <string>
#include <vector>

using namespace std;

enum {
    VC_LEN = 6,//$ A C G N T
    
    LETTER_BITS = 3, //defined
    NUMBER_BITS = 5, //8-letterBits
    NUM_POWER = 32,  //2**numberBits
    MASK = 7,        //255 >> numberBits
    
    //These used to be pre-defined, but are set up as user options now
    //BIT_POWER = 8,  //defined
    //BIN_SIZE = 256  //2**self.bitPower
};

struct bwtRange {
    uint64_t l;
    uint64_t h;
};

class BaseBWT {
protected:
    //loaded from disk
    string bwtFN;
    vector<uint8_t> bwt;
    
    //constructTotalCounts()
    vector<uint64_t> totalCounts;
    
    //constructIndexing()
    vector<uint64_t> startIndex;
    vector<uint64_t> endIndex;
    uint64_t totalSize;
    
    //these functions build all auxiliary structures required for the FM-index lookups
    void constructTotalCounts();
    void constructIndexing();
    
public:
    //constructor and destructor
    BaseBWT();
    ~BaseBWT();
    
    //basic query functions
    uint64_t countKmer(uint8_t * kmer, uint64_t kmerSize);
    
    //multi-query functions
    vector<uint64_t> countPileup_i(vector<uint8_t> seq, uint64_t kmerSize);
    
    //query sub-routines
    virtual bwtRange constrainRange(uint8_t sym, bwtRange inRange) = 0;
};

#endif