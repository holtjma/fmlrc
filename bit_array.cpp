
//C headers
#include <math.h>

//custom headers
#include "bit_array.h"

using namespace std;

BitArray::BitArray(uint64_t baLen) {
    this->numValues = ceil(baLen/64.0);
    this->ba = vector<uint64_t>(this->numValues, 0);
}

void BitArray::createIndex(uint64_t initialRank) {
    this->index = vector<uint64_t>(this->numValues+1, 0);
    
    uint64_t offset = initialRank;
    for(uint64_t x = 0; x < this->numValues; x++) {
        this->index[x] = offset;
        offset += rank64(this->ba[x]);
    }
    this->index[this->numValues] = offset;
}