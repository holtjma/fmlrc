#ifndef BIT_ARRAY_H
#define BIT_ARRAY_H

//C headers
#include <stdint.h>

//C++ headers
#include <vector>

using namespace std;

inline void setBit64(vector<uint64_t> &bitArray, uint64_t index) {
    bitArray[index >> 6] |= ((uint64_t)0x1 << (index & 0x3F));
}

//returns the number of set bits in value
inline uint64_t rank64(uint64_t value) {
    uint64_t ret = value-((value & (uint64_t)0xAAAAAAAAAAAAAAAA) >> 1);
    ret = (ret & (uint64_t)0x3333333333333333) + ((ret >> 2) & (uint64_t)0x3333333333333333);
    ret = (ret + (ret >> 4)) & (uint64_t)0x0F0F0F0F0F0F0F0F;
    return (ret * (uint64_t)0x0101010101010101) >> 56;
}

class BitArray {
private:
    uint64_t numValues;
    vector<uint64_t> ba;
    vector<uint64_t> index;
public:
    //constructor
    BitArray(uint64_t baLen);
    
    //use this to fill in the array with values
    inline void setBit(uint64_t ind) {
        setBit64(this->ba, ind);
    }
    
    //once filled in, use this to build the offsets
    void createIndex(uint64_t initialRank=0);
    
    //currently returns the number of set bits up to and including "pos"
    inline uint64_t rank(uint64_t pos) {
        //up to and including "pos"
        //return this->index[pos >> 6] + rank64(this->ba[pos >> 6] << (~pos & 0x3F));
        //up to but NOT including "pos"
        return this->index[pos >> 6] + rank64((this->ba[pos >> 6] << (~pos & 0x3F)) << 1);
    }
    
    //TODO: do we need this for our purposes? I don't think so right now
    //inline select(uint64_t rank);
};

#endif