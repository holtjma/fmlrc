
//C headers
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <sys/stat.h>

//C++ headers
#include <fstream>
#include <string>
#include <vector>

//Custom headers
#include "string_util.h"
#include "rle_bwt.h"

using namespace std;

BaseBWT::BaseBWT() {
    
}

void BaseBWT::constructTotalCounts() {
    //the setup
    this->totalCounts = vector<uint64_t>(VC_LEN, 0);
    uint8_t prevChar = 255;
    uint8_t currentChar;
    uint64_t powerMultiple = 1;
    uint64_t bwtSize = this->bwt.size();
    uint64_t currentCount;
    
    //go through each run and add the symbol counts
    for(uint64_t x = 0; x < bwtSize; x++) {
        currentChar = this->bwt[x] & MASK;
        if(currentChar == prevChar) {
            powerMultiple *= NUM_POWER;
        }
        else {
            powerMultiple = 1;
        }
        prevChar = currentChar;
        currentCount = (this->bwt[x] >> LETTER_BITS)* powerMultiple;
        this->totalCounts[currentChar] += currentCount;
    }
}

void BaseBWT::constructIndexing() {
    this->startIndex = vector<uint64_t>(VC_LEN, 0);
    this->endIndex = vector<uint64_t>(VC_LEN, 0);
    
    uint64_t pos = 0;
    for(uint64_t x = 0; x < VC_LEN; x++) {
        this->startIndex[x] = pos;
        pos += this->totalCounts[x];
        this->endIndex[x] = pos;
    }
    this->totalSize = pos;
}

BaseBWT::~BaseBWT() {
    
}


uint64_t BaseBWT::countKmer(uint8_t * kmer, uint64_t kmerSize) {
    bwtRange ret;
    ret.l = 0;
    ret.h = this->totalSize;
    
    for(int64_t x = kmerSize-1; x >= 0 && ret.l != ret.h; x--) {
        ret = this->constrainRange(kmer[x], ret);
    }
    
    return ret.h-ret.l;
}

vector<uint64_t> BaseBWT::countPileup_i(vector<uint8_t> seq, uint64_t kmerSize) {
    uint64_t seqLen = seq.size();
    if(seqLen < kmerSize) {
        return vector<uint64_t>(0);
    }
    
    uint64_t numCounts = seqLen-kmerSize+1;
    vector<uint64_t> ret = vector<uint64_t>(numCounts);
    
    vector<uint8_t> revComp = string_util::reverseComplement_i(seq);
    
    for(uint64_t x = 0; x < numCounts; x++) {
        ret[x] = this->countKmer(&seq[x], kmerSize)+this->countKmer(&revComp[seqLen-kmerSize-x], kmerSize);
    }
    return ret;
}
