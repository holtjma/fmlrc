
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

RLE_BWT::RLE_BWT(string inFN, uint8_t bitPower) {
    this->bwtFN = inFN;
    this->bitPower = bitPower;
    this->binSize = pow(2, this->bitPower);
    
    //get the bwt file size
    struct stat bwt_st;
    stat(this->bwtFN.c_str(), &bwt_st);
    
    //get ready to read the bwt
    ifstream bwtIn(this->bwtFN, ios::in | ios::binary);
    this->bwt = vector<uint8_t>(1000);
    
    //read the numpy header: http://docs.scipy.org/doc/numpy-1.10.1/neps/npy-format.html
    bwtIn.read((char*)&this->bwt[0], 16);
    uint16_t headerLen = this->bwt[8]+256*(int)this->bwt[9];
    uint16_t skipBytes = 10+headerLen;
    if(skipBytes % 16 != 0) {
        skipBytes = ((skipBytes / 16)+1)*16;
    }
    bwtIn.read((char*)&this->bwt[0], skipBytes-16);
    
    //finally allocate and load the bwt
    uint64_t bwtDiskSize = bwt_st.st_size - skipBytes;
    this->bwt = vector<uint8_t>(bwtDiskSize);
    bwtIn.read((char*)&this->bwt[0], bwtDiskSize);
    bwtIn.close();
    printf("loaded bwt with %lu compressed values\n", this->bwt.size());
    
    //first get the total symbol counts
    this->constructTotalCounts();
    
    //build some auxiliary indices
    this->constructIndexing();
    
    //now build the FM-index
    this->constructFMIndex();
}

void RLE_BWT::constructFMIndex() {
    //figure out the number of entries and pre-allocate
    //uint64_t samplingSize = (uint64_t)ceil(((float)this->totalSize+1)/this->binSize);
    uint64_t samplingSize = (uint64_t)ceil(((float)this->totalSize+1)/this->binSize)+1;
    this->fmIndex = new uint64_t*[VC_LEN];
    for(int x = 0; x < VC_LEN; x++) {
        this->fmIndex[x] = new uint64_t[samplingSize];
    }
    this->refFM = vector<uint64_t>(samplingSize, 0);
    
    uint8_t prevChar = 0;
    uint64_t totalCharCount = 0;
    uint64_t powerMultiple = 1;
    uint64_t binEnd = 0;
    uint64_t binID = 0;
    uint64_t bwtIndex = 0;
    uint64_t prevStart = 0;
    uint8_t currentChar;
    
    vector<uint64_t> countsSoFar = vector<uint64_t>(VC_LEN);
    for(int x = 0; x < VC_LEN; x++) {
        countsSoFar[x] = this->startIndex[x];
    }
    
    //go through each run in the BWT and set FM-indices as we go
    uint64_t numBytes = this->bwt.size();
    for(uint64_t x = 0; x < numBytes; x++) {
        currentChar = this->bwt[x] & MASK;
        if(currentChar == prevChar) {
            totalCharCount += (this->bwt[x] >> LETTER_BITS) * powerMultiple;
            powerMultiple *= NUM_POWER;
        }
        else {
            //first save the current FM-index entry
            while(bwtIndex+totalCharCount >= binEnd) {
                this->refFM[binID] = prevStart;
                for(int y = 0; y < VC_LEN; y++) {
                    this->fmIndex[y][binID] = countsSoFar[y];
                }
                binEnd += this->binSize;
                binID++;
            }
            
            //now add the previous
            countsSoFar[prevChar] += totalCharCount;
            bwtIndex += totalCharCount;
            
            prevChar = currentChar;
            prevStart = x;
            totalCharCount = this->bwt[x] >> LETTER_BITS;
            powerMultiple = NUM_POWER;
        }
    }
    
    while(bwtIndex+totalCharCount >= binEnd) {
        this->refFM[binID] = prevStart;
        for(int y = 0; y < VC_LEN; y++) {
            this->fmIndex[y][binID] = countsSoFar[y];
        }
        binEnd += this->binSize;
        binID++;
    }
    
    //set the last entry
    countsSoFar[prevChar] += totalCharCount;//forces countSoFar to hold the very end FM-index entry
    this->refFM[samplingSize-1] = numBytes; //need to point to the index at the end
    for(int y = 0; y < VC_LEN; y++) {
        this->fmIndex[y][samplingSize-1] = countsSoFar[y];
    }
    
    //calculate the total offsetSum
    this->offsetSum = 0;
    for(int x = 0; x < VC_LEN; x++) {
        this->offsetSum += this->fmIndex[x][0];
    }
}

RLE_BWT::~RLE_BWT() {
    for(int x = 0; x < VC_LEN; x++) {
        delete this->fmIndex[x];
    }
    delete this->fmIndex;
}

bwtRange RLE_BWT::constrainRange(uint8_t sym, bwtRange inRange) {
    //first find the low value
    uint64_t binID = inRange.l >> this->bitPower;
    uint64_t compressedIndex = this->refFM[binID];
    uint64_t bwtIndex = 0;
    for(uint64_t x = 0; x < VC_LEN; x++) {
        bwtIndex += this->fmIndex[x][binID];
    }
    bwtIndex -= this->offsetSum;
    
    bwtRange ret;
    ret.l = this->fmIndex[sym][binID];
    
    /*
     Dear future Matt,
     You have already tried using shifts (<<, >>) instead of multiplication, no effect (if anything it was worse).
     Sincerely,
     Pass Matt
     */
    uint8_t prevChar = 255;
    uint8_t currentChar;
    uint64_t prevCount = 0;
    uint64_t powerMultiple = 1;
    
    while(bwtIndex+prevCount < inRange.l) {
        currentChar = this->bwt[compressedIndex] & MASK;
        if(currentChar == prevChar) {
            prevCount += (this->bwt[compressedIndex] >> LETTER_BITS) * powerMultiple;
            powerMultiple *= NUM_POWER;
        }
        else {
            if(prevChar == sym) ret.l += prevCount;
            
            bwtIndex += prevCount;
            prevCount = this->bwt[compressedIndex] >> LETTER_BITS;
            prevChar = currentChar;
            powerMultiple = NUM_POWER;
        }
        compressedIndex++;
    }
    
    uint64_t tempC = ret.l;
    if(prevChar == sym) ret.l += inRange.l - bwtIndex;
    
    //now find the high value
    uint64_t binID_h = inRange.h >> this->bitPower;
    if(binID == binID_h) ret.h = tempC;
    else {
        compressedIndex = this->refFM[binID_h];
        bwtIndex = 0;
        for(uint64_t x = 0; x < VC_LEN; x++) {
            bwtIndex += this->fmIndex[x][binID_h];
        }
        bwtIndex -= this->offsetSum;
        
        ret.h = this->fmIndex[sym][binID_h];
        
        prevChar = 255;
        prevCount = 0;
        powerMultiple = 1;
    }
    
    while(bwtIndex+prevCount < inRange.h) {
        currentChar = this->bwt[compressedIndex] & MASK;
        if(currentChar == prevChar) {
            prevCount += (this->bwt[compressedIndex] >> LETTER_BITS) * powerMultiple;
            powerMultiple *= NUM_POWER;
        }
        else {
            if(prevChar == sym) ret.h += prevCount;
            
            bwtIndex += prevCount;
            prevCount = this->bwt[compressedIndex] >> LETTER_BITS;
            prevChar = currentChar;
            powerMultiple = NUM_POWER;
        }
        compressedIndex++;
    }
    
    if(prevChar == sym) ret.h += inRange.h - bwtIndex;
    return ret;
}
