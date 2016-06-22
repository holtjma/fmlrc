
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
#include "csa_bwt.h"

using namespace std;

CSA_BWT::CSA_BWT(string inFN, bool storeD) {
    this->bwtFN = inFN;
    
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
    bwtIn.read((char*)&this->bwt[0], skipBytes - 16);
    
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
    this->constructFMIndex(storeD);
    
    //now delete the original BWT
    this->bwt = vector<uint8_t>(0);
}

void CSA_BWT::constructFMIndex(bool storeD) {
    //figure out the number of entries and pre-allocate
    
    this->csa = vector<BitArray*>(VC_LEN);
    for(int x = 0; x < VC_LEN; x++) {
        if(x != 0 || storeD) this->csa[x] = new BitArray(this->totalSize);
        else this->csa[x] = new BitArray(1); //this is just a dummy to get deleted later
    }
    
    uint8_t prevChar = 0;
    uint64_t totalCharCount = 0;
    uint64_t powerMultiple = 1;
    //uint64_t binEnd = 0;
    //uint64_t binID = 0;
    uint64_t bwtIndex = 0;
    //uint64_t prevStart = 0;
    uint8_t currentChar;
    
    //vector<uint64_t> setCount = vector<uint64_t>(VC_LEN, 0);
    
    //go through each run in the BWT and set FM-indices as we go
    uint64_t numBytes = this->bwt.size();
    for(uint64_t x = 0; x < numBytes; x++) {
        //printf("%d %d\n", x, numBytes);
        currentChar = this->bwt[x] & MASK;
        if(currentChar == prevChar) {
            totalCharCount += (this->bwt[x] >> LETTER_BITS) * powerMultiple;
            powerMultiple *= NUM_POWER;
        }
        else {
            //first save the current FM-index entry
            if(prevChar != 0 || storeD) {
                for(uint64_t y = bwtIndex; y < bwtIndex+totalCharCount; y++) {
                    this->csa[prevChar]->setBit(y);
                    //setCount[prevChar] += 1;
                }
            }
            
            //now add the previous
            bwtIndex += totalCharCount;
            prevChar = currentChar;
            totalCharCount = this->bwt[x] >> LETTER_BITS;
            powerMultiple = NUM_POWER;
        }
    }
    
    if(prevChar != 0 || storeD) {
        for(uint64_t y = bwtIndex; y < bwtIndex+totalCharCount; y++) {
            this->csa[prevChar]->setBit(y);
        }
    }
    
    for(uint64_t x = 0; x < VC_LEN; x++) {
        if(x != 0 || storeD) this->csa[x]->createIndex(this->startIndex[x]);
    }
}

CSA_BWT::~CSA_BWT() {
    for(int x = 0; x < VC_LEN; x++) {
        delete this->csa[x];
    }
}

bwtRange CSA_BWT::constrainRange(uint8_t sym, bwtRange inRange) {
    //first find the low value
    bwtRange ret;
    ret.l = this->csa[sym]->rank(inRange.l);
    ret.h = this->csa[sym]->rank(inRange.h);
    return ret;
};