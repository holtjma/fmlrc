
//C headers
#include <stdint.h>

//C++ headers
#include <fstream>
#include <string>
#include <vector>

//my headers
#include "file_iterators.h"

using namespace std;

FastaIterator::FastaIterator(string fastaFN) {
    //open the file
    this->fastaFN = fastaFN;
    this->ifp.open(this->fastaFN);
    
    //read the first line so we're in the correct state
    this->isMoreData = (bool)getline(this->ifp, this->nextLine);
}

struct LongReadFA FastaIterator::getNextRead() {
    struct LongReadFA ret;
    ret.label = this->nextLine;
    vector<string> seqFrags = vector<string>();
    uint64_t seqLen = 0;
    
    while(getline(this->ifp, this->nextLine)) {
        if(this->nextLine[0] == '>') {
            //put the string together and return
            ret.seq.resize(seqLen);
            uint64_t currPos = 0;
            for(uint64_t x = 0; x < seqFrags.size(); x++) {
                ret.seq.replace(currPos, seqFrags[x].size(), seqFrags[x]);
                currPos += seqFrags[x].size();
            }
            return ret;
        }
        else {
            //push back a fragments
            seqFrags.push_back(this->nextLine);
            seqLen += this->nextLine.size();
        }
    }
    
    //we hit the last line
    this->isMoreData = false;
    
    //put the string together and return
    ret.seq.resize(seqLen);
    uint64_t currPos = 0;
    for(uint64_t x = 0; x < seqFrags.size(); x++) {
        ret.seq.replace(currPos, seqFrags[x].size(), seqFrags[x]);
        currPos += seqFrags[x].size();
    }
    return ret;
}

FastaWriter::FastaWriter(string fastaFN, int symsPerLine) {
    this->fastaFN = fastaFN;
    this->symsPerLine = symsPerLine;
    
    this->ofp.open(this->fastaFN);
}

FastaWriter::~FastaWriter() {
    this->ofp.close();
}

bool FastaWriter::writeRead(LongReadFA r) {
    this->ofp << r.label << "\n";
    for(uint64_t x = 0; x < r.seq.size(); x += this->symsPerLine) {
        this->ofp << r.seq.substr(x, this->symsPerLine) << "\n";
    }
    return true;
}
