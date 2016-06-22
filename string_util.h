#ifndef STRING_UTIL_H
#define STRING_UTIL_H

//C headers
#include <stdint.h>

//C++ headers
#include <string>
#include <vector>

using namespace std;

namespace string_util {
    extern vector<uint8_t> STRING_TO_INT;
    extern vector<uint8_t> INT_TO_STRING;
    extern vector<uint8_t> REV_COMP_I;
    
    //ALWAYS CALL THIS ONCE
    void initializeStringUtil();
    
    //utilities
    inline vector<uint8_t> reverseComplement_i(vector<uint8_t> seq) {
        uint64_t seqLen = seq.size();
        vector<uint8_t> ret = vector<uint8_t>(seqLen);
        for(uint64_t x = 0; x < seqLen; x++) {
            ret[x] = REV_COMP_I[seq[seqLen-x-1]];
        }
        return ret;
    };
    
    inline vector<uint8_t> stoi(string seq) {
        uint64_t seqLen = seq.size();
        vector<uint8_t> ret(seqLen);
        for(uint64_t x = 0; x < seqLen; x++) {
            ret[x] = STRING_TO_INT[seq[x]];
        }
        return ret;
    };
    
    inline string itos(vector<uint8_t> seq_i) {
        uint64_t seqLen = seq_i.size();
        string ret(seqLen, ' ');
        for(uint64_t x = 0; x < seqLen; x++) {
            ret[x] = INT_TO_STRING[seq_i[x]];
        }
        return ret;
    };
}

#endif