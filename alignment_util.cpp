
//C headers
#include <stdint.h>

//C++ headers
#include <vector>

//my custom headers
#include "alignment_util.h"

using namespace std;

//This code is essentially an adaptation of the algorithms on wikipedia: https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#C.2B.2B

uint64_t editDistance(const vector<uint8_t> &s1,  const vector<uint8_t> &s2)
{
    //modified version of https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#C.2B.2B
    const std::size_t len1 = s1.size(), len2 = s2.size();
    std::vector<uint64_t> col(len2+1), prevCol(len2+1);
    
    for (unsigned int i = 0; i < prevCol.size(); i++)
        prevCol[i] = i;
    for (unsigned int i = 0; i < len1; i++) {
        col[0] = i+1;
        for (unsigned int j = 0; j < len2; j++)
            // note that std::min({arg1, arg2, arg3}) works only in C++11,
            // for C++98 use std::min(std::min(arg1, arg2), arg3)
            //col[j+1] = std::min({ prevCol[1 + j] + 1, col[j] + 1, prevCol[j] + (s1[i]==s2[j] ? 0 : 1) });
            col[j+1] = std::min(std::min(prevCol[1+j]+1, col[j]+1), prevCol[j]+(s1[i]==s2[j] ? 0 : 1));
        col.swap(prevCol);
    }
    return prevCol[len2];
}

pair<uint64_t, uint64_t> editDistance_minimize(const vector<uint8_t> &s1,  const vector<uint8_t> &s2)
{
    //modified version of https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#C.2B.2B
    const std::size_t len1 = s1.size(), len2 = s2.size();
    std::vector<uint64_t> col(len2+1), prevCol(len2+1);
    
    for (unsigned int i = 0; i < prevCol.size(); i++)
        prevCol[i] = i;
    for (unsigned int i = 0; i < len1; i++) {
        col[0] = i+1;
        for (unsigned int j = 0; j < len2; j++)
            // note that std::min({arg1, arg2, arg3}) works only in C++11,
            // for C++98 use std::min(std::min(arg1, arg2), arg3)
            //col[j+1] = std::min({ prevCol[1 + j] + 1, col[j] + 1, prevCol[j] + (s1[i]==s2[j] ? 0 : 1) });
            col[j+1] = std::min(std::min(prevCol[1+j]+1, col[j]+1), prevCol[j]+(s1[i]==s2[j] ? 0 : 1));
        col.swap(prevCol);
    }
    
    //get the last smallest value
    uint64_t argMin = len2;
    for(int64_t x = len2-1; x >= 0; x--) {
        if(prevCol[x] < prevCol[argMin]) argMin = x;
    }
    
    //pick the smallest score
    pair<uint64_t, uint64_t> ret;// = scoreArray[oLen][mLen][choice];
    ret.first = prevCol[argMin];
    ret.second = argMin;
    return ret;
}