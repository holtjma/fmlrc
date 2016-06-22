#ifndef ALIGNMENT_UTIL_H
#define ALIGNMENT_UTIL_H

//C headers
#include <stdint.h>

//C++ headers
#include <vector>

using namespace std;

uint64_t editDistance(const vector<uint8_t> &s1,  const vector<uint8_t> &s2);
pair<uint64_t, uint64_t> editDistance_minimize(const vector<uint8_t> &s1,  const vector<uint8_t> &s2);

#endif