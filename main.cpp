
//C headers
#include <stdint.h>
#include <stdio.h>
#include <unistd.h>

//C++ headers
#include <algorithm>
#include <cassert>
#include <numeric>
#include <vector>

//external custom headers
#include "CTPL/ctpl_stl.h"

//my custom headers
#include "alignment_util.h"
#include "base_bwt.h"
#include "csa_bwt.h"
#include "file_iterators.h"
#include "rle_bwt.h"
#include "string_util.h"

struct Parameters {
    bool USE_FM_INDEX;
    uint64_t k;
    uint64_t K;
    uint64_t MIN_COUNT;
    uint64_t MAX_BRANCH_ATTEMPT_LENGTH;
    double BRANCH_BUFFER_FACTOR;
    double TAIL_BUFFER_FACTOR;
    double FRAC;
    //uint64_t MAX_TRIES;
    uint8_t FM_BIT_POWER;
    bool VERBOSE;
};

struct CorrectionResults {
    string label;
    string originalSeq;
    string correctedSeq;
    double avgBefore;
    double avgAfter;
};

//valid char info
enum {
    VALID_CHARS_LEN = 4
};
const vector<uint8_t> VALID_CHARS = {1, 2, 3, 5};

const string VERSION = "0.1.0";

uint64_t calculateMedian(vector<uint64_t> inArray, uint64_t minValue) {
    /*
     Calculates the median of the array ignoring all values < minValue
     Note: this impl doesn't average the median if there are an even number of values, just picks the lower one
     @param inArray - the vector to calculate the median of
     @param minValue - any values less than this will be ignored in the calculation
     */
    
    uint64_t arrayLen = inArray.size();
    vector<uint64_t> arrayCopy = vector<uint64_t>(arrayLen, 0);
    uint64_t l = 0;
    
    for(uint64_t x = 0; x < arrayLen; x++) {
        if(inArray[x] >= minValue) {
            arrayCopy[l] = inArray[x];
            l++;
        }
    }
    
    if(l == 0) {
        return 0;
    }
    else {
        nth_element(arrayCopy.begin(), arrayCopy.begin()+(l-1)/2, arrayCopy.begin()+l);
        return arrayCopy[(l-1)/2];
    }
}

vector<vector<uint8_t> > multiBridge(BaseBWT * rle_p, vector<uint8_t> seedKmer, vector<uint8_t> targetKmer, uint64_t tMin, uint64_t branchLim, uint64_t maxBranchLen) {
    /*
    printf("multibridge ");
    for(int x = 0; x < seedKmer.size(); x++) printf("%d", seedKmer[x]);
    printf(" ");
    for(int x = 0; x < targetKmer.size(); x++) printf("%d", targetKmer[x]);
    printf("\n");
    */
    
    //printf("tMin = %d\n", tMin);
    
    vector<vector<uint8_t> > ret = vector<vector<uint8_t> >(0);
    uint64_t kmerLen = seedKmer.size();
    
    vector<uint64_t> counts = vector<uint64_t>(4);
    
    uint64_t numBranched = 0;
    
    // cdef str currBridge
    vector<uint8_t> currBridge;
    // cdef list currBridgeList
    uint64_t currBridgeLen = 0;
    vector<uint8_t> currKmer = vector<uint8_t>(kmerLen, 4);
    vector<uint8_t> revKmer = vector<uint8_t>(kmerLen, 4);
    
    // cdef list possBridges = [(seedKmer, kmerLen)]
    vector<vector<uint8_t> > possBridges = vector<vector<uint8_t> >();
    possBridges.push_back(vector<uint8_t>(seedKmer));
    
    // cdef unsigned char * currBridge_view
    // cdef unsigned char * currKmer_view = <bytes>currKmer
    // cdef unsigned char * revKmer_view = <bytes>revKmer
    // cdef unsigned char * targetKmer_view = <bytes>targetKmer
    // #print currKmer_view[0], ord('A')
    
    // cdef unsigned long i, x
    // cdef str c
    
    uint64_t maxPos;
    
    //while we have things to explore, and we haven't explored too many, and we don't have a ridiculous number of possibilities
    while(possBridges.size() > 0 && numBranched < branchLim) {
        currBridge = possBridges.back();
        possBridges.pop_back();
        currBridgeLen = currBridge.size();
        numBranched++;
        
        for(unsigned int x = 0; x < kmerLen; x++) {
            currKmer[x] = currBridge[currBridgeLen-kmerLen+x];
            revKmer[kmerLen-x-1] = string_util::REV_COMP_I[currKmer[x]];
        }
        
        //try to extend the bridge
        while(currBridgeLen < maxBranchLen) {
            //shift the current k-mer over one in preparation for the last base toggle
            for(unsigned int x = 0; x < kmerLen-1; x++) {
                currKmer[x] = currKmer[x+1];
                revKmer[kmerLen-x-1] = revKmer[kmerLen-x-2];
            }
            
            maxPos = 0;
            
            //count and pick the highest
            //printf("counts ");
            for(int x = 0; x < VALID_CHARS_LEN; x++) {
                currKmer[kmerLen-1] = VALID_CHARS[x];
                revKmer[0] = string_util::REV_COMP_I[VALID_CHARS[x]];
                counts[x] = rle_p->countKmer(&currKmer[0], kmerLen)+rle_p->countKmer(&revKmer[0], kmerLen);
                //printf("%d ", counts[x]);
                if(counts[x] > counts[maxPos]) maxPos = x;
            }
            //printf("\n");
            
            //make sure the highest is high enough for us to consider it
            if(counts[maxPos] >= tMin) {
                currBridge.push_back(4);
                
                if(possBridges.size() < branchLim) {
                    for(unsigned int x = 0; x < VALID_CHARS_LEN; x++) {
                        if(x != maxPos && counts[x] >= tMin) {
                            //add the ones we aren't exploring right now if they're high enough
                            currBridge[currBridgeLen] = VALID_CHARS[x];
                            possBridges.push_back(vector<uint8_t>(currBridge.begin(), currBridge.end()));
                        }
                    }
                }
                else {
                    //printf("exit A\n");
                    return vector<vector<uint8_t> >();
                }
                
                //now really add the symbol
                currBridge[currBridgeLen] = VALID_CHARS[maxPos];
                currBridgeLen++;
                
                currKmer[kmerLen-1] = VALID_CHARS[maxPos];
                revKmer[0] = string_util::REV_COMP_I[VALID_CHARS[maxPos]];
            }
            else {
                //our BEST doesn't pass the threshold on this path, stop following
                //print currBridge, counts, tMin
                break;
            }
            
            if(equal(targetKmer.begin(), targetKmer.end(), currKmer.begin())) {
                ret.push_back(currBridge);
                if(ret.size() >= branchLim) {
                    //printf("exit B\n");
                    return vector<vector<uint8_t> >();
                }
            }
        }
    }
    
    if(numBranched < branchLim) {
        return ret;
    }
    else {
        //printf("exit C\n");
        return vector<vector<uint8_t> >();
    }
}

vector<vector<uint8_t> > shortAssemble(BaseBWT * rle_p, vector<uint8_t> seedKmer, uint64_t tMin, uint64_t branchLim, uint64_t maxBranchLen) {
    vector<vector<uint8_t> > ret = vector<vector<uint8_t> >(0);
    uint64_t kmerLen = seedKmer.size();
    
    vector<uint64_t> counts = vector<uint64_t>(4);
    
    uint64_t numBranched = 0;
    
    vector<uint8_t> currBridge;
    uint64_t currBridgeLen = 0;
    vector<uint8_t> currKmer = vector<uint8_t>(kmerLen, 4);
    vector<uint8_t> revKmer = vector<uint8_t>(kmerLen, 4);
    
    vector<vector<uint8_t> > possBridges = vector<vector<uint8_t> >();
    possBridges.push_back(vector<uint8_t>(seedKmer));
    
    uint64_t maxPos;
    
    //while we have things to explore, and we haven't explored too many, and we don't have a ridiculous number of possibilities
    while(possBridges.size() > 0 && numBranched < branchLim) {
        currBridge = possBridges.back();
        possBridges.pop_back();
        currBridgeLen = currBridge.size();
        numBranched++;
        
        for(unsigned int x = 0; x < kmerLen; x++) {
            currKmer[x] = currBridge[currBridgeLen-kmerLen+x];
            revKmer[kmerLen-x-1] = string_util::REV_COMP_I[currKmer[x]];
        }
        
        //try to extend the bridge
        while(currBridgeLen < maxBranchLen) {
            //shift the current k-mer over one in preparation for the last base toggle
            for(unsigned int x = 0; x < kmerLen-1; x++) {
                currKmer[x] = currKmer[x+1];
                revKmer[kmerLen-x-1] = revKmer[kmerLen-x-2];
            }
            
            maxPos = 0;
            
            //count and pick the highest
            //printf("counts ");
            for(int x = 0; x < VALID_CHARS_LEN; x++) {
                currKmer[kmerLen-1] = VALID_CHARS[x];
                revKmer[0] = string_util::REV_COMP_I[VALID_CHARS[x]];
                counts[x] = rle_p->countKmer(&currKmer[0], kmerLen)+rle_p->countKmer(&revKmer[0], kmerLen);
                //printf("%d ", counts[x]);
                if(counts[x] > counts[maxPos]) maxPos = x;
            }
            //printf("\n");
            
            //make sure the highest is high enough for us to consider it
            if(counts[maxPos] >= tMin) {
                currBridge.push_back(4);
                
                if(possBridges.size() < branchLim) {
                    for(unsigned int x = 0; x < VALID_CHARS_LEN; x++) {
                        if(x != maxPos && counts[x] >= tMin) {
                            //add the ones we aren't exploring right now if they're high enough
                            currBridge[currBridgeLen] = VALID_CHARS[x];
                            possBridges.push_back(vector<uint8_t>(currBridge.begin(), currBridge.end()));
                        }
                    }
                }
                else {
                    //printf("exit A\n");
                    return vector<vector<uint8_t> >();
                }
                
                //now really add the symbol
                currBridge[currBridgeLen] = VALID_CHARS[maxPos];
                currBridgeLen++;
                
                currKmer[kmerLen-1] = VALID_CHARS[maxPos];
                revKmer[0] = string_util::REV_COMP_I[VALID_CHARS[maxPos]];
            }
            else {
                //our BEST doesn't pass the threshold on this path, stop following
                //print currBridge, counts, tMin
                break;
            }
        }
        
        if(currBridgeLen == maxBranchLen) {
            ret.push_back(currBridge);
            if(ret.size() >= branchLim) {
                return vector<vector<uint8_t> >();
            }
        }
    }
    
    //make sure we didn't go overboard
    if(numBranched < branchLim) {
        return ret;
    }
    else {
        return vector<vector<uint8_t> >();
    }
}

//this structure is used to store any corrections we find
struct Correction {
    uint64_t start;
    uint64_t end;
    vector<uint8_t> seq;
};

vector<uint8_t> correctionPass(BaseBWT * rle_p, vector<uint8_t> seq_i, Parameters myParams, uint64_t kmerSize) {
    /*
     This function performs a single pass of the correction algorithm
     @param rle_p - a pointer to the BWT that represents our DBG
     @param seq_i - the read sequence we are correcting in vector<uint8_t> form
     @param myParams - parameters that are modified by the user to effect these methods
     @param kmerSize - the 'k' in k-mer
     @return - a vector<uint8_t> representation of the string
     */
    
    //this is the only parameter that is dynamic right now
    //uint64_t BRANCH_LIMIT = 2*kmerSize;
    uint64_t BRANCH_LIMIT = 10*kmerSize;
    
    vector<uint64_t> pu = rle_p->countPileup_i(seq_i, kmerSize);
    
    double nzMed = calculateMedian(pu, myParams.MIN_COUNT);
    if(nzMed < myParams.MIN_COUNT) {
        //basically if our median is super low, we have no chance of fixing it
        return seq_i;
    }
    
    //try to dynamically set the threshold, but make sure its at least MIN_COUNT
    uint64_t thresh = (uint64_t)(myParams.FRAC * nzMed);
    if(thresh < myParams.MIN_COUNT) thresh = myParams.MIN_COUNT;
    
    //prep for the actual corrections now
    int64_t prevFound = -1;
    vector<uint8_t> seedKmer = vector<uint8_t>(kmerSize);
    vector<uint8_t> targetKmer = vector<uint8_t>(kmerSize);
    uint64_t maxBranchLength;
    vector<vector<uint8_t> > bridgePoints;
    vector<vector<uint8_t> > bridgePoints_ed = vector<vector<uint8_t> >();
    
    vector<Correction> correctionsList = vector<Correction>(0);
    Correction newCorr;
    
    uint64_t x = 0;
    uint64_t puSize = pu.size();
    
    while(x < puSize) {
        if(pu[x] < thresh) {
            prevFound = x-1;
            
            //find the next index that is above the threshold
            while(x < puSize && pu[x] < thresh) {
                x++;
            }
            
            if(prevFound == -1 && x < puSize) {
                //handle the head case
                maxBranchLength = (uint64_t)(myParams.TAIL_BUFFER_FACTOR*(x+kmerSize));
                if(maxBranchLength <= myParams.MAX_BRANCH_ATTEMPT_LENGTH) {
                    //get the first found k-mer and reverse complement it
                    seedKmer.assign(seq_i.begin()+x, seq_i.begin()+x+kmerSize);
                    seedKmer = string_util::reverseComplement_i(seedKmer);
                    
                    //now assemble out from it
                    bridgePoints = shortAssemble(rle_p, seedKmer, thresh, BRANCH_LIMIT, maxBranchLength);
                    
                    //remember to rev comp this also
                    vector<uint8_t> orig = string_util::reverseComplement_i(vector<uint8_t>(seq_i.begin(), seq_i.begin()+x+kmerSize));
                    
                    vector<pair<uint64_t, uint64_t> > edScores = vector<pair<uint64_t, uint64_t> >(bridgePoints.size());
                    uint64_t minScore = 0xFFFFFFFFFFFFFFFF;
                    for(uint64_t y = 0; y < bridgePoints.size(); y++) {
                        edScores[y] = editDistance_minimize(orig, bridgePoints[y]);
                        if(edScores[y].first < minScore) minScore = edScores[y].first;
                    }
                    
                    bridgePoints_ed.clear();
                    for(uint64_t y = 0; y < bridgePoints.size(); y++) {
                        if(edScores[y].first == minScore) bridgePoints_ed.push_back(vector<uint8_t>(bridgePoints[y].begin(), bridgePoints[y].begin()+edScores[y].second));
                    }
                    
                    if(bridgePoints_ed.size() == 0) {
                        //do nothing, we didn't find anything good
                    }
                    else if(bridgePoints_ed.size() == 1)
                    {
                        //one bridge with smallest edit distance
                        newCorr.start = 0;
                        newCorr.end = x+kmerSize;
                        newCorr.seq = string_util::reverseComplement_i(bridgePoints_ed[0]);
                        correctionsList.push_back(newCorr);
                    }
                    else {
                        //multiple with same edit distance, look at overall counts
                        uint64_t maxCount = 0;
                        uint64_t maxID = 0;
                        vector<uint64_t> edPU;
                        uint64_t summation;
                        for(uint64_t y = 0; y < bridgePoints_ed.size(); y++) {
                            edPU = rle_p->countPileup_i(bridgePoints_ed[y], kmerSize);
                            summation = 0;
                            summation = accumulate(edPU.begin(), edPU.end(), summation);
                            if(summation > maxCount) {
                                maxCount = summation;
                                maxID = y;
                            }
                        }
                        
                        //now save it
                        newCorr.start = 0;
                        newCorr.end = x+kmerSize;
                        newCorr.seq = string_util::reverseComplement_i(bridgePoints_ed[maxID]);
                        correctionsList.push_back(newCorr);
                    }
                }
            }
            else if(prevFound >= 0 && x < puSize) {
                //handle a bridging case
                /*
                for(uint64_t y = 0; y < kmerSize; y++) {
                    seedKmer[y] = seq_i[prevFound+y];
                    targetKmer[y] = seq_i[x+y];
                }
                */
                seedKmer.assign(seq_i.begin()+prevFound, seq_i.begin()+prevFound+kmerSize);
                targetKmer.assign(seq_i.begin()+x, seq_i.begin()+x+kmerSize);
                maxBranchLength = (uint64_t)(myParams.BRANCH_BUFFER_FACTOR*(x-prevFound+kmerSize));
                
                //printf("testing %d %d\n", maxBranchLength, myParams.MAX_BRANCH_ATTEMPT_LENGTH);
                
                if(maxBranchLength < myParams.MAX_BRANCH_ATTEMPT_LENGTH) {
                    //try forward first
                    bridgePoints = multiBridge(rle_p, seedKmer, targetKmer, thresh, BRANCH_LIMIT, maxBranchLength);
                    
                    //try reverse complement if we failed
                    if(bridgePoints.size() == 0) {
                        bridgePoints = multiBridge(rle_p, string_util::reverseComplement_i(targetKmer), string_util::reverseComplement_i(seedKmer), thresh, BRANCH_LIMIT, maxBranchLength);
                        
                        //make sure to fix the results here
                        for(unsigned int y = 0; y < bridgePoints.size(); y++) {
                            bridgePoints[y] = string_util::reverseComplement_i(bridgePoints[y]);
                        }
                    }
                    
                    //printf("bp size: %d\n", bridgePoints.size());
                    if(bridgePoints.size() == 0) {
                        //no bridges found
                        //calculate a midpoint
                        uint64_t midPoint = (uint64_t)((prevFound+x+kmerSize)/2.0);
                        maxBranchLength = (uint64_t)(myParams.TAIL_BUFFER_FACTOR*(midPoint-prevFound));
                        
                        if(maxBranchLength < myParams.MAX_BRANCH_ATTEMPT_LENGTH) {
                            //try to extend from the left to the middle
                            bridgePoints = shortAssemble(rle_p, seedKmer, thresh, BRANCH_LIMIT, maxBranchLength);
                            vector<uint8_t> orig = vector<uint8_t>(seq_i.begin()+prevFound, seq_i.begin()+midPoint);
                            
                            //calculate the minimized edit distances
                            vector<pair<uint64_t, uint64_t> > edScores = vector<pair<uint64_t, uint64_t> >(bridgePoints.size());
                            uint64_t minScore = 0xFFFFFFFFFFFFFFFF;
                            for(uint64_t y = 0; y < bridgePoints.size(); y++) {
                                edScores[y] = editDistance_minimize(orig, bridgePoints[y]);
                                if(edScores[y].first < minScore) minScore = edScores[y].first;
                            }
                            
                            //clip the strings by the minimized length
                            bridgePoints_ed.clear();
                            for(uint64_t y = 0; y < bridgePoints.size(); y++) {
                                if(edScores[y].first == minScore) bridgePoints_ed.push_back(vector<uint8_t>(bridgePoints[y].begin(), bridgePoints[y].begin()+edScores[y].second));
                            }
                            
                            //if(bridgePoints_ed.size() > 0 && minScore > (midPoint-prevFound)*.4) printf("big ED: %llu %llu\n", minScore, midPoint-prevFound);
                            
                            //TODO: make .4 a constant
                            if(bridgePoints_ed.size() == 0 || minScore > (midPoint-prevFound)*.4) {
                                //do nothing, we didn't find anything good
                            }
                            else if(bridgePoints_ed.size() == 1)
                            {
                                //one bridge with smallest edit distance
                                newCorr.start = prevFound;
                                newCorr.end = midPoint;
                                newCorr.seq = bridgePoints_ed[0];
                                correctionsList.push_back(newCorr);
                                //printf("left to mid found\n");
                            }
                            else {
                                //multiple with same edit distance, look at overall counts
                                uint64_t maxCount = 0;
                                uint64_t maxID = 0;
                                vector<uint64_t> edPU;
                                uint64_t summation;
                                for(uint64_t y = 0; y < bridgePoints_ed.size(); y++) {
                                    edPU = rle_p->countPileup_i(bridgePoints_ed[y], kmerSize);
                                    summation = 0;
                                    summation = accumulate(edPU.begin(), edPU.end(), summation);
                                    if(summation > maxCount) {
                                        maxCount = summation;
                                        maxID = y;
                                    }
                                }
                                
                                //now save it
                                newCorr.start = prevFound;
                                newCorr.end = midPoint;
                                newCorr.seq = bridgePoints_ed[maxID];
                                correctionsList.push_back(newCorr);
                                //printf("left to mid found\n");
                            }
                            
                            //try to extend from the right to the middle
                            vector<uint8_t> revTarget = string_util::reverseComplement_i(targetKmer);
                            
                            //now assemble out from it
                            bridgePoints = shortAssemble(rle_p, revTarget, thresh, BRANCH_LIMIT, maxBranchLength);
                            
                            //remember to rev comp this also
                            orig = string_util::reverseComplement_i(vector<uint8_t>(seq_i.begin()+midPoint, seq_i.begin()+x+kmerSize));
                            
                            edScores = vector<pair<uint64_t, uint64_t> >(bridgePoints.size());
                            minScore = 0xFFFFFFFFFFFFFFFF;
                            for(uint64_t y = 0; y < bridgePoints.size(); y++) {
                                edScores[y] = editDistance_minimize(orig, bridgePoints[y]);
                                if(edScores[y].first < minScore) minScore = edScores[y].first;
                            }
                            
                            bridgePoints_ed.clear();
                            for(uint64_t y = 0; y < bridgePoints.size(); y++) {
                                if(edScores[y].first == minScore) bridgePoints_ed.push_back(vector<uint8_t>(bridgePoints[y].begin(), bridgePoints[y].begin()+edScores[y].second));
                            }
                            
                            //TODO: make .4 a constant
                            if(bridgePoints_ed.size() == 0 || minScore > (midPoint-prevFound)*.4) {
                                //do nothing, we didn't find anything good
                            }
                            else if(bridgePoints_ed.size() == 1)
                            {
                                //one bridge with smallest edit distance
                                newCorr.start = midPoint;
                                newCorr.end = x+kmerSize;
                                newCorr.seq = string_util::reverseComplement_i(bridgePoints_ed[0]);
                                correctionsList.push_back(newCorr);
                                //printf("right to mid found\n");
                            }
                            else {
                                //multiple with same edit distance, look at overall counts
                                uint64_t maxCount = 0;
                                uint64_t maxID = 0;
                                vector<uint64_t> edPU;
                                uint64_t summation;
                                for(uint64_t y = 0; y < bridgePoints_ed.size(); y++) {
                                    edPU = rle_p->countPileup_i(bridgePoints_ed[y], kmerSize);
                                    summation = 0;
                                    summation = accumulate(edPU.begin(), edPU.end(), summation);
                                    if(summation > maxCount) {
                                        maxCount = summation;
                                        maxID = y;
                                    }
                                }
                                
                                //now save it
                                newCorr.start = midPoint;
                                newCorr.end = x+kmerSize;
                                newCorr.seq = string_util::reverseComplement_i(bridgePoints_ed[maxID]);
                                correctionsList.push_back(newCorr);
                                //printf("right to mid found\n");
                            }
                        }
                    }
                    else if(bridgePoints.size() == 1) {
                        //one bridge found, add it on
                        newCorr.start = prevFound;
                        newCorr.end = x+kmerSize;
                        newCorr.seq = bridgePoints[0];
                        correctionsList.push_back(newCorr);
                    }
                    else {
                        //multiple bridges found, pick the best one by edit distance
                        vector<uint8_t> orig = vector<uint8_t>(seq_i.begin()+prevFound, seq_i.begin()+x+kmerSize);
                        
                        vector<uint64_t> edScores = vector<uint64_t>(bridgePoints.size());
                        uint64_t minScore = 0xFFFFFFFFFFFFFFFF;
                        for(uint64_t y = 0; y < bridgePoints.size(); y++) {
                            edScores[y] = editDistance(orig, bridgePoints[y]);
                            if(edScores[y] < minScore) minScore = edScores[y];
                        }
                        
                        bridgePoints_ed.clear();
                        for(uint64_t y = 0; y < bridgePoints.size(); y++) {
                            if(edScores[y] == minScore) bridgePoints_ed.push_back(bridgePoints[y]);
                        }
                        if(bridgePoints_ed.size() == 1)
                        {
                            //one bridge with smallest edit distance
                            newCorr.start = prevFound;
                            newCorr.end = x+kmerSize;
                            newCorr.seq = bridgePoints_ed[0];
                            correctionsList.push_back(newCorr);
                        }
                        else {
                            //multiple with same edit distance, look at overall counts
                            uint64_t maxCount = 0;
                            uint64_t maxID = 0;
                            vector<uint64_t> edPU;
                            uint64_t summation;
                            for(uint64_t y = 0; y < bridgePoints_ed.size(); y++) {
                                edPU = rle_p->countPileup_i(bridgePoints_ed[y], kmerSize);
                                summation = 0;
                                summation = accumulate(edPU.begin(), edPU.end(), summation);
                                if(summation > maxCount) {
                                    maxCount = summation;
                                    maxID = y;
                                }
                            }
                            
                            //now save it
                            newCorr.start = prevFound;
                            newCorr.end = x+kmerSize;
                            newCorr.seq = bridgePoints_ed[maxID];
                            correctionsList.push_back(newCorr);
                        }
                    }
                    
                    //if we found a bridge, no need to keep trying
                    //if(bridgePoints.size() > 0) break;
                }
            }
        }
        else {
            //the counts were okay, no correction needed here
            x++;
        }
    }
    
    x = puSize;
    
    //use the tail factor for the buffer
    maxBranchLength = (uint64_t)(myParams.TAIL_BUFFER_FACTOR*(x-prevFound+kmerSize));
    if(maxBranchLength <= myParams.MAX_BRANCH_ATTEMPT_LENGTH && pu[puSize-1] < thresh && prevFound >= 0) {
        //copy the seed k-mer
        seedKmer.assign(seq_i.begin()+prevFound, seq_i.begin()+prevFound+kmerSize);
        bridgePoints = shortAssemble(rle_p, seedKmer, thresh, BRANCH_LIMIT, maxBranchLength);
        
        vector<uint8_t> orig = vector<uint8_t>(seq_i.begin()+prevFound, seq_i.end());
        
        //calculate the minimized edit distances
        vector<pair<uint64_t, uint64_t> > edScores = vector<pair<uint64_t, uint64_t> >(bridgePoints.size());
        uint64_t minScore = 0xFFFFFFFFFFFFFFFF;
        for(uint64_t y = 0; y < bridgePoints.size(); y++) {
            edScores[y] = editDistance_minimize(orig, bridgePoints[y]);
            if(edScores[y].first < minScore) minScore = edScores[y].first;
        }
        
        //clip the strings by the minimized length
        bridgePoints_ed.clear();
        for(uint64_t y = 0; y < bridgePoints.size(); y++) {
            if(edScores[y].first == minScore) bridgePoints_ed.push_back(vector<uint8_t>(bridgePoints[y].begin(), bridgePoints[y].begin()+edScores[y].second));
        }
        
        if(bridgePoints_ed.size() == 0) {
            //do nothing, we didn't find anything good
        }
        else if(bridgePoints_ed.size() == 1)
        {
            //one bridge with smallest edit distance
            newCorr.start = prevFound;
            newCorr.end = seq_i.size();
            newCorr.seq = bridgePoints_ed[0];
            correctionsList.push_back(newCorr);
        }
        else {
            //multiple with same edit distance, look at overall counts
            uint64_t maxCount = 0;
            uint64_t maxID = 0;
            vector<uint64_t> edPU;
            uint64_t summation;
            for(uint64_t y = 0; y < bridgePoints_ed.size(); y++) {
                edPU = rle_p->countPileup_i(bridgePoints_ed[y], kmerSize);
                summation = 0;
                summation = accumulate(edPU.begin(), edPU.end(), summation);
                if(summation > maxCount) {
                    maxCount = summation;
                    maxID = y;
                }
            }
            
            //now save it
            newCorr.start = prevFound;
            newCorr.end = seq_i.size();
            newCorr.seq = bridgePoints_ed[maxID];
            correctionsList.push_back(newCorr);
        }
    }
    
    //go through and insert the corrections in reverse
    vector<uint8_t> ret = vector<uint8_t>(seq_i.begin(), seq_i.end());
    Correction c;
    for(int64_t x = correctionsList.size()-1; x >=0; x--) {
        //get the modification
        c = correctionsList[x];
        
        //delete what was in the range before
        ret.erase(ret.begin()+c.start, ret.begin()+c.end);
        
        //insert the new values
        ret.insert(ret.begin()+c.start, c.seq.begin(), c.seq.end());
    }
    
    return ret;
}

CorrectionResults correctRead_job(int id, BaseBWT * rle_p, LongReadFA inputRead, Parameters myParams) {
    //prep the return value
    CorrectionResults ret;
    ret.label = inputRead.label;
    ret.originalSeq = inputRead.seq;
    
    //1 - translate string to vector<uint64_t>
    vector<uint8_t> seq_i = string_util::stoi(inputRead.seq);
    
    //2 - correct with small k
    vector<uint8_t> corrected_k = correctionPass(rle_p, seq_i, myParams, myParams.k);
    /*
    while(seq_i.size() != corrected_k.size() || !equal(seq_i.begin(), seq_i.end(), corrected_k.begin())) {
        //printf("looping here\n");
        seq_i = vector<uint8_t>(corrected_k);
        corrected_k = correctionPass(rle_p, seq_i, myParams, myParams.k);
    }
    */
    
    if(myParams.k == myParams.K) {
        //3a - k = K, skip second pass
        //4 - translate vector<uint64_t> to string
        ret.correctedSeq = string_util::itos(corrected_k);
        
        if(myParams.VERBOSE) {
            //seq_i = string_util::stoi(inputRead.seq);
            vector<uint64_t> c1 = rle_p->countPileup_i(seq_i, myParams.k);
            vector<uint64_t> c2 = rle_p->countPileup_i(corrected_k, myParams.k);
            ret.avgBefore = accumulate(c1.begin(), c1.end(), 0.0)/c1.size();
            ret.avgAfter = accumulate(c2.begin(), c2.end(), 0.0)/c2.size();
        }
        else {
            ret.avgBefore = 0;
            ret.avgAfter = 0;
        }
    }
    else {
        //3b - correct with big K
        vector<uint8_t> corrected_K = correctionPass(rle_p, corrected_k, myParams, myParams.K);
        
        //4 - translate vector<uint64_t> to string
        ret.correctedSeq = string_util::itos(corrected_K);
        /*
        while(seq_i.size() != corrected_K.size() || !equal(seq_i.begin(), seq_i.end(), corrected_K.begin())) {
            //printf("looping here 2\n");
            seq_i = vector<uint8_t>(corrected_K);
            corrected_K = correctionPass(rle_p, seq_i, myParams, myParams.K);
        }
        */
        
        if(myParams.VERBOSE) {
            //seq_i = string_util::stoi(inputRead.seq);
            vector<uint64_t> c1 = rle_p->countPileup_i(seq_i, myParams.k);
            vector<uint64_t> c2 = rle_p->countPileup_i(corrected_K, myParams.k);
            ret.avgBefore = accumulate(c1.begin(), c1.end(), 0.0)/c1.size();
            ret.avgAfter = accumulate(c2.begin(), c2.end(), 0.0)/c2.size();
        }
        else {
            ret.avgBefore = 0;
            ret.avgAfter = 0;
        }
    }
    
    //for(int x = 0; x < c2.size(); x++) printf("%llu ", c2[x]);
    //printf("\n");
    
    //5 - return result
    return ret;
}

int main(int argc, char* argv[]) {
    
    //////////////////////////////////////////////////////////
    //DEFAULT PARAMETERS
    Parameters myParams;
    myParams.USE_FM_INDEX = false;              //if enabled, we will use the RLE_BWT, else CSA_BWT
    myParams.k = 21;                            //small k-mer
    myParams.K = 59;                            //big K-mer
    myParams.MIN_COUNT = 5;                     //threshold for counting, overrides FRAC*<median of read counts>
    myParams.FRAC = 0.1;                        //the factor applied to the median to determine a dynamic threshold
    myParams.MAX_BRANCH_ATTEMPT_LENGTH = 10000; //maximum length of a gap that we will try to cross; longer can mean more CPU usage
    myParams.BRANCH_BUFFER_FACTOR = 1.3;        //the factor applied to any bridge gap to allow for insertions
    myParams.TAIL_BUFFER_FACTOR = 1.05;         //the factor applied to any head/tail gap to allow for insertions
    //myParams.MAX_TRIES = 1;                   //DEPRECATED
    myParams.FM_BIT_POWER = 8;                  //the in-memory FM-index samples at 2^FM_BIT_POWER; smaller = faster access but more memory
    myParams.VERBOSE = false;                   //if true, information about each read is dumped in the output
    
    uint64_t poolSize = 10000;                  //the number of jobs waiting to be processed at any given point in time; smaller may lead to lower process utilization but also less memory
    int numThreads = 1;                         //the number of concurrent correction threads
    uint64_t beginID = 0;                       //0-indexed id of the first read to process (default: beginning)
    uint64_t endID = 0xFFFFFFFFFFFFFFFF;        //0-indexed id of the last read to process (default: all reads)
    //////////////////////////////////////////////////////////
    
    char opt;
    bool helpRequest = false;
    while((opt = getopt(argc, argv, "hvk:K:p:b:e:m:f:iF:V")) != -1) {
        if(opt == 'h') helpRequest = true;
        else if(opt == 'v') {
            printf("fmlrc version %s\n", VERSION.c_str());
            return 0;
        }
        else if(opt == 'k') myParams.k = atoi(optarg);
        else if(opt == 'K') myParams.K = atoi(optarg);
        else if(opt == 'p') numThreads = atoi(optarg);
        else if(opt == 'b') beginID = atoi(optarg);
        else if(opt == 'e') endID = atoi(optarg);
        else if(opt == 'm') myParams.MIN_COUNT = atoi(optarg);
        else if(opt == 'f') myParams.FRAC = atof(optarg);
        else if(opt == 'V') myParams.VERBOSE = true;
        //MAX_BRANCH_ATTEMPT_LENGTH
        //BRANCH_BUFFER_FACTOR
        //TAIL_BUFFER_FACTOR
        else if(opt == 'i') myParams.USE_FM_INDEX = true;
        else if(opt == 'F') myParams.FM_BIT_POWER = atoi(optarg);
        else printf("UNHANDLED OPTION: %d %c %s\n", optind, opt, optarg);
    }
    
    if(argc-optind < 3 || helpRequest) {
        printf("Usage:   fmlrc [options] <comp_msbwt.npy> <long_reads.fa> <corrected_reads.fa>\n");
        printf("Options: -h        print help menu\n");
        printf("         -v        print version number and exit\n");
        printf("         -k INT    small k-mer size (default: 21)\n");
        printf("         -K INT    large K-mer size (default: 59), set K=k for single pass\n");
        printf("         -p INT    number of correction threads\n");
        printf("         -b INT    index of read to start with (default: 0)\n");
        printf("         -e INT    index of read to end with (default: end of file)\n");
        printf("         -m INT    absolute minimum count to consider a path (default: 5)\n");
        printf("         -f FLOAT  dynamic minimum fraction of median to consider a path (default: .10)\n");
        printf("         -i        build a sampled FM-index instead of bit arrays\n");
        printf("         -F INT    FM-index is sampled every 2**<-F> values (default: 8); requires -i\n");
        printf("         -V        verbose output\n");
        return 0;
    }
    
    if(beginID > endID) {
        printf("ERROR: parameter -b must be less than or equal to parameter -e\n");
        return 1;
    }
    
    if(myParams.FRAC < 0  || myParams.FRAC > 1) {
        printf("ERROR: parameter -f must be within the range [0, 1]\n");
        return 1;
    }
    
    //we need to always call this once
    string_util::initializeStringUtil();
    
    //load the BWT into memory
    char * bwtFN = argv[optind];
    BaseBWT * rle;// = new CSA_BWT(bwtFN, myParams.FM_BIT_POWER);
    
    if(myParams.USE_FM_INDEX) rle = new RLE_BWT(bwtFN, myParams.FM_BIT_POWER);
    else rle = new CSA_BWT(bwtFN, false); //THE false MEANS WE PROMISE NOT TO QUERY '$'
    
    //open the fasta file for reading
    char * longReadFN = argv[optind+1];
    FastaIterator fi(longReadFN);
    
    //open the output fasta file for writing
    char * correctedReadFN = argv[optind+2];
    FastaWriter fw(correctedReadFN);
    
    //now we need to set up our pool and stuff
    ctpl::thread_pool myPool(numThreads);
    
    //skip however many reads we were told to skip
    if(beginID > 0) printf("Skipping %llu reads...\n", beginID);
    uint64_t skippedReadCount = 0;
    while(skippedReadCount < beginID && fi.isMore()) {
        fi.getNextRead();
        skippedReadCount++;
    }
    
    uint64_t jobsToProcess = endID - beginID;
    uint64_t jobsLoaded = 0;
    uint64_t jobsCompleted = 0;
    
    //preload the first <poolSize> jobs
    vector<std::future<CorrectionResults> > results(poolSize);
    LongReadFA inputRead;
    for(uint64_t x = 0; x < poolSize && fi.isMore() && jobsLoaded < jobsToProcess; x++) {
        inputRead = fi.getNextRead();
        results[x] = myPool.push(correctRead_job, rle, inputRead, myParams);
        jobsLoaded++;
    }
    
    //now load the jobs as they empty out
    uint64_t currJobSlot = 0;
    CorrectionResults currResults;
    LongReadFA outRead;
    while(fi.isMore() && jobsLoaded < jobsToProcess) {
        //get the results
        currResults = results[currJobSlot].get();
        outRead.label = currResults.label;
        outRead.seq = currResults.correctedSeq;
        fw.writeRead(outRead);
        if(myParams.VERBOSE) printf("%llu: avg change %lf -> %lf\n", beginID+jobsCompleted, currResults.avgBefore, currResults.avgAfter);
        jobsCompleted++;
        
        //load the next job in
        inputRead = fi.getNextRead();
        results[currJobSlot] = myPool.push(correctRead_job, rle, inputRead, myParams);
        jobsLoaded++;
        
        //increment the slot we are looking at, looping around if necessary
        currJobSlot++;
        if(currJobSlot == poolSize){
            currJobSlot = 0;
            if(!myParams.VERBOSE) printf("Processed %llu reads\n", jobsCompleted);
        }
    }
    
    //now we just wait on the remaining jobs to finish
    while(jobsCompleted < jobsLoaded) {
        //get the results
        currResults = results[currJobSlot].get();
        outRead.label = currResults.label;
        outRead.seq = currResults.correctedSeq;
        fw.writeRead(outRead);
        if(myParams.VERBOSE) printf("%llu: avg change %lf -> %lf\n", beginID+jobsCompleted, currResults.avgBefore, currResults.avgAfter);
        jobsCompleted++;
        
        //increment the slot we are looking at, looping around if necessary
        currJobSlot++;
        if(currJobSlot == poolSize){
            currJobSlot = 0;
            if(!myParams.VERBOSE) printf("Processed %llu reads\n", jobsCompleted);
        }
    }
    
    //this is the only thing to clean up
    delete rle;
    printf("Finished processing reads [%llu, %llu)\n", beginID, beginID+jobsCompleted);
    
    return 0;
}