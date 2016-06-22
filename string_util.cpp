
//C headers
#include <stdint.h>

//C++ headers
#include <vector>

//Custom headers
#include "string_util.h"

namespace string_util {
    //this is externed in the header file, so declare it here
    vector<uint8_t> STRING_TO_INT;
    vector<uint8_t> INT_TO_STRING;
    vector<uint8_t> REV_COMP_I;
    
    void initializeStringUtil() {
        //set everything initially to 'N'
        STRING_TO_INT = vector<uint8_t>(256, 4);
        
        //now set the specific values
        STRING_TO_INT['$'] = 0;
        STRING_TO_INT['A'] = 1;
        STRING_TO_INT['C'] = 2;
        STRING_TO_INT['G'] = 3;
        STRING_TO_INT['N'] = 4;
        STRING_TO_INT['T'] = 5;
        
        //set the lower-case ones also
        STRING_TO_INT['a'] = 1;
        STRING_TO_INT['c'] = 2;
        STRING_TO_INT['g'] = 3;
        STRING_TO_INT['n'] = 4;
        STRING_TO_INT['t'] = 5;
        
        //do the reverse array also
        INT_TO_STRING = vector<uint8_t>(6);
        INT_TO_STRING[0] = '$';
        INT_TO_STRING[1] = 'A';
        INT_TO_STRING[2] = 'C';
        INT_TO_STRING[3] = 'G';
        INT_TO_STRING[4] = 'N';
        INT_TO_STRING[5] = 'T';
        
        //initialize the reverse-complement arrays
        REV_COMP_I = vector<uint8_t>(6);
        REV_COMP_I[0] = 0;//$$
        REV_COMP_I[1] = 5;//AT
        REV_COMP_I[2] = 3;//CG
        REV_COMP_I[3] = 2;//GC
        REV_COMP_I[4] = 4;//NN
        REV_COMP_I[5] = 1;//TA
    }
}