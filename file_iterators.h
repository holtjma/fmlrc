#ifndef FILE_ITERATORS_H
#define FILE_ITERATORS_H

//C headers
#include <stdint.h>

//C++ headers
#include <fstream>
#include <string>

using namespace std;

struct LongReadFA {
    string label;
    string seq;
};

class FastaIterator {
private:
    bool isMoreData;
    string fastaFN;
    string nextLine;
    ifstream ifp;
public:
    //constructor
    FastaIterator(string fastaFN);
    
    //funcs that matter
    inline bool isMore() { return this->isMoreData;};
    struct LongReadFA getNextRead();
};

class FastaWriter {
private:
    string fastaFN;
    int symsPerLine;
    ofstream ofp;
public:
    //constructor
    FastaWriter(string fastaFN, int symsPerLine=50);
    
    //destructor
    ~FastaWriter();
    
    //funcs that matter
    bool writeRead(LongReadFA r);
};

#endif