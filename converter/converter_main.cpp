
//C headers
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>

//C++ headers
#include <algorithm>
#include <string>
#include <sstream>

struct Parameters {
    bool FORCE_OVERWRITE;
    bool USE_STDIN;
    std::string filename;
};

const std::string VERSION = "1.0.0";

int runConverter(Parameters myParams, char * outFN) {
    FILE * inputStream;
    if(myParams.USE_STDIN) {
        inputStream = stdin;
    }
    else {
        inputStream = fopen(myParams.filename.c_str(), "r");
    }

    FILE * outputStream = fopen(outFN, "w+");

    //TODO: lots of header related items
    unsigned long BUFFER_SIZE = 1024;
    unsigned char buffer[BUFFER_SIZE]; 

    //most of the files I've seen are 80 and '\x46', I'm increasing it just in case
    unsigned long headerSize = 96;
    std::string headerHex = "\x56";
    unsigned long x;
    for(x=0; x < headerSize-1; x++) {
        buffer[x] = 32; //hex value 20 = ' '
    }
    buffer[headerSize-1] = 10; //hex value 0a = '\n'
    fwrite(buffer, 1, headerSize, outputStream);
    
    //set up the translation, default is 255
    unsigned char translator[256];
    std::fill_n(translator, 256, 255);

    std::string validSymbols = "$ACGNT";
    x = 0;
    for(char &c: validSymbols) {
        translator[(unsigned char)c] = x;
        x++;
    }

    //first read
    unsigned long readBytes = fread(buffer, 1, BUFFER_SIZE, inputStream);

    unsigned char currSym = buffer[0];
    unsigned long currCount = 0;
    unsigned char writeByte;
    unsigned long bytesWritten = 0;

    //core loop
    while(readBytes > 0) {
        for(x = 0; x < readBytes; x++) {
            if(currSym == buffer[x]) {
                currCount++;
            }
            else {
                //check for invalid character; if it's the new line symbol, we will ignore it
                if(translator[currSym] == 255) {
                    if(currSym != 10){
                        printf("UNEXPECTED SYMBOL DETECTED: char: \"%c\", hex: \"%x\"\n", currSym, currSym);
                        return 1;
                    }
                }
                else {
                    //we are at the end of the run so handle it
                    while(currCount > 0) {
                        writeByte = translator[currSym] | ((currCount & 0x1F) << 3);
                        fwrite(&writeByte, 1, 1, outputStream);
                        currCount = currCount >> 5;
                        bytesWritten += 1;
                    }
                        
                    //get the next symbol since it's valid and start a new run
                    currSym = buffer[x];
                    currCount = 1;
                }
            }
        }
        
        //get the next batch to parse through
        readBytes = fread(buffer, 1, BUFFER_SIZE, inputStream);
    }

    //handle the last run
    if(translator[currSym] == 255){
        if(currSym != 10){
            printf("UNEXPECTED SYMBOL DETECTED: char: \"%c\", hex: \"%x\"\n", currSym, currSym);
            return 1;
        }
    }
    else {
        //we are at the end of the last run so handle it
        while(currCount > 0) {
            writeByte = translator[currSym] | ((currCount & 0x1F) << 3);
            fwrite(&writeByte, 1, 1, outputStream);
            currCount = currCount >> 5;
            bytesWritten += 1;
        }
            
        //clear these
        currSym = 0;
        currCount = 0;
    }

    //we have finished the compression part, close the input file
    fclose(inputStream);
    
    //now that we know the total length, fill in the bytes for our header
    //char initialWrite[headerSize];
    //sprintf(initialWrite, "\x93NUMPY\x01\x00%s\x00{\'descr\': \'|u1\', \'fortran_order\': False, \'shape\': (%d,), }", headerHex.c_str(), bytesWritten);
    //std::ostringstream stringStream;
    //stringStream << "\x93NUMPY\x01\x00" << headerHex << "\x00{\'descr\': \'|u1\', \'fortran_order\': False, \'shape\': (" << bytesWritten << ",), }";
    //std::string initialWrite = stringStream.str();

    //have to do some special things due to the \x00 characters; might be a better way to do this
    std::string initialWrite = "\x93NUMPY\x01";
    initialWrite.push_back('\0');
    initialWrite.push_back(headerHex.c_str()[0]);
    initialWrite.push_back('\0');
    initialWrite += "{\'descr\': \'|u1\', \'fortran_order\': False, \'shape\': (";
    initialWrite += std::to_string(bytesWritten);
    initialWrite += ",), }";

    fseek(outputStream, 0, SEEK_SET);
    fwrite(initialWrite.c_str(), 1, initialWrite.length(), outputStream);
    printf("init write len: %lu\n", initialWrite.length());
    //finally close it all out
    fclose(outputStream);
    
    return 0;
}

int main(int argc, char* argv[]) {
    
    //////////////////////////////////////////////////////////
    //DEFAULT PARAMETERS
    Parameters myParams;
    myParams.FORCE_OVERWRITE = false;       //if true, this will silently overwrite the output file if it exists
    myParams.USE_STDIN = true;              //if true, we will read all input from stdin
    myParams.filename = "";                 //if the previous parameter is false, this is the filename we are reading from
    //////////////////////////////////////////////////////////
    
    char opt;
    bool helpRequest = false;
    while((opt = getopt(argc, argv, "hvfi:")) != -1) {
        if(opt == 'h') helpRequest = true;
        else if(opt == 'v') {
            printf("fmlrc-convert version %s\n", VERSION.c_str());
            return 0;
        }
        else if(opt == 'i') {
            myParams.USE_STDIN = false;
            myParams.filename = optarg;
        }else if (opt == 'f') myParams.FORCE_OVERWRITE = true;
        else printf("UNHANDLED OPTION: %d %c %s\n", optind, opt, optarg);
    }

    if(argc-optind < 1 || helpRequest) {
        printf("Usage:   fmlrc-convert [options] <out_comp_mbswt.npy>\n");
        printf("Options: -h        print help menu\n");
        printf("         -v        print version number and exit\n");
        printf("         -f        force overwrite of existing file (default: false)\n");
        printf("         -i STR    the plain text BWT file to be converted into msbwt format (default: stdin)\n");
        return 0;
    }

    //Input error checking
    char * bwtFN = argv[optind];
    struct stat buffer;
    if(!myParams.FORCE_OVERWRITE && stat(bwtFN, &buffer) == 0) {
        printf("ERROR: output file already exists, use -f to force overwrite\n");
        return 1;
    }

    if(!myParams.USE_STDIN && stat(myParams.filename.c_str(), &buffer) != 0) {
        printf("ERROR: input filename does not exist\n");
        return 1;
    }

    //lets do the things
    return runConverter(myParams, bwtFN);
}