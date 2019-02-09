
//C headers
#include <stdio.h>
#include <unistd.h>

//C++ headers
#include <string>

struct Parameters {
    bool USE_STDIN;
    std::string filename;
};

const std::string VERSION = "1.0.0";

int main(int argc, char* argv[]) {
    
    //////////////////////////////////////////////////////////
    //DEFAULT PARAMETERS
    Parameters myParams;
    myParams.USE_STDIN = true;              //if true, we will read all input from stdin
    myParams.filename = "";                 //if the previous parameter is false, this is the filename we are reading from
    //////////////////////////////////////////////////////////
    
    char opt;
    bool helpRequest = false;
    while((opt = getopt(argc, argv, "hvf:")) != -1) {
        if(opt == 'h') helpRequest = true;
        else if(opt == 'v') {
            printf("fmlrc-convert version %s\n", VERSION.c_str());
            return 0;
        }
        else if(opt == 'f') {
            myParams.USE_STDIN = false;
            myParams.filename = optarg;
        }
        else printf("UNHANDLED OPTION: %d %c %s\n", optind, opt, optarg);
    }

    if(helpRequest) {
        printf("Usage:   fmlrc-convert [options]\n");
        printf("Options: -h        print help menu\n");
        printf("         -v        print version number and exit\n");
        printf("         -f STR    the plain text BWT file to be converted into msbwt format (default: stdin)\n");
        return 0;
    }

    if(myParams.USE_STDIN) {
        printf("Using stdin\n");
    }
    else {
        printf("filename: %s\n", myParams.filename.c_str());
    }
}