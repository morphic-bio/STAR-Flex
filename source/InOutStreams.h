#ifndef INOUTSTREAMS_DEF
#define INOUTSTREAMS_DEF

#include "IncludeDefine.h"
#include SAMTOOLS_BGZF_H

class InOutStreams {
    public:
    ostream *logStdOut, *outSAM;
    ofstream logStdOutFile, outSAMfile;
    BGZF *outBAMfileUnsorted, *outBAMfileCoord, *outQuantBAMfile;
    BGZF *outBAMfileY, *outBAMfileNoY;  // Y-chromosome split BAM handles
    
    //solo tmp file for two-pass unsorted CB/UB injection
    ofstream outBAMfileUnsortedSoloTmp;

    ofstream outChimSAM, outChimJunction, logMain, logProgress, logFinal, outUnmappedReadsStream[MAX_N_MATES];
    ifstream readIn[MAX_N_MATES];

    //compilation-optional streams
    ofstream outLocalChains;

    InOutStreams();
    ~InOutStreams();
};

#endif
