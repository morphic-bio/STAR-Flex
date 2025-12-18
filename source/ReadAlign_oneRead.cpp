#include "ReadAlign.h"
#include "readLoad.h"
#include "readBarcodeLoad.h"
#include "SequenceFuns.h"
#include "ErrorWarning.h"
#include "GlobalVariables.h"
#include "libtrim/trim.h"

int ReadAlign::oneRead() {//process one read: load, map, write

    //load read name, sequence, quality from the streams into internal arrays
    int readStatus[P.readNends];

    for (uint32 im=0; im<P.readNends; im++) {
        readStatus[im] = readLoad(*(readInStream[im]), P, readLength[im], readLengthOriginal[im], readNameMates[im], Read0[im], Read1[im], Qual0[im], clipMates[im], iReadAll, readFilesIndex, readFilter, readNameExtra[im]);
        if (readStatus[im] != readStatus[0]) {//check if the end of file was reached or not for all files
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR: read files are not consistent, reached the end of the one before the other one\n";
            errOut << "SOLUTION: Check you your input files: they may be corrupted\n";
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        };
    };

    if (readStatus[0]==-1) {//finished with the stream
        return -1;
    };    
    
    // Increment read counters BEFORE trimming (so dropped reads are counted)
    statsRA.readN++;
    statsRA.readBases += readLength[0] + (P.readNmates == 2 ? readLength[1] : 0);
    
    // Cutadapt-style trimming (if enabled)
    if (P.trimCutadapt == "Yes") {
        if (P.readNmates == 2) {
        struct TrimParams params;
        trim_params_init(&params);
        params.quality_cutoff = P.trimCutadaptQuality;
        params.min_length = P.trimCutadaptMinLength;
        // Set adapters if custom (validate exactly 2 adapters)
        if (P.trimCutadaptAdapter.size() == 2 && P.trimCutadaptAdapter[0] != "-" && P.trimCutadaptAdapter[1] != "-") {
            params.adapter_r1 = P.trimCutadaptAdapter[0].c_str();
            params.adapter_r2 = P.trimCutadaptAdapter[1].c_str();
        } else if (P.trimCutadaptAdapter.size() > 0 && (P.trimCutadaptAdapter[0] != "-" || (P.trimCutadaptAdapter.size() > 1 && P.trimCutadaptAdapter[1] != "-"))) {
            // Invalid adapter specification
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR: --trimCutadaptAdapter requires exactly 2 adapter sequences (R1 and R2), separated by space\n";
            errOut << "Provided: " << P.trimCutadaptAdapter.size() << " adapter(s)\n";
            errOut << "SOLUTION: Provide both R1 and R2 adapters, or use '-' to use default TruSeq adapters\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        }
        
        struct TrimResult result1, result2;
        uint32_t len1 = (uint32_t)readLength[0];
        uint32_t len2 = (uint32_t)readLength[1];
        trim_pair(Read0[0], Qual0[0], &len1,
                  Read0[1], Qual0[1], &len2,
                  &params, &result1, &result2);
        readLength[0] = len1;
        readLength[1] = len2;
        
        // Accumulate stats using centralized helper (counting reads, not pairs)
        struct TrimStats trimStats = {0, 0, 0, 0, 0};
        trim_stats_add(&trimStats, &result1);  // R1
        trim_stats_add(&trimStats, &result2);  // R2
        
        // Copy accumulated stats to STAR's Stats class
        statsRA.trimReadsProcessed += trimStats.reads_processed;
        statsRA.trimReadsTrimmed += trimStats.reads_trimmed;
        statsRA.trimReadsTooShort += trimStats.reads_too_short;
        statsRA.trimBasesQualityTrimmed += trimStats.bases_quality_trimmed;
        statsRA.trimBasesAdapterTrimmed += trimStats.bases_adapter_trimmed;
        
        // Handle dropped reads
        if (result1.dropped || result2.dropped) {
            readFilter = 'Y';  // Fail QC - same as other filters
            // Skip mapping for this pair (readN already incremented above)
            return 0;
        }
        
        // Update original lengths after trimming
        readLengthOriginal[0] = readLength[0];
        readLengthOriginal[1] = readLength[1];
        
        // Re-convert to numeric after trimming
        convertNucleotidesToNumbers(Read0[0], Read1[0], readLength[0]);
        convertNucleotidesToNumbers(Read0[1], Read1[1], readLength[1]);
        } else if (P.readNmates == 1) {
            // Single-end trimming
            struct TrimParams params;
            trim_params_init(&params);
            params.quality_cutoff = P.trimCutadaptQuality;
            params.min_length = P.trimCutadaptMinLength;
            // Use R1 adapter for single-end
            const char* adapter = TRUSEQ_ADAPTER_R1;
            if (P.trimCutadaptAdapter.size() >= 1 && P.trimCutadaptAdapter[0] != "-") {
                adapter = P.trimCutadaptAdapter[0].c_str();
            }
            
            struct TrimResult result;
            uint32_t len1 = (uint32_t)readLength[0];
            result = trim_read(Read0[0], Qual0[0], len1, adapter, &params);
            readLength[0] = result.new_length;
            
            // Accumulate stats using centralized helper
            struct TrimStats trimStats = {0, 0, 0, 0, 0};
            trim_stats_add(&trimStats, &result);
            
            // Copy accumulated stats to STAR's Stats class
            statsRA.trimReadsProcessed += trimStats.reads_processed;
            statsRA.trimReadsTrimmed += trimStats.reads_trimmed;
            statsRA.trimReadsTooShort += trimStats.reads_too_short;
            statsRA.trimBasesQualityTrimmed += trimStats.bases_quality_trimmed;
            statsRA.trimBasesAdapterTrimmed += trimStats.bases_adapter_trimmed;
            
            // Handle dropped reads
            if (result.dropped) {
                readFilter = 'Y';
                return 0;
            }
            
            // Update length and re-convert
            readLengthOriginal[0] = readLength[0];
            convertNucleotidesToNumbers(Read0[0], Read1[0], readLength[0]);
        }
    }
    
    if (P.outFilterBySJoutStage != 2) {
        for (uint32 im=0; im<P.readNmates; im++) {//not readNends: the barcode quality will be calculated separately
            for (uint64 ix=clipMates[im][0].clippedN; ix<readLengthOriginal[im]-clipMates[im][1].clippedN; ix++) {
                qualHist[im][(uint8)Qual0[im][ix]]++;
            };
        };
    };
    
    if (P.readNmates==2) {//combine two mates together
        Lread=readLength[0]+readLength[1]+1;
        readLengthPairOriginal=readLengthOriginal[0]+readLengthOriginal[1]+1;
        if (Lread>DEF_readSeqLengthMax) {
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR in reads input: Lread of the pair = " << Lread << "   while DEF_readSeqLengthMax=" << DEF_readSeqLengthMax <<endl;
            errOut << "Read Name="<<readNameMates[0]<<endl;
            errOut << "SOLUTION: increase DEF_readSeqLengthMax in IncludeDefine.h and re-compile STAR"<<endl<<flush;
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        };

        //marker for spacer base
        Read1[0][readLength[0]]=MARK_FRAG_SPACER_BASE;
        
        //copy 2nd mate into Read1[0] & reverse-complement
        complementSeqNumbers(Read1[1],Read1[0]+readLength[0]+1,readLength[1]);//complement. Here Read1[1] is still the 2nd mate's numeric-sequence. Later Read1[1] will be reverse complement of the combined read.
        for (uint ii=0;ii<readLength[1]/2;ii++) {
            swap(Read1[0][Lread-ii-1],Read1[0][ii+readLength[0]+1]); //reverse
        };

    } else {//1 mate

        if (readStatus[0]==-1) {//finished with the stream
            return -1;
        };

        Lread=readLength[0];
        readLengthPairOriginal=readLengthOriginal[0];
        readLength[1]=0;

    };
      
    readFileType=readStatus[0];

    complementSeqNumbers(Read1[0],Read1[1],Lread); //returns complement of Reads[ii]
    for (uint ii=0;ii<Lread;ii++) {//reverse
        Read1[2][Lread-ii-1]=Read1[1][ii];
    };

    //max number of mismatches allowed for this read
    outFilterMismatchNmaxTotal=min(P.outFilterMismatchNmax, (uint) (P.outFilterMismatchNoverReadLmax*(readLength[0]+readLength[1])));

    //map the read
    if (P.pGe.gType==101) {//SpliceGraph
        mapOneReadSpliceGraph();
    } else {//all other cases - standard alignment algorithm
        mapOneRead();
    };

    peOverlapMergeMap();
    
    multMapSelect();
    
    mappedFilter();  
    
    transformGenome();//for now genome transformation happens after multimapper selection, and mapping filter

    if (!peOv.yes) {//if the alignment was not mates merged - otherwise the chimeric detection was already done
        chimericDetection();
    };

    if (P.pCh.out.bam && chimRecord) {//chimeric alignment was recorded in main BAM files, and it contains the representative portion, so non-chimeric aligmnent is not output
        return 0;
    };

    waspMap();

    #ifdef OFF_BEFORE_OUTPUT
        #warning OFF_BEFORE_OUTPUT
        return 0;
    #endif

    //write out alignments
    outputAlignments();

    {
    #ifdef DEBUG_OutputLastRead
        lastReadStream.seekp(ios::beg);
        lastReadStream << iReadAll <<" "<< readName <<endl;
    #endif
    };

    return 0;

};


