#include <sys/types.h>
#include <sys/stat.h>

#include "IncludeDefine.h"
#include "Parameters.h"
#include "SequenceFuns.h"
#include "Genome.h"
#include "Chain.h"
#include "ReadAlignChunk.h"
#include "ReadAlign.h"
#include "Stats.h"
#include "genomeGenerate.h"
#include "outputSJ.h"
#include "ThreadControl.h"
#include "GlobalVariables.h"
#include "TimeFunctions.h"
#include "ErrorWarning.h"
#include "sysRemoveDir.h"
#include "BAMfunctions.h"
#include "bamSortByCoordinate.h"
#include "Transcriptome.h"
#include "signalFromBAM.h"
#include "mapThreadsSpawn.h"
#include "SjdbClass.h"
#include "sjdbInsertJunctions.h"
#include "Variation.h"
#include "Solo.h"
#include "samHeaders.h"
#include "systemFunctions.h"
#include "ProbeListIndex.h"
#include "InlineCBCorrection.h"
#include <memory>
#include <cstdlib>
#include <cerrno>
#include <cstring>

#include "twoPassRunPass1.h"

#include "htslib/htslib/sam.h"

namespace {

bool ensureDirectoryTree(const std::string &path, mode_t mode, std::string &failedPath, int &failedErrno) {
    if (path.empty()) {
        return true;
    }
    std::string normalized = path;
    while (!normalized.empty() && normalized.back() == '/') {
        normalized.pop_back();
    }
    if (normalized.empty()) {
        return true;
    }

    size_t start = 0;
    std::string current;
    if (!normalized.empty() && normalized[0] == '/') {
        current = "/";
        start = 1;
    }

    while (start <= normalized.size()) {
        size_t slash = normalized.find('/', start);
        size_t len = (slash == std::string::npos) ? normalized.size() - start : slash - start;
        std::string token = normalized.substr(start, len);
        if (!token.empty()) {
            if (current.empty()) {
                current = token;
            } else if (current == "/") {
                current += token;
            } else {
                current += "/" + token;
            }
            if (mkdir(current.c_str(), mode) != 0 && errno != EEXIST) {
                failedPath = current;
                failedErrno = errno;
                return false;
            }
        }
        if (slash == std::string::npos) {
            break;
        }
        start = slash + 1;
    }
    return true;
}

} // namespace
#include "parametersDefault.xxd"

void usage(int usageType)
{
    cout << "Usage: STAR  [options]... --genomeDir /path/to/genome/index/   --readFilesIn R1.fq R2.fq\n";
    cout << "Spliced Transcripts Alignment to a Reference (c) Alexander Dobin, 2009-2022\n\n";
    cout << "STAR version=" << STAR_VERSION << "\n";
    cout << "STAR compilation time,server,dir=" << COMPILATION_TIME_PLACE << "\n";
    cout << "For more details see:\n";
    cout << "<https://github.com/alexdobin/STAR>\n";
    cout << "<https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf>\n";

    if (usageType == 0)
    { // brief
        cout << "\nTo list all parameters, run STAR --help\n";
    }
    else if (usageType == 1)
    { // full
        cout.write(reinterpret_cast<char *>(parametersDefault),
                   parametersDefault_len);
    };
    exit(0);
};

int main(int argInN, char *argIn[])
{
    // If no argument is given, or the first argument is either '-h' or '--help', run usage()
    if (argInN == 1)
    {
        usage(0);
    }
    else if (argInN == 2 && (strcmp("-h", argIn[1]) == 0 || strcmp("--help", argIn[1]) == 0))
    {
        usage(1);
    };

    time(&g_statsAll.timeStart);

    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////// Parameters
    Parameters P; // all parameters
    P.inputParameters(argInN, argIn);

    *(P.inOut->logStdOut) << "\t" << P.commandLine << '\n';
    *(P.inOut->logStdOut) << "\tSTAR version: " << STAR_VERSION << "   compiled: " << COMPILATION_TIME_PLACE << '\n';
    *(P.inOut->logStdOut) << timeMonthDayTime(g_statsAll.timeStart) << " ..... started STAR run\n"
                          << flush;

    // runMode
    if (P.runMode == "alignReads" || P.runMode == "soloCellFiltering")
    {
        // continue
    }
    else if (P.runMode == "genomeGenerate")
    {
        { // normal genome generation
            Genome genomeMain(P, P.pGe);
            genomeMain.genomeGenerate();
        };

        if (P.pGe.transform.type > 0)
        { // generate original genome, in addition to the transfomed generated above
            P.pGe.transform.type = 0;
            P.pGe.transform.typeString = "None";
            P.pGe.transform.vcfFile = "-";
            P.pGe.gDir += "/OriginalGenome/";
            Genome genomeOrig(P, P.pGe);
            genomeOrig.genomeGenerate();
        };

        sysRemoveDir(P.outFileTmp);
        P.inOut->logMain << "DONE: Genome generation, EXITING\n"
                         << flush;
        exit(0);
    }
    else if (P.runMode == "liftOver")
    {
        for (uint ii = 0; ii < P.pGe.gChainFiles.size(); ii++)
        {
            Chain chain(P, P.pGe.gChainFiles.at(ii));
            chain.liftOverGTF(P.pGe.sjdbGTFfile, P.outFileNamePrefix + "GTFliftOver_" + to_string(ii + 1) + ".gtf");
            P.inOut->logMain << "DONE: lift-over of GTF file, EXITING\n"
                             << flush;
            exit(0);
        };
    }
    else
    {
        P.inOut->logMain << "EXITING because of INPUT ERROR: unknown value of input parameter runMode=" << P.runMode << endl
                         << flush;
        exit(1);
    };

    // transcripome placeholder
    Transcriptome *transcriptomeMain = NULL;

    // this will execute --runMode soloCellFiltering and exit
    Solo soloCellFilter(P, *transcriptomeMain);

    ////////////////////////////////////////////////////////////////////////
    ///////////////////////////////// Genome
    Genome genomeMain(P, P.pGe);
    genomeMain.genomeLoad();

    if (P.pGe.transform.outYes) {
        genomeMain.Var = new Variation(P, genomeMain.chrStart, genomeMain.chrNameIndex, false);//no variation for mapGen, only for genOut
        genomeMain.genomeOut.g->Var = new Variation(P, genomeMain.genomeOut.g->chrStart, genomeMain.genomeOut.g->chrNameIndex, P.var.yes);
    } else {
        genomeMain.Var = new Variation(P, genomeMain.chrStart, genomeMain.chrNameIndex, P.var.yes);
    };

    SjdbClass sjdbLoci;

    if (P.sjdbInsert.pass1) {
        Genome genomeMain1 = genomeMain; // not sure if I need to create the copy - genomeMain1 below should not be changed
        sjdbInsertJunctions(P, genomeMain, genomeMain1, sjdbLoci);
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////START
    if (P.runThreadN > 1)
    {
        g_threadChunks.threadArray = new pthread_t[P.runThreadN];
        pthread_mutex_init(&g_threadChunks.mutexInRead, NULL);
        pthread_mutex_init(&g_threadChunks.mutexOutSAM, NULL);
        pthread_mutex_init(&g_threadChunks.mutexOutBAM1, NULL);
        pthread_mutex_init(&g_threadChunks.mutexOutUnmappedFastx, NULL);
        pthread_mutex_init(&g_threadChunks.mutexOutFilterBySJout, NULL);
        pthread_mutex_init(&g_threadChunks.mutexStats, NULL);
        pthread_mutex_init(&g_threadChunks.mutexBAMsortBins, NULL);
        pthread_mutex_init(&g_threadChunks.mutexError, NULL);
    };

    g_statsAll.progressReportHeader(P.inOut->logProgress);

    /////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////// 2-pass 1st pass
    twoPassRunPass1(P, genomeMain, transcriptomeMain, sjdbLoci);

    if (P.quant.yes)
    { // load transcriptome
        transcriptomeMain = new Transcriptome(P);
    };

    // Pre-initialize inline CB correction whitelist (needed before mapping)
    if (P.pSolo.inlineCBCorrection && P.pSolo.cbWLyes && !P.pSolo.cbWLstr.empty()) {
        InlineCBCorrection::initializeWhitelist(P.pSolo);
        P.inOut->logMain << "[INLINE-CB-INIT] size=" << P.pSolo.cbWLstr.size()
                         << " exact_map=" << InlineCBCorrection::exactMapSize()
                         << " variant_map=" << InlineCBCorrection::variantMapSize()
                         << " variant_collisions=" << InlineCBCorrection::variantCollisionSize()
                         << " collision_max_fanout=" << InlineCBCorrection::variantCollisionMaxFanout()
                         << "\n";
    }

    // Initialize global gene→probe cache before mapping (needed by inline resolver)
    if (transcriptomeMain != nullptr && P.pSolo.type != 0) {
        if (!P.pSolo.probeListPath.empty() && P.pSolo.probeListPath != "-") {
            ProbeListIndex probeIdxInit;
            if (probeIdxInit.load(P.pSolo.probeListPath)) {
                SoloFeature::initGeneProbeIdx(*transcriptomeMain, &probeIdxInit);
                P.inOut->logMain << "[GENE-PROBE-INIT] Pre-mapping init done from " << P.pSolo.probeListPath << "\n";
            } else {
                P.inOut->logMain << "[GENE-PROBE-INIT] WARNING: failed to load probe list at " << P.pSolo.probeListPath << "\n";
            }
        } else {
            P.inOut->logMain << "[GENE-PROBE-INIT] Skipped pre-mapping init (no probe list path)\n";
        }
    }

    // initialize Stats
    g_statsAll.resetN();
    time(&g_statsAll.timeStartMap);
    *P.inOut->logStdOut << timeMonthDayTime(g_statsAll.timeStartMap) << " ..... started mapping\n"
                        << flush;
    g_statsAll.timeLastReport = g_statsAll.timeStartMap;

    // SAM headers
    samHeaders(P, *genomeMain.genomeOut.g, *transcriptomeMain);

    // initialize chimeric parameters here - note that chimeric parameters require samHeader
    P.pCh.initialize(&P);

    // this does not seem to work at the moment
    // P.inOut->logMain << "mlock value="<<mlockall(MCL_CURRENT|MCL_FUTURE) <<"\n"<<flush;

    // prepare chunks and spawn mapping threads
    ReadAlignChunk *RAchunk[P.runThreadN];
    for (int ii = 0; ii < P.runThreadN; ii++)
    {
        RAchunk[ii] = new ReadAlignChunk(P, genomeMain, transcriptomeMain, ii);
    };

    if (P.runRestart.type != 1)
        mapThreadsSpawn(P, RAchunk);

    if (P.outFilterBySJoutStage == 1)
    { // completed stage 1, go to stage 2
        P.inOut->logMain << "Completed stage 1 mapping of outFilterBySJout mapping\n"
                         << flush;
        outputSJ(RAchunk, P); // collapse novel junctions
        P.readFilesIndex = -1;

        P.outFilterBySJoutStage = 2;
        if (P.outBAMcoord)
        {
            for (int it = 0; it < P.runThreadN; it++)
            { // prepare the unmapped bin
                RAchunk[it]->chunkOutBAMcoord->coordUnmappedPrepareBySJout();
            };
        };

        mapThreadsSpawn(P, RAchunk);
    };

    // close some BAM files
    if (P.inOut->outBAMfileUnsorted != NULL)
    {
        bgzf_flush(P.inOut->outBAMfileUnsorted);
        bgzf_close(P.inOut->outBAMfileUnsorted);
        if (P.emitNoYBAMyes) {
            if (P.inOut->outBAMfileY != NULL) bgzf_close(P.inOut->outBAMfileY);
            if (P.inOut->outBAMfileNoY != NULL) bgzf_close(P.inOut->outBAMfileNoY);
        }
    };
    if (P.inOut->outBAMfileUnsortedSoloTmp.is_open())
    {
        P.inOut->outBAMfileUnsortedSoloTmp.flush();
        P.inOut->outBAMfileUnsortedSoloTmp.close();
    };
    if (P.inOut->outQuantBAMfile != NULL)
    {
        bgzf_flush(P.inOut->outQuantBAMfile);
        bgzf_close(P.inOut->outQuantBAMfile);
    };

    if (P.outBAMcoord && P.limitBAMsortRAM == 0)
    { // make it equal ot the genome size
        P.limitBAMsortRAM = genomeMain.nGenome + genomeMain.SA.lengthByte + genomeMain.SAi.lengthByte;
    };

    time(&g_statsAll.timeFinishMap);
    *P.inOut->logStdOut << timeMonthDayTime(g_statsAll.timeFinishMap) << " ..... finished mapping\n"
                        << flush;
    P.inOut->logMain << timeMonthDayTime(g_statsAll.timeFinishMap) << " ..... finished mapping\n"
                     << "RAM after mapping:\n"
                     << linuxProcMemory() << flush;

    // no need for genome anymore, free the memory
    genomeMain.freeMemory();

    // aggregate output junctions
    // collapse splice junctions from different threads/chunks, and output them
    if (P.runRestart.type != 1 && P.outSJ.yes)
        outputSJ(RAchunk, P);

    // solo counts
    Solo soloMain(RAchunk, P, *RAchunk[0]->chunkTr);
    
    // Skip Solo processing if inline replayer already ran (it produces MEX directly)
    if (!P.pSolo.useInlineReplayer) {
        soloMain.processAndOutput();
    } else {
        P.inOut->logMain << timeMonthDayTime() << " ..... skipping Solo processing (inline replayer already produced MEX)" << endl;
    }

    // Note: Two-pass unsorted CB/UB tag injection removed - not used in inline flex path

    if (P.quant.geCount.yes)
    { // output gene quantifications
        for (int ichunk = 1; ichunk < P.runThreadN; ichunk++)
        { // sum counts from all chunks into 0th chunk
            RAchunk[0]->chunkTr->quants->addQuants(*(RAchunk[ichunk]->chunkTr->quants));
        };
        RAchunk[0]->chunkTr->quantsOutput();
    };

    if (P.runThreadN > 1 && P.outSAMorder == "PairedKeepInputOrder")
    { // concatenate Aligned.* files
        RAchunk[0]->chunkFilesCat(P.inOut->outSAM, P.outFileTmp + "/Aligned.out.sam.chunk", g_threadChunks.chunkOutN);
    };

    bamSortByCoordinate(P, RAchunk, *genomeMain.genomeOut.g, soloMain);

    // wiggle output
    if (P.outWigFlags.yes)
    {
        *(P.inOut->logStdOut) << timeMonthDayTime() << " ..... started wiggle output\n"
                              << flush;
        P.inOut->logMain << timeMonthDayTime() << " ..... started wiggle output\n"
                         << flush;
        string wigOutFileNamePrefix = P.outFileNamePrefix + "Signal";
        signalFromBAM(P.outBAMfileCoordName, wigOutFileNamePrefix, P);
    };

    g_statsAll.writeLines(P.inOut->outChimJunction, P.pCh.outJunctionFormat, "#", STAR_VERSION + string("   ") + P.commandLine);

    g_statsAll.progressReport(P.inOut->logProgress);
    P.inOut->logProgress << "ALL DONE!\n"
                         << flush;
    P.inOut->logFinal.open((P.outFileNamePrefix + "Log.final.out").c_str());
    g_statsAll.reportFinal(P.inOut->logFinal);
    *P.inOut->logStdOut << timeMonthDayTime(g_statsAll.timeFinish) << " ..... finished successfully\n"
                        << flush;

    P.inOut->logMain << "ALL DONE!\n"
                     << flush;
    if (P.outTmpKeep == "None")
    {
        sysRemoveDir(P.outFileTmp);
    };

    P.closeReadsFiles(); // this will kill the readFilesCommand processes if necessary
    // genomeMain.~Genome(); //need explicit call because of the 'delete P.inOut' below, which will destroy P.inOut->logStdOut
    if (genomeMain.sharedMemory != NULL)
    { // need explicit call because this destructor will write to files which are deleted by 'delete P.inOut' below
        delete genomeMain.sharedMemory;
        genomeMain.sharedMemory = NULL;
    };

    delete P.inOut; // to close files

    return 0;
};
