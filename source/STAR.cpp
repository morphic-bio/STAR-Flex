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
#include "SamtoolsSorter.h"
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
#include "TranscriptQuantEC.h"
#include "LibFormatDetection.h"
#include "vb_engine.h"
#include "em_engine.h"
#include "ec_loader.h"
#include "gc_bias.h"
#include "fld_accumulator.h"
#include "TranscriptQuantOutput.h"
// Note: effective_length.h not included due to Transcriptome class name conflict
// Use wrapper function instead
#include "effective_length_wrapper.h"
#include "InlineCBCorrection.h"
#include "alignment_model.h"  // For Transcriptome and AlignmentModel
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

    // Shared libem Transcriptome sequence cache for error model (read-only, shared across threads)
    std::unique_ptr<libem::Transcriptome> libem_transcriptome;
    
    if (P.quant.yes)
    { // load transcriptome
        transcriptomeMain = new Transcriptome(P);
        
        // Load transcript sequences for error model if enabled
        if (P.quant.transcriptVB.yes && P.quant.transcriptVB.errorModelMode != "off") {
            // Determine FASTA path: P.pGe.transcriptomeFasta else transcriptome.fa
            std::string fasta_path;
            if (!P.pGe.transcriptomeFasta.empty() && P.pGe.transcriptomeFasta != "-") {
                fasta_path = P.pGe.transcriptomeFasta;
            } else {
                // Try transcriptome.fa in genome directory
                fasta_path = P.pGe.gDir + "/transcriptome.fa";
            }
            
            libem_transcriptome.reset(new libem::Transcriptome());
            if (!libem_transcriptome->loadFromFasta(fasta_path)) {
                // Failed to load - disable error model
                P.inOut->logMain << "WARNING: Failed to load transcript sequences from " << fasta_path 
                                 << ". Error model will be disabled.\n";
                libem_transcriptome.reset();
            } else {
                // Reorder by STAR transcript names to match BAM header order
                std::vector<std::string> star_names(transcriptomeMain->nTr);
                for (uint i = 0; i < transcriptomeMain->nTr; ++i) {
                    star_names[i] = transcriptomeMain->trID[i];
                }
                if (!libem_transcriptome->reorderByNames(star_names)) {
                    P.inOut->logMain << "WARNING: Failed to reorder transcript sequences to match STAR order. "
                                     << "Error model may use incorrect sequences.\n";
                }
                P.inOut->logMain << "Loaded " << libem_transcriptome->size() 
                                 << " transcript sequences for error model\n";
            }
        }
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
            uint32_t deprecatedCount = 0;
            if (probeIdxInit.load(P.pSolo.probeListPath, P.pSolo.removeDeprecated, &deprecatedCount)) {
                SoloFeature::initGeneProbeIdx(*transcriptomeMain, &probeIdxInit);
                P.inOut->logMain << "[GENE-PROBE-INIT] Pre-mapping init done from " << P.pSolo.probeListPath << "\n";
                if (P.pSolo.removeDeprecated && deprecatedCount > 0) {
                    P.inOut->logMain << "[GENE-PROBE-INIT] Removed " << deprecatedCount << " deprecated entries from probe list\n";
                }
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

    // Apply limitBAMsortRAM fallback before creating SamtoolsSorter
    if (P.outBAMcoord && P.limitBAMsortRAM == 0) {
        // make it equal to the genome size
        P.limitBAMsortRAM = genomeMain.nGenome + genomeMain.SA.lengthByte + genomeMain.SAi.lengthByte;
    }

    // Initialize samtools sorter if needed
    if (P.outBAMsortMethod == "samtools" && P.outBAMcoord) {
        g_samtoolsSorter = new SamtoolsSorter(P.limitBAMsortRAM, 
                                               P.outBAMsortingThreadNactual,
                                               P.outBAMsortTmpDir,
                                               P);
    }

    // this does not seem to work at the moment
    // P.inOut->logMain << "mlock value="<<mlockall(MCL_CURRENT|MCL_FUTURE) <<"\n"<<flush;

    // prepare chunks and spawn mapping threads
    ReadAlignChunk *RAchunk[P.runThreadN];
    for (int ii = 0; ii < P.runThreadN; ii++)
    {
        RAchunk[ii] = new ReadAlignChunk(P, genomeMain, transcriptomeMain, ii, 
                                          libem_transcriptome.get());  // Pass shared transcriptome
    };

    // === LIBRARY FORMAT DETECTION (single-threaded) ===
    if (P.quant.transcriptVB.yes && P.quant.transcriptVB.libType == "A") {
        P.inOut->logMain << "Starting library format auto-detection "
                         << "(first " << P.quant.transcriptVB.autoDetectWindow 
                         << " reads)...\n" << flush;
        
        // Create shared detector (will be accessed by TranscriptQuantEC during voting)
        P.quant.transcriptVB.libFormatDetector = new LibFormatDetector(
            P.quant.transcriptVB.autoDetectWindow);
        
        // Set detection mode BEFORE processing
        P.quant.transcriptVB.inDetectionMode = true;
        
        // Temporarily set read limit to detection window
        uint64_t originalLimit = P.readMapNumber;
        P.readMapNumber = P.quant.transcriptVB.autoDetectWindow;
        
        // Process first N reads using RAchunk[0] single-threaded
        // Voting happens INSIDE addReadAlignments() when inDetectionMode=true
        // NOTE: If total reads <= autoDetectWindow, RAchunk[0] will process all reads here
        // and then mapThreadsSpawn() will run again on RAchunk[0] with no remaining reads.
        // This is fine for datasets larger than autoDetectWindow, but avoid very small
        // read counts with auto-detect enabled to prevent potential duplicate output.
        RAchunk[0]->processChunks();
        
        // Restore original limit
        P.readMapNumber = originalLimit;
        
        // Finalize detection - HARD FAILURE if ambiguous (exits before returning)
        LibraryFormat detected = P.quant.transcriptVB.libFormatDetector->finalizeOrFail(
            P.inOut->logMain);
        
        // Store detected format as uint8_t
        P.quant.transcriptVB.detectedLibFormatId = detected.typeId();
        P.quant.transcriptVB.detectionComplete = true;
        P.quant.transcriptVB.inDetectionMode = false;  // Detection done
        
        // formatName() is defined in LibFormatDetection.h (see Section 9)
        P.inOut->logMain << "Detected library format: " << formatName(detected)
                         << " (formatID=" << static_cast<int>(detected.typeId()) << ")"
                         << " from " << P.iReadAll << " reads\n" << flush;
        
        // Clean up detector (no longer needed)
        delete P.quant.transcriptVB.libFormatDetector;
        P.quant.transcriptVB.libFormatDetector = nullptr;
        
        // Reset global counters after detection completes
        // This ensures pre-burn-in gating starts fresh for the main mapping pass (matches Salmon's behavior)
        extern std::atomic<uint64_t> global_processed_fragments;
        global_processed_fragments.store(0, std::memory_order_relaxed);
        Parameters::global_fld_obs_count.store(0, std::memory_order_relaxed);
    }

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
    
    // Transcript quantification (TranscriptVB mode)
    // TODO: Debug - check if transcriptomeMain is valid
    if (P.quant.transcriptVB.yes && transcriptomeMain != nullptr && transcriptomeMain->nTr > 0) {
        *P.inOut->logStdOut << timeMonthDayTime() << " ..... started transcript quantification\n"
                          << flush;
        P.inOut->logMain << timeMonthDayTime() << " ..... started transcript quantification\n";
        
        // 1. Merge EC tables from all threads
        TranscriptQuantEC mergedEC(transcriptomeMain->nTr, 0, "", 0, P);
        for (int ichunk = 0; ichunk < P.runThreadN; ichunk++) {
            if (RAchunk[ichunk] != nullptr && 
                RAchunk[ichunk]->RA != nullptr && 
                RAchunk[ichunk]->RA->quantEC != nullptr) {
                mergedEC.merge(*RAchunk[ichunk]->RA->quantEC);
            }
        }
        
        mergedEC.finalize();
        
        *P.inOut->logStdOut << "Merged " << mergedEC.getECTable().n_ecs << " equivalence classes from " << P.runThreadN << " threads\n"
                          << flush;
        P.inOut->logMain << "Merged " << mergedEC.getECTable().n_ecs << " equivalence classes from " << P.runThreadN << " threads\n";
        
        // Log drop statistics (per trace plan Step 3)
        P.inOut->logMain << "EC building drop statistics:\n"
                         << "  dropped_incompat: " << mergedEC.getDroppedIncompat() << "\n"
                         << "  dropped_missing_mate_fields: " << mergedEC.getDroppedMissingMateFields() << "\n"
                         << "  dropped_unknown_obs_fmt: " << mergedEC.getDroppedUnknownObsFmt() << "\n";
        
        // Note: Transcript lengths are already set per-thread in ReadAlign.cpp constructor
        // and carried through merge(), so no need to set them again here.
        
        // 2. Initialize transcript state
        TranscriptState state;
        state.resize(transcriptomeMain->nTr);
        
        // Compute effective lengths using FLD (if available) or fallback to mean fragment length
        const FLDAccumulator& observedFLD = mergedEC.getObservedFLD();
        double meanFragLen = 200.0;  // Default fallback
        bool use_fld = observedFLD.isValid();
        
        if (use_fld) {
            meanFragLen = observedFLD.getMean();
            double fragLenStdDev = observedFLD.getStdDev();
            *P.inOut->logStdOut << "Fragment length distribution: mean=" << meanFragLen 
                              << ", stddev=" << fragLenStdDev 
                              << ", fragments=" << observedFLD.getTotalFragments() << "\n"
                              << flush;
            P.inOut->logMain << "Fragment length distribution: mean=" << meanFragLen 
                           << ", stddev=" << fragLenStdDev 
                           << ", fragments=" << observedFLD.getTotalFragments() << "\n";
        }
        
        // Log GC observations if GC bias is enabled
        if (P.quant.transcriptVB.gcBias) {
            const GCFragModel& observedGC = mergedEC.getObservedGC();
            const auto& gcCounts = observedGC.getCounts();
            double totalObs = 0.0;
            for (int i = 0; i < 101; i++) {
                totalObs += gcCounts[i];
            }
            if (totalObs > 100) {
                *P.inOut->logStdOut << "GC bias: collected " << static_cast<uint64_t>(totalObs) << " fragment observations\n"
                                  << flush;
                P.inOut->logMain << "GC bias: collected " << static_cast<uint64_t>(totalObs) << " fragment observations\n";
            }
        }
        
        // Compute effective lengths using FLD PMF (if available) or simple mean-based calculation
        std::vector<double> fld_pmf;
        if (use_fld) {
            fld_pmf = observedFLD.getPMF();
        }
        
        // Build raw_lengths_int vector
        std::vector<int32_t> raw_lengths_int(transcriptomeMain->nTr);
        for (uint i = 0; i < transcriptomeMain->nTr; ++i) {
            raw_lengths_int[i] = static_cast<int32_t>(transcriptomeMain->trLen[i]);
        }
        
        // Compute effective lengths
        std::vector<double> eff_lengths;
        if (use_fld && !fld_pmf.empty()) {
            // Use EffectiveLengthCalculator via wrapper (avoids Transcriptome conflict)
            eff_lengths = computeEffectiveLengthsFromPMFWrapper(fld_pmf, raw_lengths_int);
        } else {
            // Fallback: simple mean-based calculation
            eff_lengths.resize(transcriptomeMain->nTr);
            for (uint i = 0; i < transcriptomeMain->nTr; ++i) {
                double rawLen = static_cast<double>(transcriptomeMain->trLen[i]);
                double effLen = rawLen - meanFragLen + 1.0;
                if (effLen < 1.0) effLen = 1.0;
                if (effLen > rawLen) effLen = rawLen;
                eff_lengths[i] = effLen;
            }
        }
        
        // Populate transcript state
        for (uint i = 0; i < transcriptomeMain->nTr; ++i) {
            state.names[i] = transcriptomeMain->trID[i];
            double rawLen = static_cast<double>(transcriptomeMain->trLen[i]);
            state.lengths[i] = rawLen;
            state.eff_lengths[i] = eff_lengths[i];
        }
        
        // Note: GC-corrected effective lengths would require transcript sequences
        // For now, FLD-based effective lengths are computed above
        if (P.quant.transcriptVB.gcBias) {
            *P.inOut->logStdOut << "GC bias: GC correction requires transcript sequences (not yet implemented)\n"
                              << flush;
            P.inOut->logMain << "GC bias: GC correction requires transcript sequences (not yet implemented)\n";
        }
        
        // 4. Run VB/EM quantification
        EMParams params;
        params.use_vb = P.quant.transcriptVB.vb;
        params.vb_prior = P.quant.transcriptVB.vbPrior;
        // Use defaults from EMParams (max_iters=10000, min_iters=100, tolerance=0.01)
        // Do NOT override for VB - let VB use same defaults as EM for Salmon parity
        // Thread count: use OMP default unless explicitly set (we'll pass --runThreadN externally for parity tests)
        params.threads = 0;  // 0 = use OMP default (multi-thread capable)
        
        EMResult result;
        if (params.use_vb) {
            result = run_vb(mergedEC.getECTable(), state, params);
        } else {
            result = run_em(mergedEC.getECTable(), state, params);
        }
        
        *P.inOut->logStdOut << "Quantification converged: " << (result.converged ? "yes" : "no")
                          << ", iterations: " << result.iterations << "\n"
                          << flush;
        P.inOut->logMain << "Quantification converged: " << (result.converged ? "yes" : "no")
                         << ", iterations: " << result.iterations << "\n";
        
        // 5. Write quant.sf
        writeQuantSF(result, state, P.quant.transcriptVB.outFile);
        
        // 6. Write quant.genes.sf (gene-level aggregation, Legacy mode)
        if (P.quant.transcriptVB.geneOutput) {
            int geneResult = writeQuantGeneSF(result, state, *transcriptomeMain, 
                                               P.quant.transcriptVB.outFileGene);
            
            if (geneResult == 1) {
                // File open error
                P.inOut->logMain << "ERROR: Failed to open gene output file: " 
                                 << P.quant.transcriptVB.outFileGene << "\n";
            } else {
                P.inOut->logMain << "Gene quantification written to: " 
                                 << P.quant.transcriptVB.outFileGene << "\n";
                
                if (geneResult == 2) {
                    // MissingGeneID warning
                    P.inOut->logMain << "WARNING: transcripts with MissingGeneID were aggregated "
                                     << "to a single gene entry in quant.genes.sf\n";
                }
            }
            
            // 6b. Write tximport-style gene output if enabled
            if (P.quant.transcriptVB.genesTximport) {
                int tximportResult = writeQuantGeneSFTximport(result, state, *transcriptomeMain,
                                                              P.quant.transcriptVB.outFileGeneTximport);
                
                if (tximportResult == 1) {
                    P.inOut->logMain << "ERROR: Failed to open tximport gene output file: "
                                     << P.quant.transcriptVB.outFileGeneTximport << "\n";
                } else {
                    P.inOut->logMain << "Gene quantification (tximport mode) written to: "
                                     << P.quant.transcriptVB.outFileGeneTximport << "\n";
                }
            }
        }
        
        *P.inOut->logStdOut << timeMonthDayTime() << " ..... finished transcript quantification\n"
                          << flush;
        P.inOut->logMain << timeMonthDayTime() << " ..... finished transcript quantification\n";
    }

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

    // Cleanup samtools sorter if it exists
    if (g_samtoolsSorter != nullptr) {
        delete g_samtoolsSorter;
        g_samtoolsSorter = nullptr;
    }

    delete P.inOut; // to close files

    return 0;
};
