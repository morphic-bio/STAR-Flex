#include "ReadAlign.h"
#include "SampleDetector.h"
#include "GlobalVariables.h"
#include "ErrorWarning.h"
#include "SequenceFuns.h"
#include "UmiCodec.h"
#include "solo/CbCorrector.h"

namespace {
inline const char* sanitizeQname(char** readNameMates, uint32_t readNmates, uint32_t mateIdx) {
    if (mateIdx >= readNmates) {
        return nullptr;
    }
    const char* qptr = readNameMates[mateIdx];
    if (qptr && qptr[0] == '@') {
        ++qptr;
    }
    return qptr;
}
}

void ReadAlign::outputAlignments() {
  
    outBAMbytes=0;
        
    readAnnot.reset();
    
    if (mapGen.pGe.gType==101) {//temporary
        ReadAlign::spliceGraphWriteSAM();
        return;
    };    

    ReadAlign::outFilterBySJout();//sets outFilterBySJoutPass=false if read is held for the 2nd stage of outFilterBySJout
    
    if (outFilterBySJoutPass) {//otherwise align is held for the 2nd stage of outFilterBySJout
        ////////////////////////////////////
        if (unmapType<0) {//passed mappedFilter. Unmapped reads can have nTr>0

            auto nTr1 = nTr;
            auto trOut1 = trMult[0];

            if (P.pGe.transform.outYes) {
                nTr1 = alignsGenOut.alN;
                trOut1 = alignsGenOut.alMult[0];
            };

            if (nTr1>1) {//multimappers
                statsRA.mappedReadsM++;
                unmapType = -2; //not sure if this used
            } else if (nTr1==1) {//unique mappers
                statsRA.mappedReadsU++;
                statsRA.transcriptStats(*trOut1, Lread);
            };

            if (P.pGe.transform.outSAM && (!P.twoPass.yes || P.twoPass.pass2) ) {//transform genome only on 2nd pass
                ReadAlign::recordSJ(alignsGenOut.alN, alignsGenOut.alMult, chunkOutSJ);
            } else {
                ReadAlign::recordSJ(nTr, trMult, chunkOutSJ); //this will set mateMapped
            };            
            
            ReadAlign::alignedAnnotation();
        };
        
        // Y-chromosome alignment decision: check if any alignment touches Y
        // Mode-aware logic:
        // - Single-cell/Flex: R1/R2 are NOT mates, only check current read's alignments
        // - Bulk paired-end: Check both read's alignments AND mate's alignments (via mtid in BAM output)
        // - Bulk single-end: Only check read's own alignments
        hasYAlignment_ = false;
        if (P.emitNoYBAMyes) {
            // Determine if we're in single-cell/Flex mode (no mates) vs bulk paired-end (has mates)
            bool isSingleCellOrFlex = (P.pSolo.type != 0) || P.pSolo.flexMode;
            bool hasMates = (P.readNmates == 2) && !isSingleCellOrFlex;
            
            // Use transformed genome alignments if available, otherwise use original
            uint64 nTrCheck = nTr;
            Transcript **trCheck = trMult;
            if (P.pGe.transform.outSAM && (!P.twoPass.yes || P.twoPass.pass2)) {
                nTrCheck = alignsGenOut.alN;
                trCheck = alignsGenOut.alMult;
            }
            
            if (unmapType < 0) {
                // Mapped reads: check ALL transcripts (primary + secondary/supplementary)
                // Each transcript's Chr represents where that alignment maps
                // For single-cell/Flex: only current read's alignments (R1 or R2 processed separately)
                // For bulk paired-end: alignments for both mates are in trMult, check all
                for (uint iTr = 0; iTr < nTrCheck && !hasYAlignment_; iTr++) {
                    if (trCheck[iTr] == nullptr) continue;
                    // Check if this transcript aligns to Y chromosome
                    if (mapGen.yTids.count(trCheck[iTr]->Chr)) {
                        hasYAlignment_ = true;
                        break;
                    }
                }
                
                // For bulk paired-end: also check mate's reference (mtid) if available
                // Note: mtid is set during BAM writing, so we check it there if needed
                // For now, checking all transcripts should cover both mates in paired-end mode
            } else {
                // Unmapped reads: check if mate is mapped to Y (only for bulk paired-end)
                if (hasMates) {
                    // Bulk paired-end: check all transcripts to see if mate maps to Y
                    for (uint iTr = 0; iTr < nTrCheck && !hasYAlignment_; iTr++) {
                        if (trCheck[iTr] == nullptr) continue;
                        // Check if this transcript (for the mapped mate) is on Y
                        if (mapGen.yTids.count(trCheck[iTr]->Chr)) {
                            hasYAlignment_ = true;
                            break;
                        }
                    }
                    // Also check trBest as fallback (best alignment might be for the mate)
                    if (!hasYAlignment_ && trBest != nullptr) {
                        if (mapGen.yTids.count(trBest->Chr)) {
                            hasYAlignment_ = true;
                        }
                    }
                }
                // For single-cell/Flex unmapped reads or single-end unmapped: 
                // default hasYAlignment_ = false -> route to noY
            }
        }

        //the operations below are both for mapped and unmapped reads
        soloRead->readBar->getCBandUMI(Read0, Qual0, readLengthOriginal, readNameExtra[0], readFilesIndex, readName);
        
        // Extract CB/UMI/Sample metadata from Solo structures (upstream detection)
        // This avoids re-parsing BAM tags and ensures we use the original FASTQ data
        extractedCbIdxPlus1_ = 0;
        extractedUmi24_ = 0;
        
        // Store CB sequence for Phase 2 resolution (will be looked up in BAMoutput if needed)
        // Note: For ambiguous CBs, the sequence is stored in pendingAmbiguous_ and will be
        // looked up via resolvedCbByKey_ when writing keys
        
        // Use CbCorrector for inline CB correction if available
        SoloReadBarcode *readBar = soloRead->readBar;
        static int debugCount = 0;
        const int MAX_DEBUG_LOGS = 10;
        
        if (P.pSolo.cbCorrector && !readBar->cbSeq.empty()) {
            CbMatch match = P.pSolo.cbCorrector->correct(readBar->cbSeq);
            
            // Extract UMI first (needed for ambiguous accumulation)
            uint32_t umi24 = 0;
            if (readBar->urValid && readBar->urPacked != UINT32_MAX) {
                umi24 = readBar->urPacked & 0xFFFFFFu;
            } else if (!readBar->umiSeq.empty() && readBar->umiSeq.length() == 12) {
                uint32_t encoded = encodeUMI12(readBar->umiSeq);
                if (encoded != UINT32_MAX) {
                    umi24 = encoded & 0xFFFFFFu;
                }
            }
            
            if (match.ambiguous && !match.ambiguousIdx.empty() && umi24 != 0) {
                // Ambiguous CB: accumulate UMI counts for Bayesian resolution
                // This matches process_features: capture raw CB sequence (may contain Ns),
                // quality scores, candidate whitelist indices, and UMI counts
                AmbigKey key = hashCbSeq(readBar->cbSeq);
                auto &entry = pendingAmbiguous_[key];
                
                if (entry.candidateIdx.empty()) {
                    // First time seeing this ambiguous CB: initialize candidates
                    // Store raw observed CB sequence (e.g., "ACGNT...") and quality scores
                    entry.candidateIdx = match.ambiguousIdx; // 1-based whitelist indices
                    entry.cbSeq = readBar->cbSeq; // Raw observed sequence (may contain Ns)
                    entry.cbQual = readBar->cbQual; // Phred quality scores (same length as cbSeq)
                    entry.umiCounts.reserve(32);
                    
                    // Validate quality scores match sequence length
                    if (entry.cbQual.length() != entry.cbSeq.length()) {
                        // Pad with default quality if needed (shouldn't happen, but be defensive)
                        if (entry.cbQual.length() < entry.cbSeq.length()) {
                            entry.cbQual.append(entry.cbSeq.length() - entry.cbQual.length(), 'H'); // Q39 default
                        } else {
                            entry.cbQual = entry.cbQual.substr(0, entry.cbSeq.length());
                        }
                    }
                }
                
                // Accumulate UMI count (24-bit packed UMI -> count)
                entry.umiCounts[umi24]++;
                
                // Leave as 0 for now; will be resolved after mapping completes
                // Store CB sequence for lookup when writing keys
                extractedCbIdxPlus1_ = 0;
                extractedCbSeq_ = readBar->cbSeq;
            } else if (match.whitelistIdx != 0) {
                // Immediate resolution (non-ambiguous)
                extractedCbIdxPlus1_ = match.whitelistIdx;
                extractedCbSeq_.clear(); // Not needed for resolved CBs
                cbResolutionStats_.resolvedImmediate++;
            } else {
                // No match
                extractedCbIdxPlus1_ = 0;
                extractedCbSeq_.clear();
                cbResolutionStats_.noMatch++;
            }
            
            // Debug logging for first few reads and specific failing CBs
            if (debugCount++ < MAX_DEBUG_LOGS || readBar->cbSeq == "CNGTATTTCGGGCAGT") {
                P.inOut->logMain << "CbCorrector::correct: cbSeq=" << readBar->cbSeq
                                 << ", whitelistIdx=" << match.whitelistIdx
                                 << ", hammingDist=" << (int)match.hammingDist
                                 << ", ambiguous=" << match.ambiguous
                                 << ", extractedCbIdxPlus1_=" << extractedCbIdxPlus1_ << endl;
            }
        } else {
            // Debug logging for why CbCorrector wasn't used
            if (debugCount++ < MAX_DEBUG_LOGS) {
                P.inOut->logMain << "CbCorrector NOT used: cbCorrector=" << (P.pSolo.cbCorrector ? "non-null" : "null")
                                 << ", cbSeq.empty()=" << readBar->cbSeq.empty()
                                 << ", cbSeq=" << (readBar->cbSeq.empty() ? "(empty)" : readBar->cbSeq) << endl;
            }
            
            // Fallback to original SoloReadBarcode matching logic
            // cbMatch: 0=exact match, 1=one match with 1MM, -1=no match, -3=multiple matches (not allowed)
            if ((readBar->cbMatch == 0 || readBar->cbMatch == 1) && !readBar->cbMatchInd.empty()) {
                // Exact match or single mismatch match: use first match index (0-based) + 1
                extractedCbIdxPlus1_ = static_cast<uint32_t>(readBar->cbMatchInd[0] + 1);
            }
        }
        
        // Extract UMI (24-bit packed) from SoloReadBarcode (if not already extracted above)
        if (extractedUmi24_ == 0) {
            if (readBar->urValid && readBar->urPacked != UINT32_MAX) {
                extractedUmi24_ = readBar->urPacked & 0xFFFFFFu;
            } else if (!readBar->umiSeq.empty() && readBar->umiSeq.length() == 12) {
                uint32_t encoded = encodeUMI12(readBar->umiSeq);
                if (encoded != UINT32_MAX) {
                    extractedUmi24_ = encoded & 0xFFFFFFu;
                }
            }
        }
        
        // Debug: log first few extractions
        static int extractCount = 0;
        if (extractCount++ < 5) {
            P.inOut->logMain << "ReadAlign::outputAlignments: extractedCbIdxPlus1_=" << extractedCbIdxPlus1_
                             << ", extractedUmi24_=0x" << std::hex << extractedUmi24_ << std::dec
                             << ", cbMatch=" << readBar->cbMatch
                             << ", umiSeq.length()=" << readBar->umiSeq.length() << std::endl;
        }
        
        // Detect sample index from raw R2 sequence (mate 0) before alignment/BAM conversion
        // This ensures we detect from the original FASTQ sequence, not modified BAM sequence
        detectedSampleByte_ = 0xFFu; // Default: no sample (mate 1 or no detection)
        if (sampleDetReady_ && readNmates > 0 && readLengthOriginal[0] >= 8) {
            // Mate 0 = R2 = sample read (per convention: R2,R1 order)
            // Pack raw sequence to BAM format for detection
            uint32_t seqLen = readLengthOriginal[0];
            uint32_t packedLen = (seqLen + 1) / 2;
            uint8_t *packedSeq = new uint8_t[packedLen];
            nuclPackBAM(Read0[0], reinterpret_cast<char*>(packedSeq), seqLen);
            
            // Detect sample index (returns 1-based sequential index)
            // Note: reverseStrand=false for raw R2 sequences (they're sequenced forward)
            uint32_t detectedIdx = sampleDet_->detectSampleIndex(packedSeq, seqLen, false);
            
            static int debugCount = 0;
            if (debugCount++ < 5) {
                // Log first few bytes of raw sequence for debugging
                std::string seqPreview;
                for (int i = 0; i < std::min(20, static_cast<int>(seqLen)); i++) {
                    seqPreview += Read0[0][i];
                }
                // Log probe offset from ParametersSolo (we'll need to access it)
                P.inOut->logMain << "ReadAlign::outputAlignments: seqLen=" << seqLen
                                 << ", seqPreview=" << seqPreview
                                 << ", detectedIdx=" << detectedIdx << std::endl;
            }
            
            delete[] packedSeq;
            
            if (detectedIdx > 0) {
                // Store token (5 bits) for passing to BAMoutput
                // Token registration happens automatically inside detectSampleIndex()
                detectedSampleByte_ = static_cast<uint8_t>(detectedIdx & 0x1Fu);
                static int detectCount = 0;
                if (detectCount++ < 10) {
                    P.inOut->logMain << "ReadAlign::outputAlignments: detected sample idx=" << detectedIdx
                                     << ", token=" << (unsigned int)detectedSampleByte_ << std::endl;
                }
            } else {
                static int missCount = 0;
                if (missCount++ < 10) {
                    P.inOut->logMain << "ReadAlign::outputAlignments: detection returned 0 (no match)" << std::endl;
                }
            }
        } else {
            static int skipCount = 0;
            if (skipCount++ < 5) {
                P.inOut->logMain << "ReadAlign::outputAlignments: SKIP detection - sampleDetReady_="
                                 << sampleDetReady_ << ", readNmates=" << readNmates
                                 << ", readLengthOriginal[0]=" << (readNmates > 0 ? readLengthOriginal[0] : 0) << std::endl;
            }
        }

        //transcripts: need to be run after CB/UMI are obtained to output CR/UR tags
        if ( P.quant.trSAM.yes && unmapType<0) {//Aligned.toTranscriptome output, only for mapped
            if (P.pGe.transform.outQuant) {
                quantTranscriptome(chunkTr, alignsGenOut.alN,  alignsGenOut.alMult,  alignTrAll);
            } else {
                quantTranscriptome(chunkTr, nTr, trMult,  alignTrAll);
            };
        };        
        
        // Set detected sample token in SoloReadBarcode for tag extraction in inline hash capture
        if (soloRead && soloRead->readBar) {
            soloRead->readBar->detectedSampleToken = detectedSampleByte_;
        }

        // Populate optional MAPQ/CIGAR/score on transcripts for downstream consumers (inline resolver)
        if (unmapType < 0 && nTr > 0) {
            int mapqUnique = P.outSAMmapqUnique;
            for (uint64 i = 0; i < nTr; i++) {
                if (trMult[i] == nullptr) continue;
                // Heuristic MAPQ: unique vs multimapper
                trMult[i]->mapq = (nTr == 1) ? mapqUnique : 1;
                // Alignment score / mismatch analogues for downstream ranking (AS/NM-like)
                trMult[i]->asScore = static_cast<int>(trMult[i]->maxScore);
                trMult[i]->nm = static_cast<int>(trMult[i]->nMM);
                trMult[i]->cigarString = trMult[i]->generateCigarP();
                // store chrName for probe chromosome-derived gene IDs
                trMult[i]->chrName = genOut.chrName[trMult[i]->Chr];
            }
        }

        // Store qname mapping for reject logging if enabled
        // Forward declaration - function defined in SoloReadFeature_record.cpp
        extern void storeQnameMapping(uint64_t iRead, const char* qname);
        if (readName) {
            storeQnameMapping(iReadAll, readName);
        }
        
        soloRead->record((unmapType<0 ? nTr : 0), trMult, iReadAll, readAnnot); //need to supply nTr=0 for unmapped reads

        if (P.pGe.transform.outSAM) {
            ReadAlign::writeSAM(alignsGenOut.alN, alignsGenOut.alMult, alignsGenOut.alBest);
        } else {
            ReadAlign::writeSAM(nTr, trMult, trBest); //this will set mateMapped
        };
    };    

    if (unmapType>=0) {//unmapped reads
        statsRA.unmappedAll++; //include unmapType==4, i.e. one-mate alignments of PE reads - which may have been set in writeSAM above
        ReadAlign::outReadsUnmapped(); //uses mateMapped that was set in writeSAM above
    };    
};


///////////////////////////////////////////////////////////////////////////////////////////////////
void ReadAlign::recordSJ(uint64 nTrO, Transcript **trO, OutSJ *cSJ)
{//junction output for mapped reads (i.e. passed BySJout filtering)
    if (!P.outSJ.yes)
        return; //no SJ output
    
    if ( P.outSJfilterReads=="All" || nTrO==1 ) {
        uint64 sjReadStartN=cSJ->N;
        for (uint64 iTr=0; iTr<nTrO; iTr++) {//write all transcripts junctions
            outputTranscriptSJ (*(trO[iTr]), nTrO, cSJ, sjReadStartN);
        };
    };
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ReadAlign::outFilterBySJout()
{//filtering by SJout
    outFilterBySJoutPass=true;//only false if the alignment is held for outFilterBySJoutStage. True even if unmapped
    
    if (unmapType>0 || P.outFilterBySJoutStage!=1)
        return; //unmapped, or 2nd stage
   
    for (uint iTr=0;iTr<nTr;iTr++) {//check transcript for unannotated junctions
        for (uint iex=0;iex<trMult[iTr]->nExons-1;iex++) {//check all junctions
            if (trMult[iTr]->canonSJ[iex]>=0 && trMult[iTr]->sjAnnot[iex]==0) {
                outFilterBySJoutPass=false;
                break;
            };
        };
        if (!outFilterBySJoutPass) 
            break;
    };
    
    if (!outFilterBySJoutPass) {//this read is held for further filtering BySJout, record fastq
        unmapType=-3; //the read is not conisdered mapped
        statsRA.readN--;
        statsRA.readBases -= readLength[0]+readLength[1];

        for (uint im=0;im<P.readNends;im++) {
            chunkOutFilterBySJoutFiles[im] << readNameMates[im] <<" "<< iReadAll <<" "<< readFilter <<" "<< readFilesIndex;
            if (!readNameExtra[im].empty())
                chunkOutFilterBySJoutFiles[im]<<" "<< readNameExtra[im];
            chunkOutFilterBySJoutFiles[im] <<"\n";
            chunkOutFilterBySJoutFiles[im] << Read0[im] <<"\n";
            if (readFileType==2) {//fastq
                chunkOutFilterBySJoutFiles[im] << "+\n";
                chunkOutFilterBySJoutFiles[im] << Qual0[im] <<"\n";
            };
        };
    };
    
    //SJ output for all reads, including those not passed bySJout filtering. This only needs to be at the 1st stage of BySJout filtering
    ReadAlign::recordSJ(nTr, trMult, chunkOutSJ1); //this will set mateMapped
         
};

////////////////////////////////////////////////////////////////////////////////////////////////
void ReadAlign::writeSAM(uint64 nTrOutSAM, Transcript **trOutSAM, Transcript *trBestSAM)
{
    outBAMbytes=0;
    mateMapped[0] = mateMapped[1] = false; //mateMapped = are mates present in any of the transcripts?

    if (unmapType < 0 && outFilterBySJoutPass) {//write to SAM/BAM
        
        //////////////////////////////////////////////////////////////////////////////////
        /////////////outSAMfilter
        if (P.outSAMfilter.yes) {
            if (P.outSAMfilter.KeepOnlyAddedReferences) {
                for (uint itr=0;itr<nTrOutSAM;itr++) {//check if transcripts map to chr other than added references
                    if (trOutSAM[itr]->Chr<mapGen.genomeInsertChrIndFirst) {
                        return;//no SAM output
                    };
                };
            } else if (P.outSAMfilter.KeepAllAddedReferences) {
                uint64 nTrOutSAM1=0;
                for (uint itr=0;itr<nTrOutSAM;itr++) {//check if transcripts map to chr other than added references
                    if (trOutSAM[itr]->Chr>=mapGen.genomeInsertChrIndFirst) {
                        trOutSAM[nTrOutSAM1]=trOutSAM[itr];
                        trOutSAM[nTrOutSAM1]->primaryFlag=false;
                        ++nTrOutSAM1;
                    };
                };
                if (nTrOutSAM1==0) {
                   return;//no SAM output
                } else {
                    trOutSAM[0]->primaryFlag=true;
                };
                nTrOutSAM = nTrOutSAM1;
            };
        };
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////// write SAM/BAM 
        auto nTrOutWrite=min(P.outSAMmultNmax,nTrOutSAM); //number of aligns to write to SAM/BAM files            
        
        for (uint iTr=0;iTr<nTrOutWrite;iTr++) {//write transcripts
            //mateMapped1 = true if a mate is present in this transcript
            bool mateMapped1[2]={false,false};
            mateMapped1[trOutSAM[iTr]->exons[0][EX_iFrag]]=true;
            mateMapped1[trOutSAM[iTr]->exons[trOutSAM[iTr]->nExons-1][EX_iFrag]]=true;

            if (P.outSAMbool) {//SAM output
                outBAMbytes+=outputTranscriptSAM(*(trOutSAM[iTr]), nTrOutSAM, iTr, (uint) -1, (uint) -1, 0, -1, NULL, outSAMstream);
                if (P.outSAMunmapped.keepPairs && P.readNmates>1 && ( !mateMapped1[0] || !mateMapped1[1] ) ) {//keep pairs && paired reads && one of the mates not mapped in this transcript //not readNends: this is alignment
                    outBAMbytes+= outputTranscriptSAM(*(trOutSAM[iTr]), 0, 0, (uint) -1, (uint) -1, 0, 4, mateMapped1, outSAMstream);
                };
            };

            if (P.outBAMunsorted || P.outBAMcoord) {//BAM output
                alignBAM(*(trOutSAM[iTr]), nTrOutSAM, iTr, mapGen.chrStart[trOutSAM[iTr]->Chr], (uint) -1, (uint) -1, 0, -1, NULL, P.outSAMattrOrder,outBAMoneAlign, outBAMoneAlignNbytes);

                if (P.outBAMunsorted && outBAMunsorted != NULL) {//unsorted mode
                    for (uint imate=0; imate<P.readNmates; imate++) {//output each mate //not readNends: this is alignment
                        outBAMunsorted->unsortedOneAlign(
                            outBAMoneAlign[imate],
                            outBAMoneAlignNbytes[imate],
                            (imate>0 || iTr>0) ? 0 : (outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1])*2*nTrOutWrite,
                            iReadAll,
                            (imate==0 ? detectedSampleByte_ : 0xFFu),
                            (imate==0 ? extractedCbIdxPlus1_ : 0u),
                            (imate==0 ? extractedUmi24_ : 0u),
                            (imate==0 && extractedCbIdxPlus1_ == 0 ? extractedCbSeq_ : std::string()),
                            hasYAlignment_);
                    };
                    if (P.outSAMunmapped.keepPairs && P.readNmates>1 && ( !mateMapped1[0] || !mateMapped1[1] ) ) {//keep pairs && paired reads && one of the mates not mapped in this transcript //not readNends: this is alignment
                        alignBAM(*trOutSAM[iTr], 0, 0, mapGen.chrStart[trOutSAM[iTr]->Chr], (uint) -1, (uint) -1, 0, 4, mateMapped1, P.outSAMattrOrder, outBAMoneAlign, outBAMoneAlignNbytes);
                        for (uint imate=0; imate<P.readNmates; imate++) {//output each mate //not readNends: this is alignment
                            outBAMunsorted->unsortedOneAlign(
                                outBAMoneAlign[imate],
                                outBAMoneAlignNbytes[imate],
                                (imate>0 || iTr>0) ? 0 : (outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1])*2*nTrOutWrite,
                                iReadAll,
                                (imate==0 ? detectedSampleByte_ : 0xFFu),
                                0u, 0u, std::string(), // CB/UMI not available for unmapped pairs
                                hasYAlignment_);
                        };
                    };
                };

                if (P.outBAMcoord) {//coordinate sorted
                    for (uint imate=0; imate<P.readNmates; imate++) {//output each mate //not readNends: this is alignment
                        outBAMcoord->coordOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], (iReadAll<<32) | (iTr<<8) | trOutSAM[iTr]->exons[0][EX_iFrag], hasYAlignment_);
                    };
                };
            };
        };

        /////////////////////////////////////////////////////////////////////////////////////////////
        //////// write unmapped ends
        //TODO it's better to check all transcripts in the loop above for presence of both mates
        mateMapped[trBestSAM->exons[0][EX_iFrag]] = true;
        mateMapped[trBestSAM->exons[trBestSAM->nExons-1][EX_iFrag]] = true;

        if (P.readNmates>1 && !(mateMapped[0] && mateMapped[1]) ) {//not readNends: this is alignment
            unmapType=4;
        };

        if (unmapType==4 && P.outSAMunmapped.within) {//output unmapped ends for single-end alignments of PE reads
            if (P.outSAMbool && !P.outSAMunmapped.keepPairs ) {
                outBAMbytes+= outputTranscriptSAM(*trBestSAM, 0, 0, (uint) -1, (uint) -1, 0, unmapType, mateMapped, outSAMstream);
            };

            if ( P.outBAMcoord || (P.outBAMunsorted && !P.outSAMunmapped.keepPairs) ) {//BAM output
                alignBAM(*trBestSAM, 0, 0, mapGen.chrStart[trBestSAM->Chr], (uint) -1, (uint) -1, 0, unmapType, mateMapped, P.outSAMattrOrder, outBAMoneAlign, outBAMoneAlignNbytes);
                for (uint imate=0; imate<P.readNmates; imate++) {//alignBAM output is empty for mapped mate, but still need to scan through it //not readNends: this is alignment
                    if (P.outBAMunsorted && outBAMunsorted != NULL && !P.outSAMunmapped.keepPairs) {
                        outBAMunsorted->unsortedOneAlign(outBAMoneAlign[imate],
                                                          outBAMoneAlignNbytes[imate],
                                                          imate>0 ? 0 : outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1],
                                                          iReadAll,
                                                          (imate==0 ? detectedSampleByte_ : 0xFFu),
                                                          0u, 0u, std::string(), // CB/UMI not available for unmapped pairs
                                                          hasYAlignment_);
                    };
                    if (P.outBAMcoord) {//KeepPairs option does not affect for sorted BAM since we do not want multiple entries for the same unmapped read
                        outBAMcoord->coordOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], iReadAll<<32, hasYAlignment_);
                    };
                };
            };
        };  
        
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////write completely unmapped
    } else if (unmapType>=0 && P.outSAMunmapped.within) {//output unmapped within && unmapped read && both mates unmapped
        if (P.outBAMcoord || P.outBAMunsorted || P.quant.trSAM.bamYes) {//BAM output
            alignBAM(*trBestSAM, 0, 0, mapGen.chrStart[trBestSAM->Chr], (uint) -1, (uint) -1, 0, unmapType, mateMapped, P.outSAMattrOrder, outBAMoneAlign, outBAMoneAlignNbytes);
            for (uint imate=0; imate<P.readNmates; imate++) {//output each mate //not readNends: this is alignment
                if (P.outBAMunsorted && outBAMunsorted != NULL) {
                    outBAMunsorted->unsortedOneAlign(outBAMoneAlign[imate],
                                                      outBAMoneAlignNbytes[imate],
                                                      imate>0 ? 0 : outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1],
                                                      iReadAll,
                                                      (imate==0 ? detectedSampleByte_ : 0xFFu),
                                                      (imate==0 ? extractedCbIdxPlus1_ : 0u),
                                                      (imate==0 ? extractedUmi24_ : 0u),
                                                      (imate==0 && extractedCbIdxPlus1_ == 0 ? extractedCbSeq_ : std::string()),
                                                      hasYAlignment_);
                };
                if (P.quant.trSAM.bamYes) {
                    outBAMquant->unsortedOneAlign(outBAMoneAlign[imate],
                                                  outBAMoneAlignNbytes[imate],
                                                  imate>0 ? 0 : outBAMoneAlignNbytes[0]+outBAMoneAlignNbytes[1],
                                                  iReadAll,
                                                  (imate==0 ? detectedSampleByte_ : 0xFFu),
                                                  (imate==0 ? extractedCbIdxPlus1_ : 0u),
                                                  (imate==0 ? extractedUmi24_ : 0u),
                                                  (imate==0 && extractedCbIdxPlus1_ == 0 ? extractedCbSeq_ : std::string()));
                };
                if (P.outBAMcoord) {
                    outBAMcoord->coordOneAlign(outBAMoneAlign[imate], outBAMoneAlignNbytes[imate], iReadAll<<32, hasYAlignment_);
                };
            };
        };

        if (P.outSAMbool) {//output SAM
            outBAMbytes+= outputTranscriptSAM(*trBestSAM, 0, 0, (uint) -1, (uint) -1, 0, unmapType, mateMapped, outSAMstream);
        };
    };       
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ReadAlign::outReadsUnmapped()
{
    if (P.outReadsUnmapped=="Fastx" ) {//output to fasta/q files. Include unmapType==4, i.e. one-mate alignments of PE reads
       for (uint im=0;im<P.readNends;im++) {
           chunkOutUnmappedReadsStream[im] << readNameMates[im]  <<" "<<im<<":"<< readFilter <<": "<< readNameExtra[im];
           if (P.readNmates>1) //not readNends: this is alignment
               chunkOutUnmappedReadsStream[im] <<" "<< int(mateMapped[0]) <<  int(mateMapped[1]);
           chunkOutUnmappedReadsStream[im] <<"\n";
           chunkOutUnmappedReadsStream[im] << Read0[im] <<"\n";
            if (readFileType==2) {//fastq
                chunkOutUnmappedReadsStream[im] << "+\n";
                chunkOutUnmappedReadsStream[im] << Qual0[im] <<"\n";
            };
       };
    };
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ReadAlign::spliceGraphWriteSAM()
{//temporary: SAM output for SpliceGraph
    outBAMbytes=0;
    uint64 nTrOutSAM = nTr;
    if (mapGen.genomeOut.convYes) {//convert to new genome
        nTrOutSAM=0;
        for (uint32 iTr=0; iTr<nTrOutSAM; iTr++) {//convert output transcripts into new genome
            *alignsGenOut.alMult[nTrOutSAM] = *trMult[iTr];//copy information before conversion
            if (trMult[iTr]->convertGenomeCigar(*mapGen.genomeOut.g, *alignsGenOut.alMult[nTrOutSAM])) {
                ++nTrOutSAM;
                trMult[nTrOutSAM-1] = alignsGenOut.alMult[nTrOutSAM-1]; //point to new transcsript
            };
        };
        nTrOutSAM=nTrOutSAM;
    };

    for (uint iTr=0; iTr<nTrOutSAM; iTr++) {//write all transcripts            
        outBAMbytes += outputSpliceGraphSAM(*(trMult[iTr]), nTrOutSAM, iTr, outSAMstream);
    };
};

void ReadAlign::alignedAnnotation()
{
    //TODO maybe initialize readAnnot to all empty?
    //genes
    if ( P.quant.geCount.yes ) {
        if (P.pGe.transform.outQuant) {
            chunkTr->geneCountsAddAlign(alignsGenOut.alN, alignsGenOut.alMult, readAnnot.geneExonOverlap);
        } else {
            chunkTr->geneCountsAddAlign(nTr, trMult, readAnnot.geneExonOverlap);
        };        
    };
    //solo-GeneFull
    if ( P.quant.geneFull.yes ) {
        chunkTr->geneFullAlignOverlap(nTr, trMult, P.pSolo.strand, readAnnot.annotFeatures[SoloFeatureTypes::GeneFull]);
    };   
    //solo-Gene
    if ( P.quant.gene.yes ) {
        chunkTr->classifyAlign(trMult, nTr, readAnnot);
    };
    //solo-GeneFull_ExonOverIntron
    if ( P.quant.geneFull_ExonOverIntron.yes ) {
        chunkTr->geneFullAlignOverlap_ExonOverIntron(nTr, trMult, P.pSolo.strand, readAnnot.annotFeatures[SoloFeatureTypes::GeneFull_ExonOverIntron], readAnnot.annotFeatures[SoloFeatureTypes::Gene]);
    };
    //solo-GeneFull_Ex50pAS
    if ( P.quant.geneFull_Ex50pAS.yes ) {
        chunkTr->alignExonOverlap(nTr, trMult, P.pSolo.strand, readAnnot.annotFeatures[SoloFeatureTypes::GeneFull_Ex50pAS]);
    };    
};
