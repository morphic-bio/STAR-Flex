#include "ReadAlign.h"
#include "SlamQuant.h"

namespace {
struct GenomicInterval {
    uint64_t start;
    uint64_t end;
};

inline void addInterval(std::vector<GenomicInterval>& intervals, uint64_t start, uint64_t end) {
    if (start <= end) {
        intervals.push_back({start, end});
    }
}

inline bool containsPos(const std::vector<GenomicInterval>& intervals, size_t& idx, uint64_t pos) {
    while (idx < intervals.size() && pos > intervals[idx].end) {
        ++idx;
    }
    return idx < intervals.size() && pos >= intervals[idx].start && pos <= intervals[idx].end;
}

inline bool isOppositeStrand(const Transcriptome& tr, const std::set<uint32_t>& geneIds, uint8_t readStr) {
    bool anySense = false;
    bool anyOpposite = false;
    for (uint32_t geneId : geneIds) {
        if (geneId >= tr.geStr.size()) {
            continue;
        }
        uint8_t geneStr = tr.geStr[geneId];
        if (geneStr == 0) {
            continue;
        }
        bool sense = (geneStr == static_cast<uint8_t>(readStr + 1));
        anySense = anySense || sense;
        anyOpposite = anyOpposite || !sense;
        if (anySense && anyOpposite) {
            break;
        }
    }
    if (anySense) {
        return false;
    }
    return anyOpposite;
}
} // namespace

bool ReadAlign::slamCollect(const Transcript& trOut, const std::set<uint32_t>& geneIds, double weight, bool isIntronic) {
    if (slamQuant == nullptr || weight <= 0.0) {
        return false;
    }
    if (geneIds.empty()) {
        if (!isIntronic) {
            slamQuant->diagnostics().readsZeroGenes++;
        }
        return false;
    }

    char* R = Read1[trOut.roStr == 0 ? 0 : 2];
    bool isMinus = (trOut.Str == 1);
    bool oppositeStrand = isOppositeStrand(*chunkTr, geneIds, trOut.Str);
    SlamMismatchCategory category = isIntronic ? SlamMismatchCategory::Intronic : SlamMismatchCategory::Exonic;
    SlamMismatchCategory senseCategory = isIntronic ? SlamMismatchCategory::IntronicSense : SlamMismatchCategory::ExonicSense;
    bool snpDetect = slamQuant->snpDetectEnabled() && !isIntronic;
    std::vector<uint32_t> mismatchPositions;

    uint16_t nT = 0;
    uint16_t k = 0;

    if (slamSnpMask != nullptr) {
        for (uint iex = 0; iex < trOut.nExons; ++iex) {
            uint64 gStart = trOut.exons[iex][EX_G];
            uint64 len = trOut.exons[iex][EX_L];
            for (uint64 ii = 0; ii < len; ++ii) {
                uint64 gpos = gStart + ii;
                if (slamSnpMask->contains(gpos)) {
                    slamQuant->diagnostics().readsDroppedSnpMask++;
                    return false; // discard read if it overlaps SNP mask
                }
            }
        }
    }

    std::vector<GenomicInterval> mate1Intervals;
    std::vector<GenomicInterval> mate2Intervals;
    bool hasMate = (readLength[1] > 0);
    if (hasMate) {
        mate1Intervals.reserve(trOut.nExons);
        mate2Intervals.reserve(trOut.nExons);
        uint32_t split = readLength[0];
        for (uint iex = 0; iex < trOut.nExons; ++iex) {
            uint64 gStart = trOut.exons[iex][EX_G];
            uint64 rStart = trOut.exons[iex][EX_R];
            uint64 len = trOut.exons[iex][EX_L];
            uint64 rEnd = rStart + len;
            if (rEnd <= split) {
                addInterval(mate1Intervals, gStart, gStart + len - 1);
            } else if (rStart > split) {
                addInterval(mate2Intervals, gStart, gStart + len - 1);
            } else {
                if (rStart < split) {
                    uint64 len1 = split - rStart;
                    addInterval(mate1Intervals, gStart, gStart + len1 - 1);
                    uint64 len2 = len - len1;
                    if (len2 > 0) {
                        addInterval(mate2Intervals, gStart + len1, gStart + len1 + len2 - 1);
                    }
                }
            }
        }
    }
    size_t mate1Idx = 0;
    size_t mate2Idx = 0;

    for (uint iex = 0; iex < trOut.nExons; ++iex) {
        uint64 gStart = trOut.exons[iex][EX_G];
        uint64 rStart = trOut.exons[iex][EX_R];
        uint64 len = trOut.exons[iex][EX_L];

        for (uint64 ii = 0; ii < len; ++ii) {
            uint64 gpos = gStart + ii;
            uint8_t r1 = static_cast<uint8_t>(R[rStart + ii]);
            uint8_t g1 = static_cast<uint8_t>(genOut.G[gpos]);
            if (r1 > 3 || g1 > 3) {
                continue;
            }
            if (snpDetect) {
                slamQuant->recordSnpObservation(gpos, g1 != r1);
            }
            uint64 rposRaw = rStart + ii;
            if (readLength[0] > 0 && rposRaw == readLength[0]) {
                continue; // skip spacer between mates
            }
            bool secondMate = (readLength[1] > 0 && rposRaw > readLength[0]);
            uint32_t readPos = static_cast<uint32_t>(secondMate ? rposRaw - 1 : rposRaw);
            bool overlap = false;
            if (hasMate) {
                if (secondMate) {
                    overlap = containsPos(mate1Intervals, mate1Idx, gpos);
                } else {
                    overlap = containsPos(mate2Intervals, mate2Idx, gpos);
                }
            }
            bool skipMismatch = overlap && secondMate;
            if (!skipMismatch) {
                slamQuant->addTransitionBase(category, readPos, secondMate, overlap, oppositeStrand, g1, r1, weight);
                if (!oppositeStrand) {
                    slamQuant->addTransitionBase(senseCategory, readPos, secondMate, overlap, false, g1, r1, weight);
                }
            }
            if (!isIntronic) {
                if (!isMinus) {
                    if (g1 == 3) { // T
                        ++nT;
                        if (r1 == 1) { // C
                            ++k;
                            if (snpDetect) {
                                mismatchPositions.push_back(static_cast<uint32_t>(gpos));
                            }
                        }
                    }
                } else {
                    if (g1 == 0) { // A
                        ++nT;
                        if (r1 == 2) { // G
                            ++k;
                            if (snpDetect) {
                                mismatchPositions.push_back(static_cast<uint32_t>(gpos));
                            }
                        }
                    }
                }
            }
        }
    }

    if (!isIntronic) {
        if (snpDetect && !mismatchPositions.empty()) {
            for (uint32_t geneId : geneIds) {
                slamQuant->bufferSnpRead(geneId, nT, mismatchPositions, weight);
            }
        } else {
            uint8_t k8 = static_cast<uint8_t>(k > 255 ? 255 : k);
            for (uint32_t geneId : geneIds) {
                slamQuant->addRead(geneId, nT, k8, weight);
            }
        }
    }
    slamQuant->diagnostics().readsProcessed++;
    return true;
}
