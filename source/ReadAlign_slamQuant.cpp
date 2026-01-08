#include "ReadAlign.h"
#include "SlamQuant.h"

bool ReadAlign::slamCollect(const Transcript& trOut, const std::set<uint32_t>& geneIds, double weight) {
    if (slamQuant == nullptr || weight <= 0.0) {
        return false;
    }
    if (geneIds.empty()) {
        slamQuant->diagnostics().readsZeroGenes++;
        return false;
    }

    char* R = Read1[trOut.roStr == 0 ? 0 : 2];
    bool isMinus = (trOut.Str == 1);

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
            uint64 rposRaw = rStart + ii;
            if (readLength[0] > 0 && rposRaw == readLength[0]) {
                continue; // skip spacer between mates
            }
            bool secondMate = (readLength[1] > 0 && rposRaw > readLength[0]);
            uint32_t readPos = static_cast<uint32_t>(secondMate ? rposRaw - 1 : rposRaw);
            slamQuant->addTransitionBase(readPos, secondMate, g1, r1, weight);
            if (!isMinus) {
                if (g1 == 3) { // T
                    ++nT;
                    if (r1 == 1) { // C
                        ++k;
                    }
                }
            } else {
                if (g1 == 0) { // A
                    ++nT;
                    if (r1 == 2) { // G
                        ++k;
                    }
                }
            }
        }
    }

    uint8_t k8 = static_cast<uint8_t>(k > 255 ? 255 : k);
    for (uint32_t geneId : geneIds) {
        slamQuant->addRead(geneId, nT, k8, weight);
    }
    slamQuant->diagnostics().readsProcessed++;
    return true;
}
