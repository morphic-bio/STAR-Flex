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
    double localCoverage[4] = {0.0, 0.0, 0.0, 0.0};
    double localMismatches[4][4] = {{0.0}};
    bool hasSnp = false;

    for (uint iex = 0; iex < trOut.nExons; ++iex) {
        uint64 gStart = trOut.exons[iex][EX_G];
        uint64 rStart = trOut.exons[iex][EX_R];
        uint64 len = trOut.exons[iex][EX_L];

        for (uint64 ii = 0; ii < len; ++ii) {
            uint64 gpos = gStart + ii;
            if (slamSnpMask != nullptr && slamSnpMask->contains(gpos)) {
                hasSnp = true;
                break; // discard read if it overlaps SNP mask
            }
            char r1 = R[rStart + ii];
            char g1 = genOut.G[gpos];
            if (r1 == 4 || g1 == 4) {
                continue;
            }
            if (g1 >= 0 && g1 < 4 && r1 >= 0 && r1 < 4) {
                localCoverage[static_cast<int>(g1)] += 1.0;
                if (g1 != r1) {
                    localMismatches[static_cast<int>(g1)][static_cast<int>(r1)] += 1.0;
                }
            }
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
        if (hasSnp) {
            break;
        }
    }

    if (hasSnp) {
        slamQuant->diagnostics().readsDroppedSnpMask++;
        return false;
    }

    uint8_t k8 = static_cast<uint8_t>(k > 255 ? 255 : k);
    for (uint32_t geneId : geneIds) {
        slamQuant->addRead(geneId, nT, k8, weight);
    }
    slamQuant->addTransitions(localCoverage, localMismatches, weight);
    slamQuant->diagnostics().readsProcessed++;
    return true;
}
