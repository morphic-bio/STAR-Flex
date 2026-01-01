#include "TranscriptQuantOutput.h"
#include "em_types.h"
#include "Transcriptome.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cmath>

// Helper: format double to string, trimming trailing zeros and decimal point
// Matches Salmon quant.genes.sf format: "510" not "510.000000", "295.805" not "295.805000"
// Note: Values below 1e-maxPrecision may round to 0 due to fixed precision + trimming
static std::string formatNum(double val, int maxPrecision = 6) {
    // Handle exact zeros
    if (val == 0.0) {
        return "0";
    }
    
    // Use stringstream with fixed precision, then trim
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(maxPrecision) << val;
    std::string s = oss.str();
    
    // Find decimal point
    size_t dot = s.find('.');
    if (dot != std::string::npos) {
        // Trim trailing zeros
        size_t last = s.find_last_not_of('0');
        if (last != std::string::npos && last > dot) {
            s = s.substr(0, last + 1);
        } else if (last == dot) {
            // All zeros after dot, remove the dot too
            s = s.substr(0, dot);
        }
        // Remove trailing dot if present
        if (!s.empty() && s.back() == '.') {
            s.pop_back();
        }
    }
    
    return s;
}

void writeQuantSF(const EMResult& result, const TranscriptState& state, const std::string& filename) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Error: Failed to open output file: " << filename << std::endl;
        return;
    }
    
    // Write header (Salmon-compatible format)
    out << "Name\tLength\tEffectiveLength\tTPM\tNumReads\n";
    
    // Write transcript data
    for (size_t i = 0; i < state.n; ++i) {
        out << state.names[i] << "\t"
            << std::fixed << std::setprecision(0) << state.lengths[i] << "\t"
            << std::fixed << std::setprecision(3) << state.eff_lengths[i] << "\t"
            << std::fixed << std::setprecision(6) << result.tpm[i] << "\t"
            << std::fixed << std::setprecision(3) << result.counts[i] << "\n";
    }
    
    out.close();
}

// Return codes: 0 = success, 1 = file error, 2 = had MissingGeneID
int writeQuantGeneSF(const EMResult& result, const TranscriptState& state,
                     const Transcriptome& tr, const std::string& filename) {
    // Allocate gene-level accumulators
    std::vector<double> geneTPM(tr.nGe, 0.0);
    std::vector<double> geneCounts(tr.nGe, 0.0);
    std::vector<double> geneLenWeightedTPM(tr.nGe, 0.0);     // TPM-weighted sum
    std::vector<double> geneEffLenWeightedTPM(tr.nGe, 0.0);  // TPM-weighted sum
    std::vector<double> geneLenUnweighted(tr.nGe, 0.0);     // Uniform sum (fallback when TPM == 0)
    std::vector<double> geneEffLenUnweighted(tr.nGe, 0.0);   // Uniform sum (fallback when TPM == 0)
    std::vector<uint32_t> geneTrCount(tr.nGe, 0);           // Transcript count per gene

    // Aggregate transcripts to genes
    for (size_t t = 0; t < state.n; ++t) {
        uint32_t g = tr.trGene[t];
        double tpm = result.tpm[t];
        double cnt = result.counts[t];
        
        geneTPM[g] += tpm;
        geneCounts[g] += cnt;
        
        // TPM-weighted sums (used when geneTPM > 0)
        geneLenWeightedTPM[g] += tpm * state.lengths[t];
        geneEffLenWeightedTPM[g] += tpm * state.eff_lengths[t];
        
        // Unweighted sums (used when geneTPM == 0, matching Salmon's uniform averaging)
        geneLenUnweighted[g] += state.lengths[t];
        geneEffLenUnweighted[g] += state.eff_lengths[t];
        
        geneTrCount[g]++;
    }

    // Open output file
    std::ofstream out(filename);
    if (!out.is_open()) {
        return 1;  // File error - caller logs via logMain
    }

    // Write header
    out << "Name\tLength\tEffectiveLength\tTPM\tNumReads\n";

    bool hasMissingGeneID = false;
    
    for (uint32_t g = 0; g < tr.nGe; ++g) {
        double len, effLen;

        if (geneTPM[g] > 0) {
            // Primary: TPM-weighted average for genes with TPM > 0 (matches Salmon)
            len = geneLenWeightedTPM[g] / geneTPM[g];
            effLen = geneEffLenWeightedTPM[g] / geneTPM[g];
        } else if (geneTrCount[g] > 0) {
            // Fallback: Uniform average when geneTPM == 0 (matches Salmon's aggregateEstimatesToGeneLevel)
            // This replaces the previous count-weighted fallback to match Salmon behavior
            len = geneLenUnweighted[g] / geneTrCount[g];
            effLen = geneEffLenUnweighted[g] / geneTrCount[g];
        } else {
            // Edge case: gene with no transcripts (shouldn't happen)
            len = 0.0;
            effLen = 0.0;
        }

        // Track MissingGeneID for caller to log
        if (tr.geID[g] == "MissingGeneID") {
            hasMissingGeneID = true;
        }

        // Output with trailing-zero trimming to match Salmon format
        // Per-column precision: Length(3), EffectiveLength(4), TPM(6), NumReads(3)
        out << tr.geID[g] << "\t"
            << formatNum(len, 3) << "\t"
            << formatNum(effLen, 4) << "\t"
            << formatNum(geneTPM[g], 6) << "\t"
            << formatNum(geneCounts[g], 3) << "\n";
    }

    out.close();
    return hasMissingGeneID ? 2 : 0;
}
