#include "TranscriptQuantOutput.h"
#include "em_types.h"
#include "Transcriptome.h"
#include <fstream>
#include <iomanip>
#include <iostream>

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

