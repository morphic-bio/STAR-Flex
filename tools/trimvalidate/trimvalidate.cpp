/**
 * trimvalidate: Standalone CLI tool for cutadapt-parity trimming validation
 * 
 * Reads paired-end FASTQ files, applies trimming using libtrim, and outputs
 * trimmed FASTQ files for comparison with Trim Galore/cutadapt.
 * 
 * Usage:
 *   ./trimvalidate -1 input_R1.fastq -2 input_R2.fastq \
 *                  -o1 trimmed_R1.fastq -o2 trimmed_R2.fastq \
 *                  [--quality 20] [--length 20] \
 *                  [--adapter-r1 SEQ] [--adapter-r2 SEQ]
 */

#include "../../source/libtrim/trim.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cstdlib>

void printUsage(const char* progName) {
    std::cerr << "Usage: " << progName << " -1 <R1.fastq> -2 <R2.fastq> "
              << "-o1 <out_R1.fastq> -o2 <out_R2.fastq> [options]\n";
    std::cerr << "\nOptions:\n";
    std::cerr << "  --quality <int>     Quality threshold (default: 20)\n";
    std::cerr << "  --length <int>      Minimum read length (default: 20)\n";
    std::cerr << "  --adapter-r1 <seq>  R1 adapter sequence (default: TruSeq R1)\n";
    std::cerr << "  --adapter-r2 <seq>  R2 adapter sequence (default: TruSeq R2)\n";
    std::cerr << "\n";
}

struct FastqRecord {
    std::string name;
    std::string seq;
    std::string qual;
};

bool readFastqRecord(std::ifstream& in, FastqRecord& rec) {
    if (!std::getline(in, rec.name)) return false;
    if (rec.name.empty() || rec.name[0] != '@') return false;
    
    if (!std::getline(in, rec.seq)) return false;
    
    std::string plus;
    if (!std::getline(in, plus)) return false;
    if (plus[0] != '+') return false;
    
    if (!std::getline(in, rec.qual)) return false;
    
    return true;
}

void writeFastqRecord(std::ofstream& out, const FastqRecord& rec) {
    out << rec.name << "\n";
    out << rec.seq << "\n";
    out << "+\n";
    out << rec.qual << "\n";
}

int main(int argc, char** argv) {
    std::string r1_in, r2_in, r1_out, r2_out;
    struct TrimParams params;
    trim_params_init(&params);
    
    // Parse arguments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-1") == 0 && i + 1 < argc) {
            r1_in = argv[++i];
        } else if (strcmp(argv[i], "-2") == 0 && i + 1 < argc) {
            r2_in = argv[++i];
        } else if (strcmp(argv[i], "-o1") == 0 && i + 1 < argc) {
            r1_out = argv[++i];
        } else if (strcmp(argv[i], "-o2") == 0 && i + 1 < argc) {
            r2_out = argv[++i];
        } else if (strcmp(argv[i], "--quality") == 0 && i + 1 < argc) {
            params.quality_cutoff = (uint8_t)atoi(argv[++i]);
        } else if (strcmp(argv[i], "--length") == 0 && i + 1 < argc) {
            params.min_length = (uint32_t)atoi(argv[++i]);
        } else if (strcmp(argv[i], "--adapter-r1") == 0 && i + 1 < argc) {
            params.adapter_r1 = argv[++i];
        } else if (strcmp(argv[i], "--adapter-r2") == 0 && i + 1 < argc) {
            params.adapter_r2 = argv[++i];
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printUsage(argv[0]);
            return 0;
        } else {
            std::cerr << "Unknown option: " << argv[i] << "\n";
            printUsage(argv[0]);
            return 1;
        }
    }
    
    if (r1_in.empty() || r2_in.empty() || r1_out.empty() || r2_out.empty()) {
        std::cerr << "Error: Missing required arguments\n";
        printUsage(argv[0]);
        return 1;
    }
    
    // Open files
    std::ifstream in1(r1_in);
    std::ifstream in2(r2_in);
    std::ofstream out1(r1_out);
    std::ofstream out2(r2_out);
    
    if (!in1.is_open() || !in2.is_open()) {
        std::cerr << "Error: Cannot open input files\n";
        return 1;
    }
    if (!out1.is_open() || !out2.is_open()) {
        std::cerr << "Error: Cannot open output files\n";
        return 1;
    }
    
    // Process pairs
    uint64_t pairs_processed = 0;
    uint64_t pairs_dropped = 0;
    
    FastqRecord rec1, rec2;
    while (readFastqRecord(in1, rec1) && readFastqRecord(in2, rec2)) {
        pairs_processed++;
        
        // Convert to char buffers (with space for trimming)
        const size_t max_len = 1000;  // Reasonable max read length
        char seq1[max_len + 1], qual1[max_len + 1];
        char seq2[max_len + 1], qual2[max_len + 1];
        
        uint32_t len1 = rec1.seq.length();
        uint32_t len2 = rec2.seq.length();
        
        if (len1 > max_len || len2 > max_len) {
            std::cerr << "Warning: Read length exceeds " << max_len << ", skipping\n";
            continue;
        }
        
        memcpy(seq1, rec1.seq.c_str(), len1);
        memcpy(qual1, rec1.qual.c_str(), len1);
        memcpy(seq2, rec2.seq.c_str(), len2);
        memcpy(qual2, rec2.qual.c_str(), len2);
        
        // Trim
        struct TrimResult result1, result2;
        trim_pair(seq1, qual1, &len1, seq2, qual2, &len2, &params, &result1, &result2);
        
        // Write if not dropped
        if (!result1.dropped && !result2.dropped) {
            // Update records with trimmed sequences
            rec1.seq.assign(seq1, len1);
            rec1.qual.assign(qual1, len1);
            rec2.seq.assign(seq2, len2);
            rec2.qual.assign(qual2, len2);
            
            writeFastqRecord(out1, rec1);
            writeFastqRecord(out2, rec2);
        } else {
            pairs_dropped++;
        }
    }
    
    std::cerr << "Processed " << pairs_processed << " pairs\n";
    std::cerr << "Dropped " << pairs_dropped << " pairs (below min length)\n";
    
    return 0;
}
