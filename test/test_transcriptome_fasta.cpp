/**
 * Transcriptome FASTA Generation - Test Documentation
 * 
 * The TranscriptomeFasta module is tested via STAR integration tests, not a
 * separate external harness. This is because initializing Genome and GTF objects
 * outside STAR's normal pipeline is complex and error-prone.
 * 
 * Running STAR with --genomeGenerateTranscriptome Yes exercises the exact same
 * code path that production uses, eliminating any possibility of code drift.
 * 
 * ## Running Tests
 * 
 *   # Run all tests (synthetic, default-path, cellranger)
 *   ./test/run_transcriptome_generation.sh --all
 * 
 *   # Run specific test modes
 *   ./test/run_transcriptome_generation.sh --synthetic      # Basic tests with small fixtures
 *   ./test/run_transcriptome_generation.sh --default-path   # Test ${genomeDir}/transcriptome.fa
 *   ./test/run_transcriptome_generation.sh --cellranger     # GENCODE chr21+chr22 subset
 *   ./test/run_transcriptome_generation.sh --parity         # Full GENCODE v44 parity test
 * 
 * ## Test Fixtures
 * 
 *   test/fixtures/transcriptome/
 *   ├── small_genome.fa          # Synthetic genome (chr1: 8000bp, chr2: 120bp)
 *   ├── small_genes.gtf          # GTF with 5 test transcripts
 *   ├── expected_transcriptome.fa # Expected output for synthetic test
 *   └── gencode_subset/          # Real GENCODE chr21+chr22 data
 *       ├── genome.fa            # 95MB subset of GRCh38
 *       ├── genes.gtf            # 110k lines from GENCODE v44
 *       └── expected_transcriptome.fa # gffread-generated expected output
 * 
 * ## Test Cases Covered
 * 
 *   1. TR_SINGLE_PLUS   - Single exon, positive strand
 *   2. TR_MULTI_PLUS    - Multiple exons, positive strand (tests concatenation)
 *   3. TR_MULTI_MINUS   - Multiple exons, negative strand (tests reverse complement)
 *   4. TR_LONG          - Long transcript (tests 70bp line wrapping)
 *   5. TR_MISSING       - Exon on non-existent chromosome (tests skipping with warning)
 *   6. CellRanger test  - 7,703 real transcripts with CellRanger biotype filtering
 * 
 * ## Implementation Details
 * 
 *   - Module: source/TranscriptomeFasta.cpp
 *   - Called from: source/Genome_genomeGenerate.cpp after GTF parsing
 *   - Ordering: Matches transcriptInfo.tab (trStart, trEnd, trId) for Salmon parity
 *   - Atomic write: Uses .tmp file + rename() for crash safety
 * 
 * ## Parameters
 * 
 *   --genomeGenerateTranscriptome Yes/No      Enable transcriptome FASTA generation
 *   --genomeGenerateTranscriptomeFasta <path> Custom output path (default: ${genomeDir}/transcriptome.fa)
 *   --genomeGenerateTranscriptomeOverwrite Yes/No  Overwrite existing file
 */

#include <iostream>

int main() {
    std::cout << "Transcriptome FASTA Test Documentation\n";
    std::cout << "=======================================\n\n";
    std::cout << "Tests are run via STAR integration, not this binary.\n\n";
    std::cout << "Run the integration tests:\n\n";
    std::cout << "  # All tests:\n";
    std::cout << "  ./test/run_transcriptome_generation.sh --all\n\n";
    std::cout << "  # Specific modes:\n";
    std::cout << "  ./test/run_transcriptome_generation.sh --synthetic      # Basic tests\n";
    std::cout << "  ./test/run_transcriptome_generation.sh --default-path   # Default path test\n";
    std::cout << "  ./test/run_transcriptome_generation.sh --cellranger     # CellRanger-style test\n";
    std::cout << "  ./test/run_transcriptome_generation.sh --parity         # Full GENCODE v44\n\n";
    std::cout << "The TranscriptomeFasta module is compiled into STAR, so running STAR\n";
    std::cout << "with --genomeGenerateTranscriptome Yes IS the integration test.\n";
    return 0;
}
