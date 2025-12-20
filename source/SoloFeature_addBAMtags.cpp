#include "SoloFeature.h"
#include "SequenceFuns.h"
#include "BAMfunctions.h"
#include "ErrorWarning.h"
#include <cstring>

// Legacy overload: extracts iRead from trailing 8 bytes (legacy bin sorter convention with Y-bit encoding)
void SoloFeature::addBAMtags(char *&bam0, uint32 &size0, char *bam1) {
    // Extract iReadWithY from trailing 8 bytes (legacy convention: bit63=Y-flag)
    uint64_t iReadWithY = *(uint64_t*)(bam0 + size0);
    // Mask Y-bit and extract readId for compatibility with samtools path
    uint64_t iRead = (iReadWithY & ~(1ULL << 63));
    addBAMtags(bam0, size0, bam1, iRead);
}

// Main implementation: explicit iRead for samtools sorter (no Y-bit encoding)
//
// Policy: CB+UB together (requireCbUbTogether)
// - When true (default): UB requires CB, CB lookup failures are fatal if UB is requested.
//   This is an intentional policy to ensure data consistency and may be relaxed later.
// - When false: Allow UB-only or CB-only injection without hard-failure.
//   CB lookup is only guarded when needCB is true, and UB decode failures are non-fatal.
//
// To relax the policy: set requireCbUbTogether = false below.
void SoloFeature::addBAMtags(char *&bam0, uint32 &size0, char *bam1, uint64_t iRead) {
    // Policy flag: require CB+UB together (default: true)
    // Set to false to allow UB-only or CB-only injection
    static const bool requireCbUbTogether = true;
    
    // Early return if tags not requested
    if (!pSolo.samAttrYes) {
        return;
    }
    
    // Check if CB or UB are requested
    bool needCB = P.outSAMattrPresent.CB;
    bool needUB = P.outSAMattrPresent.UB;
    if (!needCB && !needUB) {
        return;
    }
    
    // Policy enforcement: if requireCbUbTogether and UB is requested, require CB
    if (requireCbUbTogether && needUB && !needCB) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: UB tag requested but CB tag not requested\n";
        errOut << "SOLUTION: Request both CB and UB tags, or set requireCbUbTogether = false in SoloFeature_addBAMtags.cpp";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        return;
    }
    
    // Hard error: Check packedReadInfo is initialized
    if (packedReadInfo.data.empty()) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: PackedReadInfo not initialized for CB/UB injection\n";
        errOut << "SOLUTION: This should not happen - please report this error";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        return;
    }
    
    // Extract readId (samtools path: no Y-bit encoding, just shift)
    // iRead format: bits[63:32]=readId, bits[31:0]=other data
    uint32_t readId = static_cast<uint32_t>(iRead >> 32);
    
    // Hard error: Check readId in range
    if (readId >= packedReadInfo.data.size()) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: readId out of range for PackedReadInfo\n";
        errOut << "readId=" << readId << " packedReadInfo.data.size()=" << packedReadInfo.data.size() << "\n";
        errOut << "SOLUTION: This should not happen - please report this error";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        return;
    }
    
    // Lookup from PackedReadInfo
    uint32_t cbIdx = packedReadInfo.getCB(readId);
    uint32_t umiPacked = packedReadInfo.getUMI(readId);
    uint8_t status = packedReadInfo.getStatus(readId);
    
    // Status-aware tag emission (matches recordReadInfo semantics)
    // status==0: skip both CB and UB (missing CB) - no error
    // status==1: emit both CB and UB (valid)
    // status==2: emit CB only, skip UB (invalid UMI)
    if (status == 0) {
        // Policy: if requireCbUbTogether and UB is requested, missing CB is fatal
        if (requireCbUbTogether && needUB) {
            ostringstream errOut;
            errOut << "EXITING because of fatal ERROR: UB tag requested but CB is missing (status==0)\n";
            errOut << "SOLUTION: This is enforced by requireCbUbTogether policy. Set requireCbUbTogether = false to allow UB-only injection.";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
            return;
        }
        return;  // No tags to add
    }
    
    // Decode CB with bounds check
    // Policy: if requireCbUbTogether and UB is requested, CB lookup failures are fatal
    // Otherwise, only guard CB lookup when needCB is true
    bool cbLookupRequired = needCB || (requireCbUbTogether && needUB);
    if (status != 0 && cbLookupRequired && cbIdx >= pSolo.cbWLstr.size()) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: cbIdx out of range for whitelist\n";
        errOut << "cbIdx=" << cbIdx << " pSolo.cbWLstr.size()=" << pSolo.cbWLstr.size() << "\n";
        if (requireCbUbTogether && needUB && !needCB) {
            errOut << "NOTE: CB lookup required because UB is requested and requireCbUbTogether=true\n";
        }
        errOut << "SOLUTION: This should not happen - please report this error";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        return;
    }
    
    string cbStr;
    bool emitCB = false;
    if (needCB && status != 0 && cbIdx < pSolo.cbWLstr.size()) {
        cbStr = pSolo.cbWLstr[cbIdx];
        emitCB = true;
    }
    
    // Decode UMI using length-aware decoder (only if status == 1)
    // Policy: if requireCbUbTogether is false, allow UB decode failures to be non-fatal
    string ubStr;
    bool emitUB = false;
    if (status == 1 && needUB && pSolo.umiL > 0 && pSolo.umiL <= 16) {
        // Use convertNuclInt64toString from SequenceFuns.h
        // Do NOT skip when umiPacked == 0 (valid all-A UMI)
        ubStr = convertNuclInt64toString(umiPacked, pSolo.umiL);
        emitUB = true;
    }
    
    // If no tags to add, leave bam0/size0 unchanged
    if (!emitCB && !emitUB) {
        return;
    }
    
    // Block_size sanity check: verify size0 matches block_size + 4
    uint32_t blockSizeFromHeader = *(uint32_t*)bam0;
    uint32_t expectedSize = blockSizeFromHeader + 4;
    if (size0 != expectedSize) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: BAM record size mismatch\n";
        errOut << "size0=" << size0 << " expected (block_size+4)=" << expectedSize << " block_size=" << blockSizeFromHeader << "\n";
        errOut << "SOLUTION: This indicates corrupted BAM data - please report this error";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        return;
    }
    
    // Manual raw tag append (no htslib bam_aux_append - unsafe on non-owned memory)
    // BAM format: [block_size:uint32][core:32 bytes][data:variable]
    // block_size = size of core + data (doesn't include the 4-byte size field itself)
    // data includes: QNAME, CIGAR, SEQ, QUAL, aux tags
    
    // Read current block_size (already validated above)
    uint32_t blockSize = blockSizeFromHeader;
    
    // Calculate tag sizes
    uint32_t tagSize = 0;
    if (emitCB) {
        tagSize += 3 + cbStr.size() + 1;  // CB:Z:<cbStr>\0
    }
    if (emitUB) {
        tagSize += 3 + ubStr.size() + 1;  // UB:Z:<ubStr>\0
    }
    
    // Calculate new sizes
    uint32_t newBlockSize = blockSize + tagSize;
    uint32_t newSize = newBlockSize + 4;  // +4 for block_size field
    
    // Hard error on size overflow
    if (newSize > BAM_ATTR_MaxSize) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: BAM record exceeds maximum size after tag injection\n";
        errOut << "newSize=" << newSize << " BAM_ATTR_MaxSize=" << BAM_ATTR_MaxSize << "\n";
        errOut << "SOLUTION: This should not happen - please report this error";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        return;
    }
    
    // Copy existing record to bam1
    memcpy(bam1, bam0, size0);
    
    // Update block_size
    uint32_t* blockSizePtr = reinterpret_cast<uint32_t*>(bam1);
    *blockSizePtr = newBlockSize;
    
    // Append tags at end of existing data (after block_size + core + existing data)
    // Current data ends at: bam0 + 4 + blockSize
    // New data ends at: bam1 + 4 + newBlockSize
    char* tagPtr = bam1 + 4 + blockSize;
    
    if (emitCB) {
        // Append CB tag: CB:Z:<cbStr>\0
        tagPtr[0] = 'C';
        tagPtr[1] = 'B';
        tagPtr[2] = 'Z';
        memcpy(tagPtr + 3, cbStr.c_str(), cbStr.size() + 1);
        tagPtr += 3 + cbStr.size() + 1;
    }
    
    if (emitUB) {
        // Append UB tag: UB:Z:<ubStr>\0
        tagPtr[0] = 'U';
        tagPtr[1] = 'B';
        tagPtr[2] = 'Z';
        memcpy(tagPtr + 3, ubStr.c_str(), ubStr.size() + 1);
        tagPtr += 3 + ubStr.size() + 1;
    }
    
    // Update bam0/size0 to point to new record
    bam0 = bam1;
    size0 = newSize;
}

