#ifndef SAMTOOLS_SORTER_H
#define SAMTOOLS_SORTER_H

#include "IncludeDefine.h"
#include "Parameters.h"
#include <vector>
#include <string>
#include <pthread.h>
#include <cstdint>
#include <climits>
#include <cassert>

// =============================================================================
// SORTING SPECIFICATION (REQUIRED):
// When numeric keys are equal (tid, pos, flag, mtid, mpos, isize), ordering
// MUST be deterministic by QNAME byte comparison. This is a non-negotiable
// requirement for reproducibility and is part of STAR's coordinate sort spec.
// =============================================================================

// SortKey for deterministic coordinate sorting (samtools-equivalent)
// All fields derived from raw BAM buffer without parsing to bam1_t
struct SortKey {
    int32_t tid;        // refID (-1 for unmapped -> INT32_MAX)
    int32_t pos;        // position (-1 for unmapped -> INT32_MAX)
    uint16_t flag;      // BAM flag
    int32_t mtid;       // mate refID
    int32_t mpos;       // mate position
    int32_t isize;      // insert size
    // Note: QNAME comparison done separately for tie-breaking (not stored in key)
    
    // Numeric comparison (QNAME compared separately if all numeric fields equal)
    bool operator<(const SortKey& other) const {
        if (tid != other.tid) return tid < other.tid;
        if (pos != other.pos) return pos < other.pos;
        if (flag != other.flag) return flag < other.flag;
        if (mtid != other.mtid) return mtid < other.mtid;
        if (mpos != other.mpos) return mpos < other.mpos;
        return isize < other.isize;
    }
    
    bool operator==(const SortKey& other) const {
        return tid == other.tid && pos == other.pos && flag == other.flag &&
               mtid == other.mtid && mpos == other.mpos && isize == other.isize;
    }
};

// Structure to hold BAM record data with metadata
struct BAMRecord {
    char* data;
    uint32_t size;
    bool hasY;
    uint64_t iRead;  // readId (bits[63:32]=readId, bits[31:0]=other data, no Y-bit encoding)
    SortKey key;  // Sort key computed from raw BAM buffer
    
    BAMRecord() : data(nullptr), size(0), hasY(false), iRead(0) {
        key.tid = INT32_MAX;
        key.pos = INT32_MAX;
        key.flag = 0;
        key.mtid = -1;
        key.mpos = -1;
        key.isize = 0;
    }
    
    // Move constructor
    BAMRecord(BAMRecord&& other) noexcept 
        : data(other.data), size(other.size), hasY(other.hasY), iRead(other.iRead), key(other.key) {
        other.data = nullptr;
        other.size = 0;
        other.iRead = 0;
    }
    
    // Move assignment
    BAMRecord& operator=(BAMRecord&& other) noexcept {
        if (this != &other) {
            if (data) delete[] data;
            data = other.data;
            size = other.size;
            hasY = other.hasY;
            iRead = other.iRead;
            key = other.key;
            other.data = nullptr;
            other.size = 0;
            other.iRead = 0;
        }
        return *this;
    }
    
    // Delete copy constructor and assignment to prevent accidental copies
    BAMRecord(const BAMRecord& other) = delete;
    BAMRecord& operator=(const BAMRecord& other) = delete;
    
    ~BAMRecord() {
        if (data) delete[] data;
    }
};

// Unified QNAME comparison helper (used by both BAMRecordComparator and HeapLess)
namespace SamtoolsSorterHelpers {
    // REQUIRED: Compare QNAME bytes from raw BAM buffers for tie-breaking.
    // Called when all numeric SortKey fields are equal.
    // This is mandatory for deterministic ordering - not optional.
    inline int compareQName(const char* bamA, const char* bamB) {
        const uint8_t* a8 = reinterpret_cast<const uint8_t*>(bamA + 4);
        const uint8_t* b8 = reinterpret_cast<const uint8_t*>(bamB + 4);
        
        // Extract l_qname from core (offset 8 from start of core = offset 8 from bam8)
        uint8_t l_qname_a = a8[8];
        uint8_t l_qname_b = b8[8];
        
        if (l_qname_a == 0 || l_qname_b == 0) {
            return (l_qname_a < l_qname_b) ? -1 : ((l_qname_a > l_qname_b) ? 1 : 0);
        }
        
        // QNAME starts at offset 36 from bamData (4 block_size + 32 core)
        const char* qname_a = reinterpret_cast<const char*>(bamA + 36);
        const char* qname_b = reinterpret_cast<const char*>(bamB + 36);
        
        // Compare QNAME bytes (excluding trailing NUL if present)
        size_t len_a = l_qname_a;
        size_t len_b = l_qname_b;
        if (len_a > 0 && qname_a[len_a - 1] == 0) len_a--;
        if (len_b > 0 && qname_b[len_b - 1] == 0) len_b--;
        
        size_t min_len = (len_a < len_b) ? len_a : len_b;
        int cmp = memcmp(qname_a, qname_b, min_len);
        if (cmp != 0) return cmp;
        return (len_a < len_b) ? -1 : ((len_a > len_b) ? 1 : 0);
    }
}

// Comparator for sorting BAM records by SortKey + QNAME (samtools-equivalent)
struct BAMRecordComparator {
    bool operator()(const BAMRecord& a, const BAMRecord& b) const {
        // First compare numeric SortKey
        if (a.key < b.key) return true;
        if (b.key < a.key) return false;
        // If all numeric fields equal, compare QNAME bytes
        return SamtoolsSorterHelpers::compareQName(a.data, b.data) < 0;
    }
};

class SamtoolsSorter {
public:
    // Initialize with memory limit and thread count
    SamtoolsSorter(uint64_t maxRAM, int nThreads, const string& tmpDir, Parameters& P);
    
    // Thread-safe: called from coordOneAlign replacement path
    // Takes raw BAM record bytes (same format as coordOneAlign input)
    // iRead: readId (bits[63:32]=readId, bits[31:0]=other data, no Y-bit encoding)
    // hasY: separate Y-chromosome flag (not encoded in iRead)
    void addRecord(const char* bamData, uint32_t bamSize, uint64_t iRead, bool hasY);
    
    // Called after all threads complete; performs sort + merge
    void finalize();
    
    // Iterate sorted records for output (post-sort streaming via k-way merge)
    // Returns false when no more records
    // Caller does NOT own the returned data pointer - it's valid until next call
    // iRead: returns readId (no Y-bit encoding - hasY is separate)
    bool nextRecord(const char** bamData, uint32_t* bamSize, uint64_t* iRead, bool* hasY);
    
    // Cleanup temp files
    ~SamtoolsSorter();

private:
    uint64_t maxRAM_;
    int nThreads_;
    string tmpDir_;
    Parameters& P_;
    
    // Thread-safe buffer
    pthread_mutex_t bufferMutex_;
    std::vector<BAMRecord> records_;
    uint64_t currentRAM_;
    
    // Spill files for when memory limit is exceeded
    std::vector<string> spillFiles_;
    pthread_mutex_t spillFileCounterMutex_;  // Protects spillFileCounter_ from concurrent access
    int spillFileCounter_;
    
    // Forward declaration - SpillFileReader defined in .cpp file
    struct SpillFileReader;
    std::vector<SpillFileReader*> spillReaders_;
    bool finalized_;
    
    // Min-heap for k-way merge with QNAME tie-breaking
    // source_id: -1 = in-memory records, >=0 = spill file index
    struct HeapEntry {
        SortKey key;
        int sourceId;        // -1 = in-memory, >=0 spill file index
        size_t recordIdx;   // index into records_ when sourceId == -1
        const char* bamPtr;  // pointer to current BAM record for tie-break
    };
    
    // Comparator for min-heap with REQUIRED QNAME tie-breaking.
    // INVARIANT: bamPtr must be non-null for any entry in the heap.
    // When numeric keys are equal, QNAME byte comparison determines order.
    struct HeapLess {
        bool operator()(const HeapEntry& a, const HeapEntry& b) const {
            // Compare numeric SortKey first
            if (a.key < b.key) return false;  // a comes before b
            if (b.key < a.key) return true;   // b comes before a
            // Numeric keys equal - QNAME tie-break is REQUIRED
            // INVARIANT: bamPtr must be valid for all heap entries
            assert(a.bamPtr != nullptr && "HeapEntry::bamPtr must not be null");
            assert(b.bamPtr != nullptr && "HeapEntry::bamPtr must not be null");
            return SamtoolsSorterHelpers::compareQName(a.bamPtr, b.bamPtr) > 0;  // min-heap: smaller QNAME has higher priority
        }
    };
    
    std::vector<HeapEntry> mergeHeap_;
    HeapLess heapLess_;  // Comparator instance
    
    // Helper functions
    SortKey computeSortKey(const char* bamData) const;
    void sortAndSpill();
    void initializeKWayMerge();
    void cleanupSpillFiles();
};

#endif // SAMTOOLS_SORTER_H
