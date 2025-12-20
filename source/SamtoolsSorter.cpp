#include "SamtoolsSorter.h"
#include "BAMfunctions.h"
#include "ErrorWarning.h"
#include "Parameters.h"
#include SAMTOOLS_BGZF_H
#include SAMTOOLS_SAM_H
#include <algorithm>
#include <fstream>
#include <cstring>
#include <pthread.h>
#include <queue>
#include <climits>

// SpillFileReader implementation for k-way merge
struct SamtoolsSorter::SpillFileReader {
    ifstream stream;
    BAMRecord currentRecord;
    bool hasRecord;
    string filename;
    int sourceId;
    
    SpillFileReader(const string& fname, int id) : hasRecord(false), filename(fname), sourceId(id) {
        stream.open(fname.c_str(), std::ios::binary);
    }
    
    bool readNext() {
        if (!stream.good()) {
            hasRecord = false;
            return false;
        }
        uint32_t size;
        uint8_t hasYFlag;
        SortKey key;
        
        stream.read(reinterpret_cast<char*>(&size), sizeof(uint32_t));
        if (stream.eof()) {
            hasRecord = false;
            return false;
        }
        
        stream.read(reinterpret_cast<char*>(&hasYFlag), sizeof(uint8_t));
        if (!stream.good()) {
            hasRecord = false;
            return false;
        }
        
        // Read SortKey fields explicitly to avoid padding/ABI issues
        // Note: This is a temp-only format, not a stable ABI
        stream.read(reinterpret_cast<char*>(&key.tid), sizeof(int32_t));
        stream.read(reinterpret_cast<char*>(&key.pos), sizeof(int32_t));
        stream.read(reinterpret_cast<char*>(&key.flag), sizeof(uint16_t));
        stream.read(reinterpret_cast<char*>(&key.mtid), sizeof(int32_t));
        stream.read(reinterpret_cast<char*>(&key.mpos), sizeof(int32_t));
        stream.read(reinterpret_cast<char*>(&key.isize), sizeof(int32_t));
        if (!stream.good()) {
            hasRecord = false;
            return false;
        }
        
        currentRecord.size = size;
        currentRecord.hasY = (hasYFlag != 0);
        currentRecord.key = key;
        if (currentRecord.data) delete[] currentRecord.data;
        currentRecord.data = new char[size];
        stream.read(currentRecord.data, size);
        
        if (stream.good()) {
            hasRecord = true;
            return true;
        } else {
            hasRecord = false;
            return false;
        }
    }
    
    ~SpillFileReader() {
        if (stream.is_open()) stream.close();
    }
};

SamtoolsSorter::SamtoolsSorter(uint64_t maxRAM, int nThreads, const string& tmpDir, Parameters& P)
    : maxRAM_(maxRAM), nThreads_(nThreads), tmpDir_(tmpDir), P_(P),
      currentRAM_(0), spillFileCounter_(0), finalized_(false)
{
    records_.reserve(100000); // Pre-allocate space
    pthread_mutex_init(&bufferMutex_, nullptr);
    pthread_mutex_init(&spillFileCounterMutex_, nullptr);
}

SamtoolsSorter::~SamtoolsSorter() {
    cleanupSpillFiles();
    pthread_mutex_destroy(&bufferMutex_);
    pthread_mutex_destroy(&spillFileCounterMutex_);
}


SortKey SamtoolsSorter::computeSortKey(const char* bamData) const {
    SortKey key;
    
    // Parse BAM core header directly from raw buffer
    // bam32 points to bamData: bam32[0]=block_size, bam32[1]=tid, bam32[2]=pos, ...
    const uint32_t* bam32 = reinterpret_cast<const uint32_t*>(bamData);
    
    // Extract core fields using correct offsets
    // bam32[1]=tid, bam32[2]=pos
    int32_t refID = static_cast<int32_t>(bam32[1]);  // tid
    int32_t pos = static_cast<int32_t>(bam32[2]);    // pos
    
    // Unmapped reads: set to INT32_MAX so they sort last
    if (refID == -1) {
        key.tid = INT32_MAX;
    } else {
        key.tid = refID;
    }
    
    if (pos == -1) {
        key.pos = INT32_MAX;
    } else {
        key.pos = pos;
    }
    
    // Extract flag from bam32[4] (upper 16 bits)
    // bam32[3]=bin_mq_nl (ignore), bam32[4]=flag_nc (flag = high 16 bits)
    key.flag = static_cast<uint16_t>(bam32[4] >> 16);
    
    // Extract mtid, mpos, isize using correct offsets
    // bam32[5]=l_seq (ignore), bam32[6]=mtid, bam32[7]=mpos, bam32[8]=isize
    key.mtid = static_cast<int32_t>(bam32[6]);  // mtid
    key.mpos = static_cast<int32_t>(bam32[7]);  // mpos
    key.isize = static_cast<int32_t>(bam32[8]); // isize
    
    return key;
}

void SamtoolsSorter::addRecord(const char* bamData, uint32_t bamSize, bool hasY) {
    if (bamSize == 0) return;
    
    BAMRecord record;
    record.size = bamSize;
    record.hasY = hasY;
    record.data = new char[bamSize];
    memcpy(record.data, bamData, bamSize);
    
    // Compute sort key from raw buffer
    record.key = computeSortKey(bamData);
    
    // Quick check without lock for common case
    bool needSpill = false;
    {
        pthread_mutex_lock(&bufferMutex_);
        records_.push_back(std::move(record));
        // Estimate RAM usage: record data + overhead
        currentRAM_ += bamSize + sizeof(BAMRecord) + 64; // 64 bytes overhead per record
        
        // Check if we need to spill to disk
        if (currentRAM_ > maxRAM_ && maxRAM_ > 0 && !records_.empty()) {
            needSpill = true;
        }
        pthread_mutex_unlock(&bufferMutex_);
    }
    
    // Spill outside the lock to avoid blocking writers
    if (needSpill) {
        sortAndSpill();
    }
}

void SamtoolsSorter::sortAndSpill() {
    // Lock only to extract records, then release for sorting/writing
    std::vector<BAMRecord> recordsToSpill;
    {
        pthread_mutex_lock(&bufferMutex_);
        if (records_.empty()) {
            pthread_mutex_unlock(&bufferMutex_);
            return;
        }
        // Move all records out (no copy due to move semantics)
        recordsToSpill = std::move(records_);
        records_.clear();
        currentRAM_ = 0;
        pthread_mutex_unlock(&bufferMutex_);
    }
    
    // Sort records by SortKey + QNAME (outside lock)
    std::sort(recordsToSpill.begin(), recordsToSpill.end(), BAMRecordComparator());
    
    // Write sorted chunk to temp file with SortKey fields serialized explicitly
    // Protect spillFileCounter_ from concurrent access
    int spillFileNum;
    {
        pthread_mutex_lock(&spillFileCounterMutex_);
        spillFileNum = spillFileCounter_++;
        pthread_mutex_unlock(&spillFileCounterMutex_);
    }
    string spillFile = tmpDir_ + "/samtools_sort_spill_" + to_string(spillFileNum) + ".dat";
    
    ofstream spillStream(spillFile.c_str(), std::ios::binary);
    if (!spillStream) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: could not open spill file: " << spillFile << "\n";
        errOut << "SOLUTION: check that the disk is not full, increase the max number of open files with Linux command ulimit -n before running STAR";
        exitWithError(errOut.str(), std::cerr, P_.inOut->logMain, EXIT_CODE_PARAMETER, P_);
    }
    
    // Write records to spill file with metadata
    // Format: [bamSize:uint32][hasY:uint8][key.tid:int32][key.pos:int32][key.flag:uint16][key.mtid:int32][key.mpos:int32][key.isize:int32][bamData:bytes]
    // Note: Serializing fields explicitly to avoid struct padding/ABI issues (temp-only format)
    for (auto& record : recordsToSpill) {
        uint32_t size = record.size;
        uint8_t hasYFlag = record.hasY ? 1 : 0;
        spillStream.write(reinterpret_cast<const char*>(&size), sizeof(uint32_t));
        spillStream.write(reinterpret_cast<const char*>(&hasYFlag), sizeof(uint8_t));
        spillStream.write(reinterpret_cast<const char*>(&record.key.tid), sizeof(int32_t));
        spillStream.write(reinterpret_cast<const char*>(&record.key.pos), sizeof(int32_t));
        spillStream.write(reinterpret_cast<const char*>(&record.key.flag), sizeof(uint16_t));
        spillStream.write(reinterpret_cast<const char*>(&record.key.mtid), sizeof(int32_t));
        spillStream.write(reinterpret_cast<const char*>(&record.key.mpos), sizeof(int32_t));
        spillStream.write(reinterpret_cast<const char*>(&record.key.isize), sizeof(int32_t));
        spillStream.write(record.data, record.size);
    }
    
    spillStream.close();
    
    // Store spill file name (will be used for k-way merge)
    // Use spillFileCounterMutex_ since we're already protecting counter access
    {
        pthread_mutex_lock(&spillFileCounterMutex_);
        spillFiles_.push_back(spillFile);
        pthread_mutex_unlock(&spillFileCounterMutex_);
    }
}

void SamtoolsSorter::finalize() {
    if (finalized_) return;
    
    pthread_mutex_lock(&bufferMutex_);
    
    // Sort remaining in-memory records
    if (!records_.empty()) {
        std::sort(records_.begin(), records_.end(), BAMRecordComparator());
    }
    
    // Initialize k-way merge with spill files
    initializeKWayMerge();
    
    finalized_ = true;
    pthread_mutex_unlock(&bufferMutex_);
}

void SamtoolsSorter::initializeKWayMerge() {
    mergeHeap_.clear();
    
    // Add in-memory records to heap (if any)
    if (!records_.empty()) {
        HeapEntry entry;
        entry.key = records_[0].key;
        entry.sourceId = -1;  // -1 indicates in-memory source
        entry.recordIdx = 0;   // First record in sorted vector
        entry.bamPtr = records_[0].data;  // Pointer to BAM buffer for QNAME comparison
        mergeHeap_.push_back(entry);
    }
    
    // Open all spill files for reading and add first record from each to heap
    for (size_t i = 0; i < spillFiles_.size(); i++) {
        SpillFileReader* reader = new SpillFileReader(spillFiles_[i], static_cast<int>(i));
        if (reader->readNext()) {
            HeapEntry entry;
            entry.key = reader->currentRecord.key;
            entry.sourceId = reader->sourceId;
            entry.recordIdx = 0;  // Unused for spill files
            entry.bamPtr = reader->currentRecord.data;  // Pointer to BAM buffer for QNAME comparison
            mergeHeap_.push_back(entry);
            spillReaders_.push_back(reader);
        } else {
            delete reader;
        }
    }
    
    // Build min-heap using explicit comparator with QNAME tie-breaking
    std::make_heap(mergeHeap_.begin(), mergeHeap_.end(), heapLess_);
}

bool SamtoolsSorter::nextRecord(const char** bamData, uint32_t* bamSize, bool* hasY) {
    if (!finalized_) {
        finalize();
    }
    
    // K-way merge using min-heap with QNAME tie-breaking in comparator
    while (!mergeHeap_.empty()) {
        // Get smallest key from heap (QNAME tie-breaking handled by comparator)
        std::pop_heap(mergeHeap_.begin(), mergeHeap_.end(), heapLess_);
        HeapEntry top = mergeHeap_.back();
        mergeHeap_.pop_back();
        
        // Get the record for this entry
        BAMRecord* record = nullptr;
        if (top.sourceId == -1) {
            // In-memory source - use recordIdx to get the correct record
            if (top.recordIdx >= records_.size()) {
                continue;  // Invalid index, try next heap entry
            }
            record = &records_[top.recordIdx];
        } else {
            // Spill file source
            if (top.sourceId >= static_cast<int>(spillReaders_.size())) {
                continue;  // Invalid source, try next
            }
            SpillFileReader* reader = spillReaders_[top.sourceId];
            if (!reader->hasRecord) {
                continue;  // Reader exhausted, try next
            }
            record = &reader->currentRecord;
        }
        
        if (record == nullptr) {
            continue;
        }
        
        // Advance the source for the record we're returning
        if (top.sourceId == -1) {
            // Advance to next in-memory record
            size_t nextIdx = top.recordIdx + 1;
            
            // Push next in-memory record if available
            if (nextIdx < records_.size()) {
                HeapEntry nextEntry;
                nextEntry.key = records_[nextIdx].key;
                nextEntry.sourceId = -1;
                nextEntry.recordIdx = nextIdx;
                nextEntry.bamPtr = records_[nextIdx].data;  // Update bamPtr
                mergeHeap_.push_back(nextEntry);
                std::push_heap(mergeHeap_.begin(), mergeHeap_.end(), heapLess_);
            }
        } else {
            // Read next record from this spill file
            SpillFileReader* reader = spillReaders_[top.sourceId];
            if (reader->readNext()) {
                // Push next record from this source
                HeapEntry nextEntry;
                nextEntry.key = reader->currentRecord.key;
                nextEntry.sourceId = reader->sourceId;
                nextEntry.recordIdx = 0;  // Unused for spill files
                nextEntry.bamPtr = reader->currentRecord.data;  // Update bamPtr
                mergeHeap_.push_back(nextEntry);
                std::push_heap(mergeHeap_.begin(), mergeHeap_.end(), heapLess_);
            } else {
                // File exhausted, will be cleaned up later
            }
        }
        
        // Return the record
        *bamData = record->data;
        *bamSize = record->size;
        *hasY = record->hasY;
        return true;
    }
    
    return false;
}

void SamtoolsSorter::cleanupSpillFiles() {
    // Close and delete all spill file readers
    for (auto* reader : spillReaders_) {
        delete reader;
    }
    spillReaders_.clear();
    
    // Delete spill files
    for (const string& spillFile : spillFiles_) {
        remove(spillFile.c_str());
    }
    spillFiles_.clear();
}
