#ifndef SLAM_DUMP_H
#define SLAM_DUMP_H

#include "SlamReadBuffer.h"
#include <cstdint>
#include <string>
#include <vector>

struct SlamDumpMetadata {
    uint32_t version = 1;
    uint64_t nReads = 0;
    double errorRate = 0.0;
    double convRate = 0.0;
    std::vector<std::string> geneIds;
    std::vector<std::string> geneNames;
    std::vector<std::string> chrNames;
    std::vector<uint64_t> chrStart;
};

// Write dump file from one or more buffers.
// maxReads=0 means no limit.
bool writeSlamDump(const std::string& path,
                   const SlamDumpMetadata& meta,
                   const std::vector<const SlamReadBuffer*>& buffers,
                   uint64_t maxReads,
                   std::string* err);

// Read dump file into metadata + reads.
bool readSlamDump(const std::string& path,
                  SlamDumpMetadata* meta,
                  std::vector<SlamBufferedRead>* reads,
                  std::string* err);

#endif // SLAM_DUMP_H
