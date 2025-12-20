#include "SamtoolsSorter.h"
#include "Parameters.h"
#include <sys/stat.h>
#include <cstring>
#include <vector>
#include <string>
#include <iostream>
#include <dirent.h>
#include <algorithm>
#include <cstdlib>
#include <unistd.h>

static std::vector<char> make_bam_record(int32_t tid, int32_t pos,
                                         const std::string& qname,
                                         uint16_t flag,
                                         int32_t mtid, int32_t mpos, int32_t isize) {
    uint8_t l_qname = static_cast<uint8_t>(qname.size() + 1); // include NUL
    uint16_t n_cigar = 1;
    int32_t l_seq = 1;

    uint32_t bin = 0;
    uint8_t mapq = 0;
    uint32_t bin_mq_nl = (bin << 16) | (mapq << 8) | l_qname;
    uint32_t flag_nc = (static_cast<uint32_t>(flag) << 16) | n_cigar;

    uint32_t data_len = 32 + l_qname + 4 * n_cigar + ((l_seq + 1) / 2) + l_seq;
    std::vector<char> buf(4 + data_len, 0);

    uint32_t* p32 = reinterpret_cast<uint32_t*>(buf.data());
    p32[0] = data_len;
    p32[1] = static_cast<uint32_t>(tid);
    p32[2] = static_cast<uint32_t>(pos);
    p32[3] = bin_mq_nl;
    p32[4] = flag_nc;
    p32[5] = static_cast<uint32_t>(l_seq);
    p32[6] = static_cast<uint32_t>(mtid);
    p32[7] = static_cast<uint32_t>(mpos);
    p32[8] = static_cast<uint32_t>(isize);

    char* p = buf.data() + 36; // 4 + 32
    memcpy(p, qname.c_str(), qname.size());
    p[qname.size()] = '\0';
    p += l_qname;

    uint32_t cigar = (1u << 4) | 0u; // 1M
    memcpy(p, &cigar, 4);
    p += 4;

    p[0] = 1 << 4; // seq 'A'
    p += 1;

    p[0] = 30; // qual
    return buf;
}

static std::string get_qname(const char* bam) {
    const uint8_t* core = reinterpret_cast<const uint8_t*>(bam + 4);
    uint8_t l_qname = core[8];
    const char* q = bam + 36;
    if (l_qname > 0 && q[l_qname - 1] == 0) l_qname--;
    return std::string(q, q + l_qname);
}

int main() {
    Parameters P;
    // Use unique temp dir per test run
    P.outBAMsortTmpDir = "/tmp/samtools_sorter_test_" + std::to_string(getpid());
    
    // Clean up any existing test directory
    system(("rm -rf " + P.outBAMsortTmpDir).c_str());
    mkdir(P.outBAMsortTmpDir.c_str(), 0755);

    // Force spills with tiny RAM
    SamtoolsSorter sorter(64, 1, P.outBAMsortTmpDir, P);

    // Same numeric key, different QNAME
    auto r1 = make_bam_record(0, 10, "readB", 0, -1, -1, 0);
    auto r2 = make_bam_record(0, 10, "readA", 0, -1, -1, 0);
    // Unmapped (should be last)
    auto r3 = make_bam_record(-1, -1, "zzz", 0, -1, -1, 0);

    sorter.addRecord(r1.data(), r1.size(), 0ULL, false);
    sorter.addRecord(r2.data(), r2.size(), 1ULL, false);
    sorter.addRecord(r3.data(), r3.size(), 2ULL, false);
    
    // Check if spill files were created (with tiny maxRAM, should trigger)
    DIR* dir = opendir(P.outBAMsortTmpDir.c_str());
    int spillCount = 0;
    if (dir != nullptr) {
        struct dirent* entry;
        while ((entry = readdir(dir)) != nullptr) {
            std::string name(entry->d_name);
            if (name.find("samtools_sort_spill_") == 0 && name.find(".dat") != std::string::npos) {
                spillCount++;
            }
        }
        closedir(dir);
    }
    
    sorter.finalize();

    std::vector<std::string> out;
    std::vector<uint64_t> iReads;
    const char* bam = nullptr;
    uint32_t size = 0;
    uint64_t iRead = 0;
    bool hasY = false;
    while (sorter.nextRecord(&bam, &size, &iRead, &hasY)) {
        out.push_back(get_qname(bam));
        iReads.push_back(iRead);
    }

    if (out.size() != 3 ||
        out[0] != "readA" ||
        out[1] != "readB" ||
        out[2] != "zzz") {
        std::cerr << "FAIL: unexpected order\n";
        return 1;
    }
    
    // Verify iRead round-trips through spill/merge unchanged
    if (iReads.size() != 3) {
        std::cerr << "FAIL: expected 3 iRead values, got " << iReads.size() << "\n";
        return 1;
    }
    // Check that iRead values match what we put in (order may differ due to sorting)
    std::vector<uint64_t> expectedReads = {0ULL, 1ULL, 2ULL};
    std::sort(iReads.begin(), iReads.end());
    std::sort(expectedReads.begin(), expectedReads.end());
    if (iReads != expectedReads) {
        std::cerr << "FAIL: iRead round-trip mismatch\n";
        std::cerr << "Expected: ";
        for (auto r : expectedReads) std::cerr << r << " ";
        std::cerr << "\nGot: ";
        for (auto r : iReads) std::cerr << r << " ";
        std::cerr << "\n";
        return 1;
    }
    
    // Verify spill occurred (with maxRAM=64, should have created at least one spill file)
    if (spillCount == 0) {
        std::cerr << "FAIL: No spill files created (expected with maxRAM=64)\n";
        system(("rm -rf " + P.outBAMsortTmpDir).c_str());
        return 1;
    }
    
    // Clean up temp directory
    system(("rm -rf " + P.outBAMsortTmpDir).c_str());

    std::cout << "OK\n";
    return 0;
}

