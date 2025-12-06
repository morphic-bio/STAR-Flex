#include "ProbeListIndex.h"
#include <fstream>
#include <sstream>

bool ProbeListIndex::load(const std::string &path) {
    geneIdToIndex_.clear();
    if (path.empty() || path == "-") return true; // treat empty as not provided
    std::ifstream in(path.c_str());
    if (!in.is_open()) return false;
    std::string line;
    uint32_t lineNo = 0;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '#') continue;
        // trim
        size_t beg = line.find_first_not_of(" \t\r\n");
        size_t end = line.find_last_not_of(" \t\r\n");
        if (beg == std::string::npos) continue;
        std::string geneId = line.substr(beg, end - beg + 1);
        lineNo++;
        // 1-based index stored; guard against overflow but defer to geneIndex15()
        geneIdToIndex_[geneId] = lineNo;
    }
    return true;
}


