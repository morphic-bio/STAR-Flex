#ifndef CODE_ProbeListIndex
#define CODE_ProbeListIndex

#include <string>
#include <unordered_map>

class ProbeListIndex {
public:
    // Load mapping from a probe_list.txt file (one gene ID per line; 1-based index)
    // Lines starting with '#' are ignored. Empty lines are skipped.
    bool load(const std::string &path);

    // Returns 0 if not found or if index exceeds 15-bit cap (32767)
    inline uint16_t geneIndex15(const std::string &geneId) const {
        auto it = geneIdToIndex_.find(geneId);
        if (it == geneIdToIndex_.end()) return 0;
        uint32_t idx = it->second;
        return idx <= 0x7FFFu ? static_cast<uint16_t>(idx) : 0;
    }

    inline bool empty() const { return geneIdToIndex_.empty(); }

private:
    std::unordered_map<std::string, uint32_t> geneIdToIndex_;
};

#endif


