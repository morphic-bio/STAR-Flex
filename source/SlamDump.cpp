#include "SlamDump.h"
#include <fstream>
#include <cstring>

namespace {
constexpr char kMagic[8] = {'S','L','A','M','D','U','M','P'};

void writeString(std::ofstream& out, const std::string& s) {
    uint32_t len = static_cast<uint32_t>(s.size());
    out.write(reinterpret_cast<const char*>(&len), sizeof(len));
    if (len > 0) {
        out.write(s.data(), len);
    }
}

bool readString(std::ifstream& in, std::string* s, std::string* err) {
    uint32_t len = 0;
    if (!in.read(reinterpret_cast<char*>(&len), sizeof(len))) {
        if (err) *err = "Failed to read string length";
        return false;
    }
    s->assign(len, '\0');
    if (len > 0 && !in.read(&(*s)[0], len)) {
        if (err) *err = "Failed to read string payload";
        return false;
    }
    return true;
}
} // namespace

bool writeSlamDump(const std::string& path,
                   const SlamDumpMetadata& meta,
                   const std::vector<const SlamReadBuffer*>& buffers,
                   uint64_t maxReads,
                   std::string* err) {
    std::ofstream out(path.c_str(), std::ios::binary);
    if (!out.good()) {
        if (err) *err = "Failed to open dump file for writing: " + path;
        return false;
    }

    // Header
    out.write(kMagic, sizeof(kMagic));
    uint32_t version = meta.version;
    uint32_t flags = 0; // reserved
    uint32_t nGenes = static_cast<uint32_t>(meta.geneIds.size());
    uint32_t nChrom = static_cast<uint32_t>(meta.chrNames.size());
    uint64_t nReads = 0;
    out.write(reinterpret_cast<const char*>(&version), sizeof(version));
    out.write(reinterpret_cast<const char*>(&flags), sizeof(flags));
    out.write(reinterpret_cast<const char*>(&nGenes), sizeof(nGenes));
    out.write(reinterpret_cast<const char*>(&nChrom), sizeof(nChrom));
    // placeholder for nReads
    std::streampos nReadsPos = out.tellp();
    out.write(reinterpret_cast<const char*>(&nReads), sizeof(nReads));
    out.write(reinterpret_cast<const char*>(&meta.errorRate), sizeof(meta.errorRate));
    out.write(reinterpret_cast<const char*>(&meta.convRate), sizeof(meta.convRate));

    // Gene metadata
    for (size_t i = 0; i < meta.geneIds.size(); ++i) {
        writeString(out, meta.geneIds[i]);
        if (i < meta.geneNames.size()) {
            writeString(out, meta.geneNames[i]);
        } else {
            writeString(out, meta.geneIds[i]);
        }
    }

    // Chrom metadata
    for (size_t i = 0; i < meta.chrNames.size(); ++i) {
        writeString(out, meta.chrNames[i]);
        uint64_t start = (i < meta.chrStart.size()) ? meta.chrStart[i] : 0;
        out.write(reinterpret_cast<const char*>(&start), sizeof(start));
    }

    // Reads
    uint64_t written = 0;
    for (const auto* buffer : buffers) {
        if (buffer == nullptr) continue;
        for (const SlamBufferedRead& r : buffer->reads()) {
            if (maxReads > 0 && written >= maxReads) break;

            writeString(out, r.readName);
            out.write(reinterpret_cast<const char*>(&r.readLength0), sizeof(r.readLength0));
            out.write(reinterpret_cast<const char*>(&r.readLength1), sizeof(r.readLength1));
            out.write(reinterpret_cast<const char*>(&r.isMinus), sizeof(r.isMinus));
            out.write(reinterpret_cast<const char*>(&r.oppositeStrand), sizeof(r.oppositeStrand));
            out.write(reinterpret_cast<const char*>(&r.isIntronic), sizeof(r.isIntronic));
            out.write(reinterpret_cast<const char*>(&r.fileIndex), sizeof(r.fileIndex));
            out.write(reinterpret_cast<const char*>(&r.weight), sizeof(r.weight));

            uint32_t geneCount = static_cast<uint32_t>(r.geneIds.size());
            out.write(reinterpret_cast<const char*>(&geneCount), sizeof(geneCount));
            if (geneCount > 0) {
                out.write(reinterpret_cast<const char*>(r.geneIds.data()), geneCount * sizeof(uint32_t));
            }

            uint32_t posCount = static_cast<uint32_t>(r.positions.size());
            out.write(reinterpret_cast<const char*>(&posCount), sizeof(posCount));
            for (const SlamBufferedPosition& p : r.positions) {
                out.write(reinterpret_cast<const char*>(&p.readPos), sizeof(p.readPos));
                out.write(reinterpret_cast<const char*>(&p.genomicPos), sizeof(p.genomicPos));
                out.write(reinterpret_cast<const char*>(&p.refBase), sizeof(p.refBase));
                out.write(reinterpret_cast<const char*>(&p.readBase), sizeof(p.readBase));
                out.write(reinterpret_cast<const char*>(&p.qual), sizeof(p.qual));
                out.write(reinterpret_cast<const char*>(&p.secondMate), sizeof(p.secondMate));
                out.write(reinterpret_cast<const char*>(&p.overlap), sizeof(p.overlap));
            }

            ++written;
        }
        if (maxReads > 0 && written >= maxReads) break;
    }

    // finalize nReads
    out.seekp(nReadsPos);
    out.write(reinterpret_cast<const char*>(&written), sizeof(written));
    out.close();
    return true;
}

bool readSlamDump(const std::string& path,
                  SlamDumpMetadata* meta,
                  std::vector<SlamBufferedRead>* reads,
                  std::string* err) {
    if (meta == nullptr || reads == nullptr) {
        if (err) *err = "Null output pointers for readSlamDump";
        return false;
    }
    std::ifstream in(path.c_str(), std::ios::binary);
    if (!in.good()) {
        if (err) *err = "Failed to open dump file for reading: " + path;
        return false;
    }

    char magic[8] = {0};
    if (!in.read(magic, sizeof(magic)) || std::memcmp(magic, kMagic, sizeof(kMagic)) != 0) {
        if (err) *err = "Invalid dump magic";
        return false;
    }
    uint32_t version = 0, flags = 0, nGenes = 0, nChrom = 0;
    uint64_t nReads = 0;
    if (!in.read(reinterpret_cast<char*>(&version), sizeof(version)) ||
        !in.read(reinterpret_cast<char*>(&flags), sizeof(flags)) ||
        !in.read(reinterpret_cast<char*>(&nGenes), sizeof(nGenes)) ||
        !in.read(reinterpret_cast<char*>(&nChrom), sizeof(nChrom)) ||
        !in.read(reinterpret_cast<char*>(&nReads), sizeof(nReads)) ||
        !in.read(reinterpret_cast<char*>(&meta->errorRate), sizeof(meta->errorRate)) ||
        !in.read(reinterpret_cast<char*>(&meta->convRate), sizeof(meta->convRate))) {
        if (err) *err = "Failed to read dump header";
        return false;
    }
    meta->version = version;
    meta->nReads = nReads;

    meta->geneIds.clear();
    meta->geneNames.clear();
    meta->chrNames.clear();
    meta->chrStart.clear();
    meta->geneIds.reserve(nGenes);
    meta->geneNames.reserve(nGenes);
    for (uint32_t i = 0; i < nGenes; ++i) {
        std::string gid, gname;
        if (!readString(in, &gid, err) || !readString(in, &gname, err)) {
            return false;
        }
        meta->geneIds.push_back(std::move(gid));
        meta->geneNames.push_back(std::move(gname));
    }

    meta->chrNames.reserve(nChrom);
    meta->chrStart.reserve(nChrom);
    for (uint32_t i = 0; i < nChrom; ++i) {
        std::string cname;
        if (!readString(in, &cname, err)) {
            return false;
        }
        uint64_t start = 0;
        if (!in.read(reinterpret_cast<char*>(&start), sizeof(start))) {
            if (err) *err = "Failed to read chrStart";
            return false;
        }
        meta->chrNames.push_back(std::move(cname));
        meta->chrStart.push_back(start);
    }

    reads->clear();
    reads->reserve(static_cast<size_t>(nReads));
    for (uint64_t i = 0; i < nReads; ++i) {
        SlamBufferedRead r;
        if (!readString(in, &r.readName, err)) return false;
        if (!in.read(reinterpret_cast<char*>(&r.readLength0), sizeof(r.readLength0)) ||
            !in.read(reinterpret_cast<char*>(&r.readLength1), sizeof(r.readLength1)) ||
            !in.read(reinterpret_cast<char*>(&r.isMinus), sizeof(r.isMinus)) ||
            !in.read(reinterpret_cast<char*>(&r.oppositeStrand), sizeof(r.oppositeStrand)) ||
            !in.read(reinterpret_cast<char*>(&r.isIntronic), sizeof(r.isIntronic)) ||
            !in.read(reinterpret_cast<char*>(&r.fileIndex), sizeof(r.fileIndex)) ||
            !in.read(reinterpret_cast<char*>(&r.weight), sizeof(r.weight))) {
            if (err) *err = "Failed to read read header";
            return false;
        }
        uint32_t geneCount = 0;
        if (!in.read(reinterpret_cast<char*>(&geneCount), sizeof(geneCount))) {
            if (err) *err = "Failed to read geneCount";
            return false;
        }
        r.geneIds.resize(geneCount);
        if (geneCount > 0 &&
            !in.read(reinterpret_cast<char*>(r.geneIds.data()), geneCount * sizeof(uint32_t))) {
            if (err) *err = "Failed to read geneIds";
            return false;
        }
        uint32_t posCount = 0;
        if (!in.read(reinterpret_cast<char*>(&posCount), sizeof(posCount))) {
            if (err) *err = "Failed to read posCount";
            return false;
        }
        r.positions.resize(posCount);
        for (uint32_t p = 0; p < posCount; ++p) {
            SlamBufferedPosition bp;
            if (!in.read(reinterpret_cast<char*>(&bp.readPos), sizeof(bp.readPos)) ||
                !in.read(reinterpret_cast<char*>(&bp.genomicPos), sizeof(bp.genomicPos)) ||
                !in.read(reinterpret_cast<char*>(&bp.refBase), sizeof(bp.refBase)) ||
                !in.read(reinterpret_cast<char*>(&bp.readBase), sizeof(bp.readBase)) ||
                !in.read(reinterpret_cast<char*>(&bp.qual), sizeof(bp.qual)) ||
                !in.read(reinterpret_cast<char*>(&bp.secondMate), sizeof(bp.secondMate)) ||
                !in.read(reinterpret_cast<char*>(&bp.overlap), sizeof(bp.overlap))) {
                if (err) *err = "Failed to read position record";
                return false;
            }
            r.positions[p] = bp;
        }
        reads->push_back(std::move(r));
    }
    return true;
}
