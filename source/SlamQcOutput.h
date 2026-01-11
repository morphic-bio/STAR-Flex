#ifndef SLAM_QC_OUTPUT_H
#define SLAM_QC_OUTPUT_H

#include "SlamVarianceAnalysis.h"
#include <string>
#include <cstdint>

// Write QC JSON output with per-position variance stats
bool writeSlamQcJson(const SlamVarianceAnalyzer& analyzer, 
                     const std::string& outputPath,
                     uint32_t fileIndex,
                     const std::string& trimScope,
                     int trim5p,
                     int trim3p,
                     uint64_t readsAnalyzed);

// Write QC HTML report using Plotly CDN
bool writeSlamQcHtml(const std::string& jsonPath,
                     const std::string& htmlPath,
                     uint32_t fileIndex);

#endif // SLAM_QC_OUTPUT_H
