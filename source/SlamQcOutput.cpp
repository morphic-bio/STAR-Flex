#include "SlamQcOutput.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>

bool writeSlamQcJson(const SlamVarianceAnalyzer& analyzer,
                     const std::string& outputPath,
                     uint32_t fileIndex,
                     const std::string& trimScope,
                     int trim5p,
                     int trim3p,
                     uint64_t readsAnalyzed) {
    std::ofstream out(outputPath.c_str());
    if (!out.good()) {
        return false;
    }
    
    const auto& stats = analyzer.getStats();
    
    // Find max position
    uint32_t maxPos = 0;
    for (const auto& kv : stats) {
        if (kv.first > maxPos) {
            maxPos = kv.first;
        }
    }
    
    out << "{\n";
    out << "  \"version\": \"1.0\",\n";
    out << "  \"file_index\": " << fileIndex << ",\n";
    out << "  \"trim_scope\": \"" << trimScope << "\",\n";
    out << "  \"trim5p\": " << trim5p << ",\n";
    out << "  \"trim3p\": " << trim3p << ",\n";
    out << "  \"reads_analyzed\": " << readsAnalyzed << ",\n";
    out << "  \"positions\": [\n";
    
    bool first = true;
    for (uint32_t pos = 0; pos <= maxPos; ++pos) {
        auto it = stats.find(pos);
        if (it == stats.end()) {
            continue;
        }
        
        if (!first) {
            out << ",\n";
        }
        first = false;
        
        const auto& s = it->second;
        double meanQual = s.meanQual();
        double stddevQual = s.stddevQual();
        double meanTcRate = s.meanTcRate();
        
        out << "    {\n";
        out << "      \"position\": " << pos << ",\n";
        out << "      \"read_count\": " << s.readCount << ",\n";
        out << "      \"t_count\": " << s.tCount << ",\n";
        out << "      \"tc_count\": " << s.tcCount << ",\n";
        out << "      \"mean_qual\": " << std::fixed << std::setprecision(2) << meanQual << ",\n";
        out << "      \"stddev_qual\": " << std::fixed << std::setprecision(2) << stddevQual << ",\n";
        out << "      \"variance_qual\": " << std::fixed << std::setprecision(4) << s.varianceQual() << ",\n";
        out << "      \"mean_tc_rate\": " << std::fixed << std::setprecision(4) << meanTcRate << "\n";
        out << "    }";
    }
    
    out << "\n  ]\n";
    out << "}\n";
    
    return out.good();
}

bool writeSlamQcHtml(const std::string& jsonPath,
                     const std::string& htmlPath,
                     uint32_t fileIndex) {
    std::ofstream out(htmlPath.c_str());
    if (!out.good()) {
        return false;
    }
    
    std::string jsonFileName = jsonPath;
    size_t lastSlash = jsonPath.find_last_of("/\\");
    if (lastSlash != std::string::npos) {
        jsonFileName = jsonPath.substr(lastSlash + 1);
    }
    
    out << "<!DOCTYPE html>\n";
    out << "<html>\n";
    out << "<head>\n";
    out << "  <title>SLAM QC Report - File " << fileIndex << "</title>\n";
    out << "  <script src=\"https://cdn.plot.ly/plotly-latest.min.js\"></script>\n";
    out << "  <style>\n";
    out << "    body { font-family: Arial, sans-serif; margin: 20px; }\n";
    out << "    .plot-container { margin: 20px 0; }\n";
    out << "  </style>\n";
    out << "</head>\n";
    out << "<body>\n";
    out << "  <h1>SLAM Quality Control Report</h1>\n";
    out << "  <p>File Index: " << fileIndex << "</p>\n";
    out << "  <div id=\"plots\"></div>\n";
    out << "  <script>\n";
    out << "    fetch('" << jsonFileName << "')\n";
    out << "      .then(response => response.json())\n";
    out << "      .then(data => {\n";
    out << "        const positions = data.positions.map(p => p.position);\n";
    out << "        const meanQual = data.positions.map(p => p.mean_qual);\n";
    out << "        const stddevQual = data.positions.map(p => p.stddev_qual);\n";
    out << "        const varianceQual = data.positions.map(p => p.variance_qual);\n";
    out << "        const meanTcRate = data.positions.map(p => p.mean_tc_rate);\n";
    out << "        \n";
    out << "        const trace1 = {\n";
    out << "          x: positions,\n";
    out << "          y: meanQual,\n";
    out << "          type: 'scatter',\n";
    out << "          mode: 'lines+markers',\n";
    out << "          name: 'Mean Quality',\n";
    out << "          yaxis: 'y'\n";
    out << "        };\n";
    out << "        \n";
    out << "        const trace2 = {\n";
    out << "          x: positions,\n";
    out << "          y: varianceQual,\n";
    out << "          type: 'scatter',\n";
    out << "          mode: 'lines+markers',\n";
    out << "          name: 'Quality Variance',\n";
    out << "          yaxis: 'y2'\n";
    out << "        };\n";
    out << "        \n";
    out << "        const layout = {\n";
    out << "          title: 'Quality Statistics by Position',\n";
    out << "          xaxis: { title: 'Read Position' },\n";
    out << "          yaxis: { title: 'Mean Quality', side: 'left' },\n";
    out << "          yaxis2: { title: 'Variance', side: 'right', overlaying: 'y' },\n";
    out << "          hovermode: 'closest'\n";
    out << "        };\n";
    out << "        \n";
    out << "        Plotly.newPlot('plots', [trace1, trace2], layout);\n";
    out << "      });\n";
    out << "  </script>\n";
    out << "</body>\n";
    out << "</html>\n";
    
    return out.good();
}
