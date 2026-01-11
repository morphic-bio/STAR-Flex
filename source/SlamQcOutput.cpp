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
                     uint64_t readsAnalyzed,
                     const SlamVarianceTrimResult* trimResult,
                     const std::string& trimSource) {
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
    out << "  \"version\": \"1.2\",\n";
    out << "  \"algorithm\": \"segmented_regression\",\n";
    out << "  \"file_index\": " << fileIndex << ",\n";
    out << "  \"trim_scope\": \"" << trimScope << "\",\n";
    if (!trimSource.empty()) {
        out << "  \"trim_source\": \"" << trimSource << "\",\n";
    }
    out << "  \"trim5p\": " << trim5p << ",\n";
    out << "  \"trim3p\": " << trim3p << ",\n";
    out << "  \"reads_analyzed\": " << readsAnalyzed << ",\n";
    
    // Include segmented regression info if available
    if (trimResult != nullptr) {
        out << "  \"segmented_regression\": {\n";
        out << "    \"breakpoint_b1\": " << trimResult->kneeBin5p << ",\n";
        out << "    \"breakpoint_b2\": " << trimResult->kneeBin3p << ",\n";
        out << "    \"total_sse\": " << std::fixed << std::setprecision(6) << trimResult->totalSSE << ",\n";
        out << "    \"mode\": \"" << trimResult->mode << "\",\n";
        out << "    \"segment1\": {\n";
        out << "      \"slope\": " << std::fixed << std::setprecision(6) << trimResult->seg1.slope << ",\n";
        out << "      \"intercept\": " << std::fixed << std::setprecision(6) << trimResult->seg1.intercept << ",\n";
        out << "      \"sse\": " << std::fixed << std::setprecision(6) << trimResult->seg1.sse << "\n";
        out << "    },\n";
        out << "    \"segment2\": {\n";
        out << "      \"slope\": " << std::fixed << std::setprecision(6) << trimResult->seg2.slope << ",\n";
        out << "      \"intercept\": " << std::fixed << std::setprecision(6) << trimResult->seg2.intercept << ",\n";
        out << "      \"sse\": " << std::fixed << std::setprecision(6) << trimResult->seg2.sse << "\n";
        out << "    },\n";
        out << "    \"segment3\": {\n";
        out << "      \"slope\": " << std::fixed << std::setprecision(6) << trimResult->seg3.slope << ",\n";
        out << "      \"intercept\": " << std::fixed << std::setprecision(6) << trimResult->seg3.intercept << ",\n";
        out << "      \"sse\": " << std::fixed << std::setprecision(6) << trimResult->seg3.sse << "\n";
        out << "    },\n";
        
        // Include smoothed curve if available
        if (!trimResult->smoothedCurve.empty()) {
            out << "    \"smoothed_stdev_curve\": [";
            for (size_t i = 0; i < trimResult->smoothedCurve.size(); ++i) {
                if (i > 0) out << ", ";
                if (std::isnan(trimResult->smoothedCurve[i])) {
                    out << "null";
                } else {
                    out << std::fixed << std::setprecision(6) << trimResult->smoothedCurve[i];
                }
            }
            out << "]\n";
        } else {
            out << "    \"smoothed_stdev_curve\": []\n";
        }
        out << "  },\n";
    }
    
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
        double stddevTcRate = s.stddevTcRate();
        
        out << "    {\n";
        out << "      \"position\": " << pos << ",\n";
        out << "      \"read_count\": " << s.readCount << ",\n";
        out << "      \"t_count\": " << s.tCount << ",\n";
        out << "      \"tc_count\": " << s.tcCount << ",\n";
        out << "      \"mean_qual\": " << std::fixed << std::setprecision(2) << meanQual << ",\n";
        out << "      \"stddev_qual\": " << std::fixed << std::setprecision(2) << stddevQual << ",\n";
        out << "      \"variance_qual\": " << std::fixed << std::setprecision(4) << s.varianceQual() << ",\n";
        out << "      \"mean_tc_rate\": " << std::fixed << std::setprecision(6) << meanTcRate << ",\n";
        out << "      \"stddev_tc_rate\": " << std::fixed << std::setprecision(6) << stddevTcRate << ",\n";
        out << "      \"variance_tc_rate\": " << std::fixed << std::setprecision(6) << s.varianceTcRate() << "\n";
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
    out << "    body { font-family: 'Segoe UI', Arial, sans-serif; margin: 20px; background: #f5f5f5; }\n";
    out << "    .container { max-width: 1200px; margin: 0 auto; }\n";
    out << "    h1 { color: #333; }\n";
    out << "    .summary { background: white; padding: 15px; border-radius: 8px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }\n";
    out << "    .summary-item { display: inline-block; margin-right: 30px; }\n";
    out << "    .summary-label { color: #666; font-size: 12px; }\n";
    out << "    .summary-value { font-size: 24px; font-weight: bold; color: #2196F3; }\n";
    out << "    .plot-container { background: white; padding: 15px; border-radius: 8px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }\n";
    out << "  </style>\n";
    out << "</head>\n";
    out << "<body>\n";
    out << "  <div class=\"container\">\n";
    out << "    <h1>SLAM Auto-Trim QC Report</h1>\n";
    out << "    <div class=\"summary\" id=\"summary\"></div>\n";
    out << "    <div class=\"plot-container\" id=\"tcStddevPlot\"></div>\n";
    out << "    <div class=\"plot-container\" id=\"qualityPlot\"></div>\n";
    out << "  </div>\n";
    out << "  <script>\n";
    out << "    fetch('" << jsonFileName << "')\n";
    out << "      .then(response => response.json())\n";
    out << "      .then(data => {\n";
    out << "        // Summary\n";
    out << "        const summaryHtml = `\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">Algorithm</div><div class=\"summary-value\">${data.algorithm || 'segmented_regression'}</div></div>\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">Trim 5'</div><div class=\"summary-value\">${data.trim5p}</div></div>\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">Trim 3'</div><div class=\"summary-value\">${data.trim3p}</div></div>\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">Reads Analyzed</div><div class=\"summary-value\">${data.reads_analyzed.toLocaleString()}</div></div>\n";
    out << "          <div class=\"summary-item\"><div class=\"summary-label\">File Index</div><div class=\"summary-value\">${data.file_index}</div></div>\n";
    out << "        `;\n";
    out << "        document.getElementById('summary').innerHTML = summaryHtml;\n";
    out << "        \n";
    out << "        const positions = data.positions.map(p => p.position);\n";
    out << "        const stddevTcRate = data.positions.map(p => p.stddev_tc_rate);\n";
    out << "        const meanQual = data.positions.map(p => p.mean_qual);\n";
    out << "        const stddevQual = data.positions.map(p => p.stddev_qual);\n";
    out << "        \n";
    out << "        // T→C Stdev plot with segmented regression\n";
    out << "        const traces1 = [{\n";
    out << "          x: positions,\n";
    out << "          y: stddevTcRate,\n";
    out << "          type: 'scatter',\n";
    out << "          mode: 'lines+markers',\n";
    out << "          name: 'T→C Stdev (raw)',\n";
    out << "          marker: { size: 4, color: '#2196F3' },\n";
    out << "          line: { width: 1 }\n";
    out << "        }];\n";
    out << "        \n";
    out << "        // Add smoothed curve and segment fits if available\n";
    out << "        if (data.segmented_regression && data.segmented_regression.smoothed_stdev_curve) {\n";
    out << "          const smoothed = data.segmented_regression.smoothed_stdev_curve;\n";
    out << "          traces1.push({\n";
    out << "            x: positions.slice(0, smoothed.length),\n";
    out << "            y: smoothed,\n";
    out << "            type: 'scatter',\n";
    out << "            mode: 'lines',\n";
    out << "            name: 'Smoothed',\n";
    out << "            line: { width: 2, color: '#333', dash: 'dash' }\n";
    out << "          });\n";
    out << "          \n";
    out << "          // Segment fits\n";
    out << "          const b1 = data.segmented_regression.breakpoint_b1;\n";
    out << "          const b2 = data.segmented_regression.breakpoint_b2;\n";
    out << "          const seg1 = data.segmented_regression.segment1;\n";
    out << "          const seg2 = data.segmented_regression.segment2;\n";
    out << "          const seg3 = data.segmented_regression.segment3;\n";
    out << "          \n";
    out << "          if (b1 > 0) {\n";
    out << "            const x1 = Array.from({length: b1}, (_, i) => i);\n";
    out << "            const y1 = x1.map(x => seg1.slope * x + seg1.intercept);\n";
    out << "            traces1.push({ x: x1, y: y1, type: 'scatter', mode: 'lines', name: 'Seg1 (5\\' artifact)', line: { width: 2, color: '#f44336' } });\n";
    out << "          }\n";
    out << "          \n";
    out << "          const x2 = Array.from({length: b2 - b1 + 1}, (_, i) => b1 + i);\n";
    out << "          const y2 = x2.map(x => seg2.slope * x + seg2.intercept);\n";
    out << "          traces1.push({ x: x2, y: y2, type: 'scatter', mode: 'lines', name: 'Seg2 (signal)', line: { width: 2, color: '#4CAF50' } });\n";
    out << "          \n";
    out << "          if (b2 < smoothed.length - 1) {\n";
    out << "            const x3 = Array.from({length: smoothed.length - b2 - 1}, (_, i) => b2 + 1 + i);\n";
    out << "            const y3 = x3.map(x => seg3.slope * x + seg3.intercept);\n";
    out << "            traces1.push({ x: x3, y: y3, type: 'scatter', mode: 'lines', name: 'Seg3 (3\\' artifact)', line: { width: 2, color: '#ff9800' } });\n";
    out << "          }\n";
    out << "        }\n";
    out << "        \n";
    out << "        const shapes = [];\n";
    out << "        if (data.trim5p > 0) {\n";
    out << "          shapes.push({ type: 'line', x0: data.trim5p, x1: data.trim5p, y0: 0, y1: 1, yref: 'paper', line: { color: 'purple', width: 2, dash: 'dot' } });\n";
    out << "        }\n";
    out << "        if (data.trim3p > 0) {\n";
    out << "          const trim3pPos = positions[positions.length - 1] - data.trim3p;\n";
    out << "          shapes.push({ type: 'line', x0: trim3pPos, x1: trim3pPos, y0: 0, y1: 1, yref: 'paper', line: { color: 'purple', width: 2, dash: 'dot' } });\n";
    out << "        }\n";
    out << "        \n";
    out << "        Plotly.newPlot('tcStddevPlot', traces1, {\n";
    out << "          title: 'T→C Rate Standard Deviation by Position (Segmented Regression)',\n";
    out << "          xaxis: { title: 'Read Position' },\n";
    out << "          yaxis: { title: 'T→C Stdev' },\n";
    out << "          shapes: shapes,\n";
    out << "          hovermode: 'closest',\n";
    out << "          annotations: [\n";
    out << "            { x: data.trim5p, y: 1, yref: 'paper', text: 'trim5p=' + data.trim5p, showarrow: false, yanchor: 'bottom' },\n";
    out << "            { x: positions[positions.length - 1] - data.trim3p, y: 1, yref: 'paper', text: 'trim3p=' + data.trim3p, showarrow: false, yanchor: 'bottom' }\n";
    out << "          ]\n";
    out << "        });\n";
    out << "        \n";
    out << "        // Quality plot\n";
    out << "        Plotly.newPlot('qualityPlot', [\n";
    out << "          { x: positions, y: meanQual, type: 'scatter', mode: 'lines+markers', name: 'Mean Quality', marker: { size: 4 } },\n";
    out << "          { x: positions, y: stddevQual, type: 'scatter', mode: 'lines+markers', name: 'Quality Stdev', yaxis: 'y2', marker: { size: 4 } }\n";
    out << "        ], {\n";
    out << "          title: 'Quality Statistics by Position',\n";
    out << "          xaxis: { title: 'Read Position' },\n";
    out << "          yaxis: { title: 'Mean Quality', side: 'left' },\n";
    out << "          yaxis2: { title: 'Quality Stdev', side: 'right', overlaying: 'y' },\n";
    out << "          shapes: shapes,\n";
    out << "          hovermode: 'closest'\n";
    out << "        });\n";
    out << "      });\n";
    out << "  </script>\n";
    out << "</body>\n";
    out << "</html>\n";
    
    return out.good();
}
