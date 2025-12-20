#ifndef CODE_TranscriptQuantOutput
#define CODE_TranscriptQuantOutput

#include "em_types.h"
#include <string>

// Write quantification results to Salmon-compatible quant.sf format
void writeQuantSF(const EMResult& result, const TranscriptState& state, const std::string& filename);

#endif // CODE_TranscriptQuantOutput

