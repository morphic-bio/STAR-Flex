#ifndef GLOBAL_VARIABLES_DEF
#define GLOBAL_VARIABLES_DEF

#include "ThreadControl.h"
#include <atomic>

extern Stats g_statsAll;
extern ThreadControl g_threadChunks;
extern std::atomic<uint64_t> g_bamRecordIndex;

#endif

