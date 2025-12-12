#include "GlobalVariables.h"
Stats g_statsAll;//global mapping statistics
ThreadControl g_threadChunks;
std::atomic<uint64_t> g_bamRecordIndex{0};

