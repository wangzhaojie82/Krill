#ifndef DATASET_LOADER_H
#define DATASET_LOADER_H

#include <cstdint>
#include <utility>
#include <vector>
#include "ExactStats.h"

using PairData = std::vector<std::pair<uint32_t, uint32_t>>;

struct LoadedDataset {
    PairData data;
    ExactStats exact_stats;
};

LoadedDataset loadDatasetECommerceEvents();
LoadedDataset loadDatasetTenrec();
LoadedDataset loadDataSetStackOverflow();

#endif

