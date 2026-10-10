#ifndef STREAM_ORDER_EXPERIMENT_H
#define STREAM_ORDER_EXPERIMENT_H

#include "DatasetLoader.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

struct StreamOrderExperimentConfig {
    std::string dataset_name;
    std::vector<float> memory_kb_list;
    float pair_witness_ratio = 0.3f;
    uint32_t cell_num = 4;
    uint32_t reward = 12;
    uint32_t penalty = 1;
    float tau_f_ratio = 0.0001f;
    double tau_s = 1024.0;
    uint32_t random_seed = 20260901;
    uint32_t block_shuffle_seed = 20260901;
    std::size_t block_size = 1000000;
};

void runStreamOrderRobustness(const LoadedDataset& loaded,
                              const StreamOrderExperimentConfig& config);

#endif

