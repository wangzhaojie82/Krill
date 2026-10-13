#ifndef COCO_SKETCH_BASELINE_H
#define COCO_SKETCH_BASELINE_H

#include <cstdint>
#include <random>
#include <utility>
#include <vector>

#include "FlowEstimate.h"

struct ExactStats;

class CocoSketch {
private:
    struct Bucket {
        uint64_t key = 0;
        uint32_t value = 0;
    };

    float memory_kb = 0.0f;
    uint32_t depth = 2;
    uint32_t width = 1;
    uint32_t reward = 12;
    uint32_t penalty = 1;
    std::vector<std::vector<Bucket>> table;
    std::vector<uint32_t> seeds;
    std::mt19937 rng;
    std::uniform_real_distribution<double> probability;

    static uint64_t combineKey(uint32_t flow_id, uint32_t element_id);
    static std::pair<uint32_t, uint32_t> splitKey(uint64_t key);
    uint32_t bucketIndex(uint64_t key, uint32_t row) const;

public:
    explicit CocoSketch(float memory_kb);

    void update(uint32_t flow_id, uint32_t element_id);

    std::vector<FlowEstimate> collectEstimatesForEvaluation(
            const ExactStats& exact_stats) const;

    void eval(const std::vector<std::pair<uint32_t, uint32_t>>& data,
              const ExactStats& exact_stats,
              const std::vector<float>& tau_f_ratio_list,
              const std::vector<double>& tau_s_list);
};

#endif
