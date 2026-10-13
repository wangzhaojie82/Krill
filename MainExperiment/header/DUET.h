#ifndef DUET_BASELINE_H
#define DUET_BASELINE_H

#include <cstdint>
#include <utility>
#include <vector>

#include "CUSketch.h"
#include "FlowEstimate.h"

struct ExactStats;

class DUET {
private:
    struct Bucket {
        uint64_t element = 0;
        uint32_t count = 0;
    };

    float memory_kb = 0.0f;
    float cu_ratio = 0.25f;
    float filter_ratio = 0.35f;
    float stable_ratio = 0.40f;
    uint32_t heavy_flow_threshold = 1000;
    uint32_t reward = 12;
    uint32_t penalty = 1;

    CUSketch cu_sketch;
    uint32_t filter_depth = 4;
    uint32_t filter_width = 1;
    std::vector<std::vector<Bucket>> filter;
    uint32_t stable_rows = 50;
    uint32_t stable_width = 1;
    std::vector<std::vector<Bucket>> stable;
    std::vector<uint32_t> filter_seeds;

    static uint64_t combineKey(uint32_t flow_id, uint32_t element_id);
    static std::pair<uint32_t, uint32_t> splitKey(uint64_t key);
    void insertIntoFilter(uint32_t flow_id, uint32_t element_id);
    void insertIntoStable(uint32_t flow_id, uint32_t element_id, uint32_t count);

public:
    explicit DUET(float memory_kb);

    void update(uint32_t flow_id, uint32_t element_id);

    std::vector<FlowEstimate> collectEstimatesForEvaluation(
            const ExactStats& exact_stats) const;

    void eval(const std::vector<std::pair<uint32_t, uint32_t>>& data,
              const ExactStats& exact_stats,
              const std::vector<float>& tau_f_ratio_list,
              const std::vector<double>& tau_s_list);
};

#endif
