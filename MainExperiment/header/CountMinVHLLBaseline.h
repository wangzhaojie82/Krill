#ifndef COUNT_MIN_VHLL_BASELINE_H
#define COUNT_MIN_VHLL_BASELINE_H

#include "CountMin.h"
#include "vHLL.h"
#include "FlowEstimate.h"
#include <cstdint>
#include <utility>
#include <unordered_map>
#include <vector>

struct ExactStats;

class CountMinVHLLBaseline {
private:
    struct Candidate {
        uint32_t flow_id = 0;
        uint32_t est_f = 0;
    };

    float memory_kb = 0.0f;
    float count_min_ratio = 0.35f;
    float heap_ratio = 0.3f;
    float vhll_ratio = 0.35f;
    uint32_t reward = 12;
    uint32_t penalty = 1;
    uint32_t heap_entry_bits = 64; // flow id + frequency

    CountMin count_min;
    vHLL cardinality_estimator;
    size_t heap_capacity = 1;
    std::vector<Candidate> heap;
    std::unordered_map<uint32_t, size_t> heap_index;

    bool candidateLess(size_t lhs, size_t rhs) const;
    void swapCandidates(size_t lhs, size_t rhs);
    void siftUp(size_t index);
    void siftDown(size_t index);
    void updateCandidate(uint32_t flow_id, uint32_t est_f);

public:
    explicit CountMinVHLLBaseline(float memory_kb);

    void update(uint32_t flow_id, uint32_t element_id);

    std::vector<FlowEstimate> collectEstimatesForEvaluation(
            const ExactStats& exact_stats);

    void eval(const std::vector<std::pair<uint32_t, uint32_t>>& data,
              const ExactStats& exact_stats,
              const std::vector<float>& tau_f_ratio_list,
              const std::vector<double>& tau_s_list);
};

#endif

