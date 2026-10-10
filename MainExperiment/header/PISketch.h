#ifndef PI_SKETCH_H
#define PI_SKETCH_H

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>
#include "FlowEstimate.h"

struct ExactStats;
struct TruthSnapshot;

class PISketch {
private:
    struct Cell {
        uint32_t flow_id = 0;
        uint32_t freq = 0;
        int64_t weight = 0;
    };

    class BloomFilter {
    private:
        uint64_t bit_num = 0;
        uint32_t hash_num = 3;
        uint32_t seed1 = 7001;
        uint32_t seed2 = 8001;
        std::vector<uint64_t> words;
        uint64_t set_bit_num = 0;

        bool setBit(uint64_t index);

    public:
        BloomFilter(uint64_t bit_num, uint32_t hash_num);
        bool queryAndInsert(uint32_t flow_id, uint32_t element_id);
        double fillRatio() const;
        double estimatedFalsePositiveRate() const;
    };

    float memory_kb = 0.0f;
    float bloom_ratio = 0.5f;
    uint32_t bloom_hash_num = 3;
    uint32_t cell_num = 4;
    uint32_t reward = 12;
    uint32_t penalty = 1;
    uint32_t flow_id_bits = 32;
    uint32_t freq_bits = 32;
    uint32_t weight_bits = 32;
    uint32_t bucket_seed = 9001;
    uint32_t bucket_num = 1;
    BloomFilter bloom_filter;
    std::vector<Cell> cells;

public:
    explicit PISketch(float memory_kb);

    void update(uint32_t flow_id, uint32_t element_id);

    std::vector<FlowEstimate> collectEstimatesForEvaluation(
            const ExactStats& exact_stats);

    void eval(const std::vector<std::pair<uint32_t, uint32_t>>& data,
              const ExactStats& exact_stats,
              const std::vector<float>& tau_f_ratio_list,
              const std::vector<double>& truth_tau_s_list);
};

#endif

