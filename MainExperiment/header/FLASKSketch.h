#ifndef FLASK_SKETCH_H
#define FLASK_SKETCH_H

#include <cstdint>
#include <utility>
#include <vector>
#include "FlowEstimate.h"

struct ExactStats;

class FLASKSketch {
private:
    struct Bucket {
        uint32_t key = 0;
        uint32_t frequency = 0;
        int32_t auxiliary = 0;
        uint8_t leading_zero = 0;
        bool occupied = false;
        double hll_sum = 0.0;
        uint32_t hll_zero_num = 0;
    };

    float memory_kb = 0.0f;
    uint32_t row_num = 4;
    uint32_t bucket_num = 1;
    uint32_t hll_register_num = 16;
    double alpha = 1.75;
    uint32_t reward = 12;
    uint32_t penalty = 1;
    uint32_t bucket_seed = 3001;
    uint32_t pair_seed = 4001;
    uint32_t hll_seed = 5001;
    uint32_t hll_index_bits = 5;
    double hll_alpha = 0.0;
    std::vector<Bucket> buckets;
    std::vector<uint8_t> hll_registers;

    size_t bucketOffset(uint32_t row, uint32_t bucket) const;
    size_t registerOffset(size_t bucket_offset, uint32_t register_index) const;
    uint32_t hash32(uint32_t value, uint32_t seed) const;
    uint32_t hashPair(uint32_t first, uint32_t second, uint32_t seed) const;
    uint8_t leadingZeros(uint32_t hash_value) const;
    uint32_t bucketIndex(uint32_t flow_id, uint32_t row) const;
    void resetHLL(size_t bucket_offset);
    void updateHLL(size_t bucket_offset, uint32_t register_index, uint8_t rank);
    uint32_t queryHLL(size_t bucket_offset) const;
    void reinitializeBucket(size_t bucket_offset,
                            uint32_t flow_id,
                            uint32_t hll_register_index,
                            uint8_t hll_rank,
                            uint8_t leading_zero);

public:
    explicit FLASKSketch(float memory_kb);

    void update(uint32_t flow_id, uint32_t element_id);
    std::pair<uint32_t, uint32_t> query(uint32_t flow_id) const;

    std::vector<FlowEstimate> collectEstimatesForEvaluation(
            const ExactStats& exact_stats) const;

    void eval(const std::vector<std::pair<uint32_t, uint32_t>>& data,
              const ExactStats& exact_stats,
              const std::vector<float>& tau_f_ratio_list,
              const std::vector<double>& tau_s_list);
};

#endif

