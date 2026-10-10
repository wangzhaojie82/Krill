#ifndef KRILL_H
#define KRILL_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>
#include "FlowEstimate.h"

struct ExactStats;
struct TruthSnapshot;

struct KrillEvalResult {
    float memory_kb = 0.0f;
    float pair_witness_ratio = 0.0f;
    uint32_t cell_num = 0;
    uint32_t reward = 0;
    uint32_t penalty = 1;
    double register_num = 0.0;
    double bucket_num = 0.0;
    double final_q = 0.0;
    double state_change_events = 0.0;
    double avg_compensation = 0.0;
    float tau_f_ratio = 0.0f;
    double tau_f = 0.0;
    double diffuse_threshold = 0.0;
    double resident_flows = 0.0;
    double ground_truth_count = 0.0;
    double true_heavy_flows = 0.0;
    double predicted_count = 0.0;
    double resident_frequency_pass = 0.0;
    double resident_dss_pass = 0.0;
    double tp = 0.0;
    double fp = 0.0;
    double fn = 0.0;
    float precision = 0.0f;
    float recall = 0.0f;
    float f1 = 0.0f;
    double frequency_are = 0.0;
    double dss_are = 0.0;
    double fn_not_resident = 0.0;
    double fn_resident = 0.0;
    double fn_frequency_loss = 0.0;
    double fn_dss_loss = 0.0;
    double update_mips = 0.0;
};

class Krill {
private:
    struct Cell {
        double dss = 0.0;
        uint32_t flow_id = 0;
        uint32_t freq = 0;
    };

    float memory_kb;
    float pair_witness_ratio;
    uint32_t cell_num;
    uint32_t reward;
    uint32_t penalty;
    uint32_t flow_id_bits = 32;
    uint32_t dss_bits = 32;
    uint32_t freq_bits = 32;
    uint32_t register_num;
    uint32_t bucket_num;
    std::vector<uint8_t> register_arr;
    std::array<double, 32> event_prob{};
    double event_probability = 1.0;
    uint64_t state_change_events = 0;
    double compensation_sum = 0.0;
    std::vector<Cell> cells;

    static uint8_t rankFromHash(uint32_t hash_value);

    static double exponentialProbabilisticDecay(uint32_t flow_id,
                                                uint32_t element_id,
                                                const Cell& min_cell,
                                                double event_weight);

public:
    explicit Krill(float memory_kb);

    Krill(float memory_kb,
                  float pair_witness_ratio,
                  uint32_t cell_num,
                  uint32_t reward,
                  uint32_t penalty);

    void update(uint32_t flow_id, uint32_t element_id);

    std::vector<FlowEstimate> collectEstimatesForEvaluation(
            const ExactStats& exact_stats);

    void eval(const std::vector<std::pair<uint32_t, uint32_t>>& data,
              const ExactStats& exact_stats,
              const std::vector<float>& tau_f_ratio_list,
              const std::vector<double>& tau_s_list);

    KrillEvalResult evaluateForOrder(
            const std::vector<std::pair<uint32_t, uint32_t>>& data,
            const ExactStats& exact_stats,
            float tau_f_ratio,
            double tau_s);

    void evalForSweep(
            const std::vector<std::pair<uint32_t, uint32_t>>& data,
            const ExactStats& exact_stats,
            const std::vector<float>& tau_f_ratio_list,
            const std::vector<double>& diffuse_threshold_list);
};

#endif



