#ifndef COMPENSATION_ABLATION_H
#define COMPENSATION_ABLATION_H

#include "AblationSketch.h"

#include <cstdint>
#include <vector>
#include <utility>

struct ExactStats;

enum class CompensationMechanism {
    Compensated,
    Uncompensated
};

struct CompensationAblationConfig {
    double total_memory_kb = 400.0;
    double pair_memory_ratio = 0.3;
    double residence_memory_ratio = 0.7;
    uint32_t cell_num = 4;
    uint32_t reward = 12;
    uint32_t penalty = 1;
    double decay_base = 1.08;
};

class CompensationAblation {
public:
    explicit CompensationAblation(const CompensationAblationConfig& config);

    std::vector<AblationResult> run(
            const std::vector<std::pair<uint32_t, uint32_t>>& data,
            const ExactStats& exact_stats,
            const std::vector<float>& tau_f_ratio_list,
            const std::vector<double>& tau_s_list,
            CompensationMechanism compensation,
            ResidenceMechanism residence) const;

    static const char* compensationName(CompensationMechanism mechanism);
    static const char* residenceName(ResidenceMechanism mechanism);

private:
    CompensationAblationConfig config_;
};

#endif

