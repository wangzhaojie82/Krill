#ifndef REPLACEMENT_ABLATION_H
#define REPLACEMENT_ABLATION_H

#include "AblationSketch.h"

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

struct ExactStats;

enum class ReplacementMechanism {
    ProbabilisticDecay,
    DirectDecay,
    DirectReplace
};

struct ReplacementAblationConfig {
    double total_memory_kb = 400.0;
    double pair_memory_ratio = 0.3;
    double residence_memory_ratio = 0.7;
    uint32_t cell_num = 4;
    uint32_t reward = 12;
    uint32_t penalty = 1;
    double decay_base = 1.08;
};

class ReplacementAblation {
public:
    explicit ReplacementAblation(const ReplacementAblationConfig& config);

    std::vector<AblationResult> run(
            const std::vector<std::pair<uint32_t, uint32_t>>& data,
            const ExactStats& exact_stats,
            const std::vector<float>& tau_f_ratio_list,
            const std::vector<double>& tau_s_list,
            ReplacementMechanism mechanism) const;

    static const char* mechanismName(ReplacementMechanism mechanism);

private:
    ReplacementAblationConfig config_;
};

#endif

