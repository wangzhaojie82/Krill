#ifndef ABLATION_SKETCH_H
#define ABLATION_SKETCH_H

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

struct ExactStats;

enum class PairMechanism {
    PairWitness,
    BloomFilter,
    Exact
};

enum class ResidenceMechanism {
    Ideal,
    FlowResidence
};

struct AblationConfig {
    PairMechanism pair_mechanism = PairMechanism::PairWitness;
    ResidenceMechanism residence_mechanism = ResidenceMechanism::FlowResidence;
    double total_memory_kb = 400.0;
    double pair_memory_ratio = 0.3;
    double residence_memory_ratio = 0.7;
    uint32_t cell_num = 4;
    uint32_t reward = 12;
    uint32_t penalty = 1;
    uint32_t bloom_hash_num = 3;
    double decay_base = 1.08;
};

struct AblationResult {
    std::string variant;
    double total_memory_kb = 0.0;
    double pair_memory_kb = 0.0;
    bool pair_memory_limited = true;
    double residence_memory_kb = 0.0;
    bool residence_memory_limited = true;
    uint32_t cell_num = 0;
    uint32_t tau_f = 0;
    double tau_s = 0.0;
    uint64_t resident_flows = 0;
    uint64_t ground_truth_flows = 0;
    uint64_t predicted_flows = 0;
    uint64_t tp = 0;
    uint64_t fp = 0;
    uint64_t fn = 0;
    double precision = 0.0;
    double recall = 0.0;
    double f1 = 0.0;
    double dss_are = 0.0;
    double frequency_are = 0.0;
    double update_mips = 0.0;
};

class AblationSketch {
public:
    explicit AblationSketch(const AblationConfig& config);

    std::vector<AblationResult> run(
            const std::vector<std::pair<uint32_t, uint32_t>>& data,
            const ExactStats& exact_stats,
            const std::vector<float>& tau_f_ratio_list,
            const std::vector<double>& tau_s_list) const;

    static const char* pairMechanismName(PairMechanism mechanism);
    static const char* residenceMechanismName(ResidenceMechanism mechanism);

private:
    AblationConfig config_;
};

#endif

