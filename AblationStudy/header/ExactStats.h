#ifndef EXACT_STATS_H
#define EXACT_STATS_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>

struct ExactStats {
    std::unordered_map<uint32_t, uint32_t> true_flow_frequency;
    std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> element_frequency;

    void update(uint32_t flow_id, uint32_t element_id) {
        ++true_flow_frequency[flow_id];
        ++element_frequency[flow_id][element_id];
    }

    uint32_t trueFrequency(uint32_t flow_id) const {
        const auto iter = true_flow_frequency.find(flow_id);
        return iter == true_flow_frequency.end() ? 0 : iter->second;
    }

    uint32_t trueCardinality(uint32_t flow_id) const {
        const auto element_iter = element_frequency.find(flow_id);
        return element_iter == element_frequency.end()
               ? 0
               : static_cast<uint32_t>(element_iter->second.size());
    }

    uint64_t totalItemCount() const {
        uint64_t total_item_num = 0;
        for (const auto& item : true_flow_frequency) {
            total_item_num += static_cast<uint64_t>(item.second);
        }
        return total_item_num;
    }

    uint64_t uniquePairCount() const {
        uint64_t unique_pair_num = 0;
        for (const auto& item : element_frequency) {
            unique_pair_num += static_cast<uint64_t>(item.second.size());
        }
        return unique_pair_num;
    }
};

struct TruthSnapshot {
    float tau_f_ratio = 0.0001f;
    uint32_t tau_f = 1;
    uint32_t reward = 24;
    uint32_t penalty = 1;
    double tau_s = 1024.0;
    uint32_t truth_count = 0;
    std::unordered_set<uint32_t> diffuse_flows;

    bool contains(uint32_t flow_id) const {
        return diffuse_flows.find(flow_id) != diffuse_flows.end();
    }
};

inline uint32_t computeFrequencyThreshold(uint64_t total_item_num, float tau_f_ratio) {
    return std::max<uint32_t>(
            1,
            static_cast<uint32_t>(std::ceil(static_cast<double>(total_item_num) *
                                            static_cast<double>(tau_f_ratio))));
}

inline double computeDiffuseSupport(uint32_t frequency,
                                    uint32_t cardinality,
                                    uint32_t reward,
                                    uint32_t penalty) {
    return static_cast<double>(reward) * static_cast<double>(cardinality) -
           static_cast<double>(penalty) *
           (static_cast<double>(frequency) - static_cast<double>(cardinality));
}

inline bool matchesDiffuseRule(uint32_t frequency,
                               double diffuse_support,
                               uint32_t tau_f,
                               double tau_s) {
    return frequency >= tau_f && diffuse_support >= tau_s;
}

inline TruthSnapshot buildTruthSnapshot(const ExactStats& exact_stats,
                                       float tau_f_ratio,
                                       uint32_t reward,
                                       uint32_t penalty,
                                       double tau_s) {
    TruthSnapshot snapshot;
    snapshot.tau_f_ratio = tau_f_ratio;
    snapshot.tau_f = computeFrequencyThreshold(exact_stats.totalItemCount(), tau_f_ratio);
    snapshot.reward = reward;
    snapshot.penalty = penalty;
    snapshot.tau_s = tau_s;
    snapshot.diffuse_flows.reserve(exact_stats.true_flow_frequency.size());

    for (const auto& item : exact_stats.true_flow_frequency) {
        const uint32_t flow_id = item.first;
        const uint32_t frequency = item.second;
        const uint32_t cardinality = exact_stats.trueCardinality(flow_id);
        const double diffuse_support =
                computeDiffuseSupport(frequency, cardinality, reward, penalty);
        if (matchesDiffuseRule(frequency, diffuse_support, snapshot.tau_f, tau_s)) {
            snapshot.diffuse_flows.insert(flow_id);
        }
    }

    snapshot.truth_count = static_cast<uint32_t>(snapshot.diffuse_flows.size());
    return snapshot;
}

#endif


