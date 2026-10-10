#include "header/ReplacementAblation.h"

#include "header/ExactStats.h"
#include "header/MurmurHash3.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

constexpr uint32_t kRegisterSeed = 2001;
constexpr uint32_t kRankSeed = 3001;
constexpr uint32_t kBucketSeed = 1001;
constexpr uint32_t kDecaySeed = 6001;
constexpr uint32_t kFlowIdBits = 32;
constexpr uint32_t kFrequencyBits = 32;
constexpr uint32_t kDssBits = 32;

uint64_t pairKey(uint32_t flow_id, uint32_t element_id) {
    return (static_cast<uint64_t>(flow_id) << 32u) |
           static_cast<uint64_t>(element_id);
}

double normalizedHash01(uint32_t hash_value) {
    return static_cast<double>(hash_value) /
           static_cast<double>(std::numeric_limits<uint32_t>::max());
}

uint8_t rankFromHash(uint32_t hash_value) {
    uint8_t rank = 1;
    while ((hash_value & 0x80000000u) == 0 && rank < 31) {
        ++rank;
        hash_value <<= 1u;
    }
    return rank;
}

uint64_t memoryBits(double memory_kb) {
    if (memory_kb <= 0.0) {
        return 1;
    }
    return static_cast<uint64_t>(std::max(1.0, memory_kb * 1024.0 * 8.0));
}

class PairWitness {
public:
    struct Update {
        bool positive_candidate = false;
        double delta = 0.0;
    };

    PairWitness(double memory_kb, uint32_t reward, uint32_t penalty)
            : reward_(reward), penalty_(penalty) {
        const uint64_t bits = memoryBits(memory_kb);
        register_num_ = static_cast<uint32_t>(std::max<uint64_t>(1, bits / 5u));
        register_arr_.assign(register_num_, 0);
        event_prob_.fill(0.0);
        event_prob_[0] = 1.0;
        for (uint32_t i = 1; i <= 30; ++i) {
            event_prob_[i] = std::ldexp(1.0, -static_cast<int>(i));
        }
    }

    Update process(uint32_t flow_id, uint32_t element_id) {
        const uint32_t pair[2] = {flow_id, element_id};
        uint32_t register_hash = 0;
        uint32_t rank_hash = 0;
        MurmurHash3_x86_32(pair, static_cast<int>(sizeof(pair)),
                           kRegisterSeed, &register_hash);
        MurmurHash3_x86_32(pair, static_cast<int>(sizeof(pair)),
                           kRankSeed, &rank_hash);

        const uint32_t index = register_hash % register_num_;
        const uint8_t rank = rankFromHash(rank_hash);
        const uint8_t old_rank = register_arr_[index];
        const bool witness = rank > old_rank;
        double delta = -static_cast<double>(penalty_);
        if (witness) {
            const double q = event_probability_;
            if (q <= 0.0) {
                throw std::runtime_error(
                        "PairWitness produced a zero witness probability.");
            }
            delta += static_cast<double>(reward_ + penalty_) / q;
            register_arr_[index] = rank;
            event_probability_ +=
                    (event_prob_[rank] - event_prob_[old_rank]) /
                    static_cast<double>(register_num_);
            event_probability_ = std::max(0.0, event_probability_);
        }
        return {witness, delta};
    }

private:
    uint32_t reward_;
    uint32_t penalty_;
    uint32_t register_num_ = 1;
    std::vector<uint8_t> register_arr_;
    std::array<double, 32> event_prob_{};
    double event_probability_ = 1.0;
};

struct Cell {
    uint32_t flow_id = 0;
    uint32_t frequency = 0;
    double dss = 0.0;
};

class FlowResidence {
public:
    FlowResidence(double memory_kb, uint32_t cell_num,
                  double decay_base, uint32_t reward, uint32_t penalty,
                  ReplacementMechanism mechanism)
            : cell_num_(std::max<uint32_t>(1, cell_num)),
              decay_base_(decay_base),
              reward_(reward),
              penalty_(penalty),
              mechanism_(mechanism) {
        const uint64_t total_bits = memoryBits(memory_kb);
        const uint64_t cell_bits = kFlowIdBits + kFrequencyBits + kDssBits;
        const uint64_t bucket_bits =
                std::max<uint64_t>(1, static_cast<uint64_t>(cell_num_) * cell_bits);
        bucket_num_ = static_cast<uint32_t>(
                std::max<uint64_t>(1, total_bits / bucket_bits));
        cells_.assign(static_cast<size_t>(bucket_num_) * cell_num_, Cell{});
    }

    void update(uint32_t flow_id, uint32_t element_id,
                const PairWitness::Update& pair_update) {
        Cell* bucket = bucketFor(flow_id);
        Cell* empty_cell = nullptr;
        Cell* min_cell = nullptr;
        for (uint32_t i = 0; i < cell_num_; ++i) {
            Cell& cell = bucket[i];
            if (cell.flow_id == flow_id) {
                ++cell.frequency;
                cell.dss += pair_update.delta;
                return;
            }
            if (cell.flow_id == 0) {
                if (empty_cell == nullptr) {
                    empty_cell = &cell;
                }
            } else if (min_cell == nullptr || cell.dss < min_cell->dss) {
                min_cell = &cell;
            }
        }

        if (empty_cell != nullptr) {
            if (!pair_update.positive_candidate) {
                return;
            }
            *empty_cell = {flow_id, 1, pair_update.delta};
            return;
        }

        if (!pair_update.positive_candidate || min_cell == nullptr) {
            return;
        }

        if (mechanism_ == ReplacementMechanism::DirectReplace) {
            *min_cell = {flow_id, 1, pair_update.delta};
            return;
        }

        const double decay = mechanism_ == ReplacementMechanism::DirectDecay
                             ? pair_update.delta
                             : probabilisticDecay(flow_id, element_id,
                                                   *min_cell,
                                                   pair_update.delta);
        min_cell->dss -= decay;
        if (min_cell->dss <= 0.0) {
            *min_cell = {flow_id, 1, pair_update.delta};
        }
    }

    uint64_t residentFlowCount() const {
        uint64_t count = 0;
        for (const Cell& cell : cells_) {
            if (cell.flow_id != 0) {
                ++count;
            }
        }
        return count;
    }

    template <typename Function>
    void forEachRecord(Function&& function) const {
        for (const Cell& cell : cells_) {
            if (cell.flow_id != 0) {
                function(cell.flow_id, cell.frequency, cell.dss);
            }
        }
    }

private:
    Cell* bucketFor(uint32_t flow_id) {
        uint32_t bucket_hash = 0;
        MurmurHash3_x86_32(&flow_id, static_cast<int>(sizeof(flow_id)),
                           kBucketSeed, &bucket_hash);
        return cells_.data() +
               static_cast<size_t>(bucket_hash % bucket_num_) * cell_num_;
    }

    double probabilisticDecay(uint32_t flow_id, uint32_t element_id,
                              const Cell& victim, double event_weight) const {
        const double value =
                std::floor(std::log2(std::max(0.0, victim.dss) + 1.0));
        const double probability = std::pow(decay_base_, -value);
        const uint32_t key[3] = {flow_id, element_id, victim.flow_id};
        uint32_t hash_value = 0;
        MurmurHash3_x86_32(key, static_cast<int>(sizeof(key)), kDecaySeed,
                           &hash_value);
        return normalizedHash01(hash_value) <= probability ? event_weight : 0.0;
    }

    uint32_t cell_num_;
    uint32_t bucket_num_ = 1;
    double decay_base_;
    uint32_t reward_;
    uint32_t penalty_;
    ReplacementMechanism mechanism_;
    std::vector<Cell> cells_;
};

struct EstimateMap {
    std::unordered_map<uint32_t, uint32_t> frequency;
    std::unordered_map<uint32_t, double> dss;
};

} // namespace

ReplacementAblation::ReplacementAblation(
        const ReplacementAblationConfig& config)
        : config_(config) {
    if (config_.reward == 0 || config_.penalty == 0) {
        throw std::invalid_argument("reward and penalty must be positive");
    }
    if (config_.cell_num == 0) {
        throw std::invalid_argument("cell_num must be positive");
    }
    if (config_.decay_base <= 1.0) {
        throw std::invalid_argument("decay_base must be greater than 1");
    }
    if (config_.total_memory_kb <= 0.0) {
        throw std::invalid_argument("total_memory_kb must be positive");
    }
    if (config_.pair_memory_ratio < 0.0 ||
        config_.residence_memory_ratio < 0.0 ||
        config_.pair_memory_ratio + config_.residence_memory_ratio > 1.0 + 1e-12) {
        throw std::invalid_argument(
                "memory ratios must be non-negative and sum to at most 1");
    }
}

const char* ReplacementAblation::mechanismName(
        ReplacementMechanism mechanism) {
    switch (mechanism) {
        case ReplacementMechanism::ProbabilisticDecay:
            return "DS-ProbDecay";
        case ReplacementMechanism::DirectDecay:
            return "DS-DirectDecay";
        case ReplacementMechanism::DirectReplace:
            return "DS-DirectReplace";
    }
    return "Unknown";
}

std::vector<AblationResult> ReplacementAblation::run(
        const std::vector<std::pair<uint32_t, uint32_t>>& data,
        const ExactStats& exact_stats,
        const std::vector<float>& tau_f_ratio_list,
        const std::vector<double>& tau_s_list,
        ReplacementMechanism mechanism) const {
    if (tau_f_ratio_list.empty() || tau_s_list.empty()) {
        return {};
    }

    PairWitness pair_witness(
            config_.total_memory_kb * config_.pair_memory_ratio,
            config_.reward, config_.penalty);
    FlowResidence residence(
            config_.total_memory_kb * config_.residence_memory_ratio,
            config_.cell_num, config_.decay_base, config_.reward,
            config_.penalty, mechanism);

    const auto update_start = std::chrono::steady_clock::now();
    for (const auto& item : data) {
        const PairWitness::Update pair_update =
                pair_witness.process(item.first, item.second);
        residence.update(item.first, item.second, pair_update);
    }
    const auto update_end = std::chrono::steady_clock::now();
    const double update_seconds =
            std::chrono::duration<double>(update_end - update_start).count();
    const double update_mips = update_seconds <= 0.0
                               ? 0.0
                               : static_cast<double>(data.size()) /
                                 update_seconds / 1000000.0;

    EstimateMap estimates;
    estimates.frequency.reserve(
            static_cast<size_t>(residence.residentFlowCount()));
    estimates.dss.reserve(static_cast<size_t>(residence.residentFlowCount()));
    residence.forEachRecord([&](uint32_t flow_id, uint32_t frequency,
                                double dss) {
        estimates.frequency[flow_id] = frequency;
        estimates.dss[flow_id] = dss;
    });

    std::vector<AblationResult> results;
    for (float tau_f_ratio : tau_f_ratio_list) {
        for (double tau_s : tau_s_list) {
            const TruthSnapshot truth = buildTruthSnapshot(
                    exact_stats, tau_f_ratio, config_.reward, config_.penalty,
                    tau_s);
            uint64_t tp = 0;
            uint64_t fp = 0;
            uint64_t predicted = 0;
            residence.forEachRecord([&](uint32_t flow_id, uint32_t frequency,
                                        double dss) {
                if (frequency < truth.tau_f || dss < tau_s) {
                    return;
                }
                ++predicted;
                if (truth.contains(flow_id)) {
                    ++tp;
                } else {
                    ++fp;
                }
            });

            const uint64_t truth_count = truth.truth_count;
            const uint64_t fn = truth_count - tp;
            double dss_are_sum = 0.0;
            double frequency_are_sum = 0.0;
            for (uint32_t flow_id : truth.diffuse_flows) {
                const uint32_t true_frequency =
                        exact_stats.trueFrequency(flow_id);
                const uint32_t true_cardinality =
                        exact_stats.trueCardinality(flow_id);
                const double true_dss = computeDiffuseSupport(
                        true_frequency, true_cardinality, config_.reward,
                        config_.penalty);
                const auto freq_it = estimates.frequency.find(flow_id);
                const auto dss_it = estimates.dss.find(flow_id);
                const double estimated_frequency =
                        freq_it == estimates.frequency.end() ? 0.0
                                                             : freq_it->second;
                const double estimated_dss =
                        dss_it == estimates.dss.end() ? 0.0 : dss_it->second;
                const double frequency_error =
                        std::abs(estimated_frequency - true_frequency);
                const double dss_error = std::abs(estimated_dss - true_dss);
                frequency_are_sum +=
                        frequency_error / static_cast<double>(true_frequency);
                if (true_dss != 0.0) {
                    dss_are_sum += dss_error / std::abs(true_dss);
                }
            }

            const double precision = tp + fp == 0
                                     ? 0.0
                                     : static_cast<double>(tp) /
                                       static_cast<double>(tp + fp);
            const double recall = truth_count == 0
                                  ? 0.0
                                  : static_cast<double>(tp) /
                                    static_cast<double>(truth_count);
            const double f1 = precision + recall == 0.0
                              ? 0.0
                              : 2.0 * precision * recall /
                                (precision + recall);
            const double denominator = static_cast<double>(truth_count);

            AblationResult result;
            result.variant = mechanismName(mechanism);
            result.total_memory_kb = config_.total_memory_kb;
            result.pair_memory_kb =
                    config_.total_memory_kb * config_.pair_memory_ratio;
            result.pair_memory_limited = true;
            result.residence_memory_kb =
                    config_.total_memory_kb * config_.residence_memory_ratio;
            result.residence_memory_limited = true;
            result.cell_num = config_.cell_num;
            result.tau_f = truth.tau_f;
            result.tau_s = tau_s;
            result.resident_flows = residence.residentFlowCount();
            result.ground_truth_flows = truth_count;
            result.predicted_flows = predicted;
            result.tp = tp;
            result.fp = fp;
            result.fn = fn;
            result.precision = precision;
            result.recall = recall;
            result.f1 = f1;
            result.dss_are = denominator == 0.0
                             ? 0.0
                             : dss_are_sum / denominator;
            result.frequency_are = denominator == 0.0
                                   ? 0.0
                                   : frequency_are_sum / denominator;
            result.update_mips = update_mips;
            results.push_back(result);
        }
    }
    return results;
}

