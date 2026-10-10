#include "header/Krill.h"

#include "header/DatasetLoader.h"
#include "header/ExactStats.h"
#include "header/MurmurHash3.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace {

double normalizedHash01(uint32_t hash_value) {
    return static_cast<double>(hash_value) /
           static_cast<double>(std::numeric_limits<uint32_t>::max());
}

TruthSnapshot makeTruthSnapshot(const ExactStats& exact_stats,
                                float tau_f_ratio,
                                uint32_t reward,
                                uint32_t penalty,
                                double tau_s) {
    const TruthSnapshot truth = buildTruthSnapshot(exact_stats,
                                                   tau_f_ratio,
                                                   reward,
                                                   penalty,
                                                   tau_s);
    std::cout << "Truth: tau_f=" << truth.tau_f
              << ", tau_S=" << truth.tau_s
              << ", reward=" << truth.reward
              << ", penalty=" << truth.penalty
              << ", truth=" << truth.truth_count
              << std::endl;
    return truth;
}

} // namespace

uint8_t Krill::rankFromHash(uint32_t hash_value) {
    uint8_t rank = 1;
    while ((hash_value & 0x80000000u) == 0 && rank < 31) {
        ++rank;
        hash_value <<= 1u;
    }
    return rank;
}

Krill::Krill(float input_memory_kb)
        : Krill(input_memory_kb, 0.3f, 4, 12, 1) {}

Krill::Krill(float input_memory_kb,
                             float input_pair_witness_ratio,
                             uint32_t input_cell_num,
                             uint32_t input_reward,
                             uint32_t input_penalty)
        : memory_kb(input_memory_kb),
          pair_witness_ratio(input_pair_witness_ratio),
          cell_num(input_cell_num),
          reward(input_reward),
          penalty(input_penalty) {
    if (cell_num == 0) {
        throw std::invalid_argument("cell_num must be positive.");
    }
    if (pair_witness_ratio <= 0.0f || pair_witness_ratio >= 1.0f) {
        throw std::invalid_argument("pair_witness_ratio must be in (0, 1).");
    }
    if (reward == 0 || penalty == 0) {
        throw std::invalid_argument("reward and penalty must be positive.");
    }

    event_prob[0] = 1.0;
    for (uint32_t i = 1; i <= 30; ++i) {
        event_prob[i] = std::ldexp(1.0, -static_cast<int>(i));
    }
    event_prob[31] = 0.0;

    const uint64_t total_bits = static_cast<uint64_t>(
            std::max<double>(1.0, static_cast<double>(memory_kb) * 1024.0 * 8.0));
    uint64_t pair_witness_bits = static_cast<uint64_t>(
            static_cast<double>(total_bits) * static_cast<double>(pair_witness_ratio));
    pair_witness_bits = std::max<uint64_t>(5, pair_witness_bits);
    if (pair_witness_bits >= total_bits) {
        pair_witness_bits = total_bits / 2;
    }
    register_num = static_cast<uint32_t>(std::max<uint64_t>(1, pair_witness_bits / 5u));
    pair_witness_bits = static_cast<uint64_t>(register_num) * 5u;

    const uint64_t sketch_bits = std::max<uint64_t>(1, total_bits - pair_witness_bits);
    const uint64_t cell_bits = static_cast<uint64_t>(flow_id_bits) +
                               static_cast<uint64_t>(freq_bits) +
                               static_cast<uint64_t>(dss_bits);
    const uint64_t bucket_bits =
            std::max<uint64_t>(1, static_cast<uint64_t>(cell_num) * cell_bits);
    bucket_num = static_cast<uint32_t>(std::max<uint64_t>(1, sketch_bits / bucket_bits));

    register_arr.assign(register_num, 0);
    cells.assign(static_cast<size_t>(bucket_num) * static_cast<size_t>(cell_num), Cell{});
}

double Krill::exponentialProbabilisticDecay(uint32_t flow_id,
                                                    uint32_t element_id,
                                                    const Cell& min_cell,
                                                    double event_weight) {
    constexpr double decay_base = 1.08;
    const double val =
            std::floor(std::log2(std::max(0.0, min_cell.dss) + 1.0));
    const double probability = std::pow(decay_base, -val);
    const uint32_t key[3] = {flow_id, element_id, min_cell.flow_id};
    uint32_t hash_value = 0;
    MurmurHash3_x86_32(key, static_cast<int>(sizeof(key)), 6001, &hash_value);
    return normalizedHash01(hash_value) <= probability ? event_weight : 0.0;
}

void Krill::update(uint32_t flow_id, uint32_t element_id) {
    constexpr uint32_t register_seed = 2001;
    constexpr uint32_t rank_seed = 3001;
    constexpr uint32_t bucket_seed = 1001;

    const uint32_t pair_key[2] = {flow_id, element_id};
    uint32_t register_hash = 0;
    uint32_t rank_hash = 0;
    MurmurHash3_x86_32(pair_key, sizeof(pair_key), register_seed,
                       &register_hash);
    MurmurHash3_x86_32(pair_key, sizeof(pair_key), rank_seed,
                       &rank_hash);

    const uint32_t register_index = register_hash % register_num;
    const uint8_t rank = rankFromHash(rank_hash);
    const uint8_t old = register_arr[register_index];
    const bool event = rank > old;

    double q_event = 0.0;
    double new_pair_compensation = 0.0;
    if (event) {
        q_event = event_probability;
        new_pair_compensation =
                static_cast<double>(reward + penalty) /
                q_event;
        register_arr[register_index] = rank;
        event_probability +=
                (event_prob[rank] - event_prob[old]) /
                static_cast<double>(register_num);
        event_probability = std::max(0.0, event_probability);
        ++state_change_events;
        compensation_sum += new_pair_compensation;
    }

    uint32_t bucket_hash = 0;
    MurmurHash3_x86_32(&flow_id, sizeof(flow_id), bucket_seed, &bucket_hash);
    const uint32_t local_cell_num = cell_num;
    Cell* bucket = cells.data() +
                   static_cast<size_t>(bucket_hash % bucket_num) *
                   static_cast<size_t>(local_cell_num);

    Cell* empty_cell = nullptr;
    Cell* min_cell = nullptr;
    for (uint32_t i = 0; i < local_cell_num; ++i) {
        Cell& cell = bucket[i];
        if (cell.flow_id == flow_id) {
            ++cell.freq;
            cell.dss -= static_cast<double>(penalty);
            if (event) {
                cell.dss += new_pair_compensation;
            }
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

    if (!event) {
        return;
    }

    const double init_dss =
            -static_cast<double>(penalty) + new_pair_compensation;
    if (empty_cell != nullptr) {
        empty_cell->flow_id = flow_id;
        empty_cell->freq = 1;
        empty_cell->dss = init_dss;
        return;
    }

    if (min_cell == nullptr) {
        return;
    }

    const double decay =
            exponentialProbabilisticDecay(flow_id, element_id, *min_cell, init_dss);
    min_cell->dss -= decay;
    if (min_cell->dss <= 0.0) {
        min_cell->flow_id = flow_id;
        min_cell->freq = 1;
        min_cell->dss = init_dss;
    }
}

std::vector<FlowEstimate> Krill::collectEstimatesForEvaluation(
        const ExactStats& exact_stats) {
    std::vector<FlowEstimate> estimates;
    estimates.reserve(cells.size());
    for (const Cell& cell : cells) {
        if (cell.flow_id == 0 || exact_stats.trueFrequency(cell.flow_id) == 0) {
            continue;
        }
        estimates.push_back({cell.flow_id,
                             cell.freq,
                             cell.dss});
    }
    return estimates;
}

void Krill::eval(const std::vector<std::pair<uint32_t, uint32_t>>& data,
                         const ExactStats& exact_stats,
                         const std::vector<float>& tau_f_ratio_list,
                         const std::vector<double>& tau_s_list) {

    if (tau_f_ratio_list.empty()) {
        std::cerr << "eval skipped: frequency threshold cannot be empty." << std::endl;
        return;
    }
    if (tau_s_list.empty()) {
        std::cerr << "eval skipped: score threshold cannot be empty." << std::endl;
        return;
    }

    const auto update_start = std::chrono::steady_clock::now();
    for (const auto& item : data) {
        update(item.first, item.second);
    }
    const auto update_end = std::chrono::steady_clock::now();
    const double update_seconds =
            std::chrono::duration<double>(update_end - update_start).count();
    const double update_mips = update_seconds <= 0.0
                               ? 0.0
                               : static_cast<double>(data.size()) /
                                 update_seconds / 1000000.0;

    struct ResidentStat {
        uint32_t flow_id = 0;
        uint32_t est_f = 0;
        double est_s = 0.0;
        uint32_t true_f = 0;
    };

    std::vector<ResidentStat> resident_stats;
    resident_stats.reserve(cells.size());
    std::unordered_map<uint32_t, uint32_t> estimated_frequency;
    estimated_frequency.reserve(cells.size());
    std::unordered_map<uint32_t, double> estimated_dss;
    estimated_dss.reserve(cells.size());
    uint32_t resident_count = 0;
    for (const auto& cell : cells) {
        if (cell.flow_id == 0) {
            continue;
        }
        ++resident_count;
        const uint32_t true_f = exact_stats.trueFrequency(cell.flow_id);
        if (true_f == 0) {
            continue;
        }
        resident_stats.push_back({cell.flow_id, cell.freq, cell.dss, true_f});
        estimated_frequency[cell.flow_id] = cell.freq;
        estimated_dss[cell.flow_id] = cell.dss;
    }

    for (float tau_f_ratio : tau_f_ratio_list) {
        for (double tau_s : tau_s_list) {
            const TruthSnapshot truth_snapshot = buildTruthSnapshot(exact_stats,
                                                                    tau_f_ratio,
                                                                    reward,
                                                                    penalty,
                                                                    tau_s);
            const uint32_t tau_f = truth_snapshot.tau_f;
            uint32_t tp = 0;
            uint32_t fp = 0;
            uint32_t predicted = 0;
            uint32_t truth_resident_total = 0;
            uint32_t fn_frequency_loss = 0;
            uint32_t fn_dss_loss = 0;
            double dss_are_sum = 0.0;
            for (const auto& cell : resident_stats) {
                const bool detected = cell.est_f >= tau_f &&
                                      cell.est_s >= tau_s;
                const bool truth = truth_snapshot.contains(cell.flow_id);
                if (truth) {
                    ++truth_resident_total;
                }
                if (detected) {
                    ++predicted;
                    if (truth) {
                        ++tp;
                    } else {
                        ++fp;
                    }
                } else if (truth) {
                    if (cell.est_f < tau_f) {
                        ++fn_frequency_loss;
                    } else if (cell.est_s < tau_s) {
                        ++fn_dss_loss;
                    }
                }
            }
            const uint32_t truth_count = truth_snapshot.truth_count;
            const uint32_t fn = truth_count - tp;
            double frequency_are_sum = 0.0;
            for (uint32_t flow_id : truth_snapshot.diffuse_flows) {
                const uint32_t true_f = exact_stats.trueFrequency(flow_id);
                const uint32_t true_c = exact_stats.trueCardinality(flow_id);
                const auto est_iter = estimated_frequency.find(flow_id);
                const uint32_t est_f =
                        est_iter == estimated_frequency.end() ? 0 : est_iter->second;
                const double abs_error = std::abs(static_cast<double>(est_f) -
                                                  static_cast<double>(true_f));
                frequency_are_sum += abs_error / static_cast<double>(true_f);

                const double true_dss = computeDiffuseSupport(true_f, true_c, reward, penalty);
                const auto dss_iter = estimated_dss.find(flow_id);
                const double est_dss = dss_iter == estimated_dss.end() ? 0.0 : dss_iter->second;
                const double dss_abs_error = std::abs(est_dss - true_dss);
                dss_are_sum += true_dss == 0.0 ? 0.0 : dss_abs_error / std::abs(true_dss);
            }
            const float precision = tp + fp == 0
                                    ? 0.0f
                                    : static_cast<float>(tp) / static_cast<float>(tp + fp);
            const float recall = tp + fn == 0
                                 ? 0.0f
                                 : static_cast<float>(tp) / static_cast<float>(tp + fn);
            const float f1 = precision + recall == 0.0f
                             ? 0.0f
                             : 2.0f * precision * recall / (precision + recall);

            KrillEvalResult result;
            result.memory_kb = memory_kb;
            result.pair_witness_ratio = pair_witness_ratio;
            result.cell_num = cell_num;
            result.reward = reward;
            result.penalty = penalty;
            result.register_num = register_num;
            result.bucket_num = bucket_num;
            result.final_q = event_probability;
            result.state_change_events = static_cast<double>(state_change_events);
            result.avg_compensation = state_change_events == 0
                                      ? 0.0
                                      : compensation_sum /
                                        static_cast<double>(state_change_events);
            result.tau_f_ratio = tau_f_ratio;
            result.tau_f = tau_f;
            result.diffuse_threshold = tau_s;
            result.resident_flows = resident_count;
            result.ground_truth_count = truth_count;
            result.predicted_count = predicted;
            result.tp = tp;
            result.fp = fp;
            result.fn = fn;
            result.precision = precision;
            result.recall = recall;
            result.f1 = f1;
            result.frequency_are = truth_count == 0
                                   ? 0.0
                                   : frequency_are_sum / static_cast<double>(truth_count);
result.dss_are = truth_count == 0
                              ? 0.0
                              : dss_are_sum / static_cast<double>(truth_count);
result.fn_not_resident = truth_count - truth_resident_total;
            result.fn_resident = truth_resident_total - tp;
            result.fn_frequency_loss = fn_frequency_loss;
            result.fn_dss_loss = fn_dss_loss;
            result.update_mips = update_mips;

            std::cout << std::fixed << std::setprecision(6)
                      << " - Krill: memory_kb=" << result.memory_kb
                      << ", pair_witness_ratio=" << result.pair_witness_ratio
                      << ", reward=" << result.reward
                      << ", penalty=" << result.penalty
                      << ", tau_f_ratio=" << result.tau_f_ratio
                      << ", tau_f=" << result.tau_f
                      << ", tau_S=" << result.diffuse_threshold
                      << ", truth=" << result.ground_truth_count
                      << ", predicted=" << result.predicted_count
                      << ", F1=" << result.f1
                      << ", precision=" << result.precision
                      << ", recall=" << result.recall
                      << ", dss_ARE=" << result.dss_are
                      << ", freq_ARE=" << result.frequency_are
                      << ", throughput_Mips=" << result.update_mips
                      << std::endl;
        }
    }
}

KrillEvalResult Krill::evaluateForOrder(
        const std::vector<std::pair<uint32_t, uint32_t>>& data,
        const ExactStats& exact_stats,
        float tau_f_ratio,
        double tau_s) {
    const auto update_start = std::chrono::steady_clock::now();
    for (const auto& item : data) {
        update(item.first, item.second);
    }
    const auto update_end = std::chrono::steady_clock::now();
    const double update_seconds =
            std::chrono::duration<double>(update_end - update_start).count();
    const double update_mips = update_seconds <= 0.0
                               ? 0.0
                               : static_cast<double>(data.size()) /
                                 update_seconds / 1000000.0;

    struct ResidentStat {
        uint32_t flow_id = 0;
        uint32_t est_f = 0;
        double est_s = 0.0;
    };

    std::vector<ResidentStat> resident_stats;
    resident_stats.reserve(cells.size());
    std::unordered_map<uint32_t, uint32_t> estimated_frequency;
    estimated_frequency.reserve(cells.size());
    std::unordered_map<uint32_t, double> estimated_dss;
    estimated_dss.reserve(cells.size());
    uint32_t resident_count = 0;
    for (const auto& cell : cells) {
        if (cell.flow_id == 0) {
            continue;
        }
        ++resident_count;
        if (exact_stats.trueFrequency(cell.flow_id) == 0) {
            continue;
        }
        resident_stats.push_back({cell.flow_id, cell.freq, cell.dss});
        estimated_frequency[cell.flow_id] = cell.freq;
        estimated_dss[cell.flow_id] = cell.dss;
    }

    const TruthSnapshot truth_snapshot = buildTruthSnapshot(exact_stats,
                                                            tau_f_ratio,
                                                            reward,
                                                            penalty,
                                                            tau_s);
    const uint32_t tau_f = truth_snapshot.tau_f;
    uint32_t tp = 0;
    uint32_t fp = 0;
    uint32_t predicted = 0;
    uint32_t truth_resident_total = 0;
    uint32_t fn_frequency_loss = 0;
    uint32_t fn_dss_loss = 0;
    double frequency_are_sum = 0.0;
    double dss_are_sum = 0.0;

    for (const auto& cell : resident_stats) {
        const bool truth = truth_snapshot.contains(cell.flow_id);
        const bool detected = cell.est_f >= tau_f && cell.est_s >= tau_s;
        if (truth) {
            ++truth_resident_total;
        }
        if (detected) {
            ++predicted;
            if (truth) {
                ++tp;
            } else {
                ++fp;
            }
        } else if (truth) {
            if (cell.est_f < tau_f) {
                ++fn_frequency_loss;
            } else if (cell.est_s < tau_s) {
                ++fn_dss_loss;
            }
        }
    }

    for (uint32_t flow_id : truth_snapshot.diffuse_flows) {
        const uint32_t true_f = exact_stats.trueFrequency(flow_id);
        const uint32_t true_c = exact_stats.trueCardinality(flow_id);
        const auto frequency_iter = estimated_frequency.find(flow_id);
        const uint32_t est_f = frequency_iter == estimated_frequency.end()
                               ? 0
                               : frequency_iter->second;
        const double frequency_abs_error =
                std::abs(static_cast<double>(est_f) - static_cast<double>(true_f));
        frequency_are_sum += frequency_abs_error / static_cast<double>(true_f);

        const double true_dss = computeDiffuseSupport(true_f,
                                                      true_c,
                                                      reward,
                                                      penalty);
        const auto dss_iter = estimated_dss.find(flow_id);
        const double est_dss = dss_iter == estimated_dss.end()
                               ? 0.0
                               : dss_iter->second;
        const double dss_abs_error = std::abs(est_dss - true_dss);
        dss_are_sum += true_dss == 0.0
                       ? 0.0
                       : dss_abs_error / std::abs(true_dss);
    }

    const uint32_t fn = truth_snapshot.truth_count - tp;
    const float precision = tp + fp == 0
                            ? 0.0f
                            : static_cast<float>(tp) /
                              static_cast<float>(tp + fp);
    const float recall = tp + fn == 0
                         ? 0.0f
                         : static_cast<float>(tp) /
                           static_cast<float>(tp + fn);
    const float f1 = precision + recall == 0.0f
                     ? 0.0f
                     : 2.0f * precision * recall / (precision + recall);

    KrillEvalResult result;
    result.memory_kb = memory_kb;
    result.pair_witness_ratio = pair_witness_ratio;
    result.cell_num = cell_num;
    result.reward = reward;
    result.penalty = penalty;
    result.register_num = register_num;
    result.bucket_num = bucket_num;
    result.final_q = event_probability;
    result.state_change_events = static_cast<double>(state_change_events);
    result.avg_compensation = state_change_events == 0
                              ? 0.0
                              : compensation_sum /
                                static_cast<double>(state_change_events);
    result.tau_f_ratio = tau_f_ratio;
    result.tau_f = tau_f;
    result.diffuse_threshold = tau_s;
    result.resident_flows = resident_count;
    result.ground_truth_count = truth_snapshot.truth_count;
    result.predicted_count = predicted;
    result.tp = tp;
    result.fp = fp;
    result.fn = fn;
    result.precision = precision;
    result.recall = recall;
    result.f1 = f1;
    result.frequency_are = truth_snapshot.truth_count == 0
                           ? 0.0
                           : frequency_are_sum /
                             static_cast<double>(truth_snapshot.truth_count);
result.dss_are = truth_snapshot.truth_count == 0
                     ? 0.0
                     : dss_are_sum /
                       static_cast<double>(truth_snapshot.truth_count);
result.fn_not_resident = truth_snapshot.truth_count - truth_resident_total;
    result.fn_resident = truth_resident_total - tp;
    result.fn_frequency_loss = fn_frequency_loss;
    result.fn_dss_loss = fn_dss_loss;
    result.update_mips = update_mips;
    return result;
}

void Krill::evalForSweep(
        const std::vector<std::pair<uint32_t, uint32_t>>& data,
        const ExactStats& exact_stats,
        const std::vector<float>& tau_f_ratio_list,
        const std::vector<double>& diffuse_threshold_list) {
    const auto update_start = std::chrono::steady_clock::now();
    for (const auto& item : data) {
        update(item.first, item.second);
    }
    const auto update_end = std::chrono::steady_clock::now();
    const double update_seconds =
            std::chrono::duration<double>(update_end - update_start).count();
    const double update_mips = update_seconds <= 0.0
                               ? 0.0
                               : static_cast<double>(data.size()) /
                                 update_seconds / 1000000.0;

    for (float tau_f_ratio : tau_f_ratio_list) {
        for (double diffuse_threshold : diffuse_threshold_list) {
            const TruthSnapshot truth_snapshot =
                    makeTruthSnapshot(exact_stats,
                                      tau_f_ratio,
                                      reward,
                                      penalty,
                                      diffuse_threshold);
            const uint32_t tau_f = truth_snapshot.tau_f;

            struct ResidentStat {
                uint32_t flow_id = 0;
                uint32_t est_f = 0;
                double est_s = 0.0;
                uint32_t true_f = 0;
            };

            std::vector<ResidentStat> resident_stats;
            resident_stats.reserve(cells.size());
            std::unordered_map<uint32_t, uint32_t> estimated_frequency;
            estimated_frequency.reserve(cells.size());
            std::unordered_map<uint32_t, double> estimated_dss;
            estimated_dss.reserve(cells.size());

            uint32_t resident_count = 0;
            for (const auto& cell : cells) {
                if (cell.flow_id == 0) {
                    continue;
                }
                ++resident_count;
                const uint32_t true_f = exact_stats.trueFrequency(cell.flow_id);
                if (true_f == 0) {
                    continue;
                }
                resident_stats.push_back({cell.flow_id, cell.freq, cell.dss, true_f});
                estimated_frequency[cell.flow_id] = cell.freq;
                estimated_dss[cell.flow_id] = cell.dss;
            }

            uint32_t tp = 0;
            uint32_t fp = 0;
            uint32_t predicted = 0;
            uint32_t truth_resident_total = 0;
            uint32_t fn_frequency_loss = 0;
            uint32_t fn_dss_loss = 0;
            double dss_are_sum = 0.0;
            for (const auto& cell : resident_stats) {
                const bool truth = truth_snapshot.contains(cell.flow_id);
                const bool detected = cell.est_f >= tau_f &&
                                      cell.est_s >= diffuse_threshold;
                if (truth) {
                    ++truth_resident_total;
                }
                if (detected) {
                    ++predicted;
                    if (truth) {
                        ++tp;
                    } else {
                        ++fp;
                    }
                } else if (truth) {
                    if (cell.est_f < tau_f) {
                        ++fn_frequency_loss;
                    } else if (cell.est_s < diffuse_threshold) {
                        ++fn_dss_loss;
                    }
                }
            }

            const uint32_t truth_count = truth_snapshot.truth_count;
            const uint32_t fn = truth_count - tp;
            double frequency_are_sum = 0.0;
            for (uint32_t flow_id : truth_snapshot.diffuse_flows) {
                const uint32_t true_f = exact_stats.trueFrequency(flow_id);
                const uint32_t true_c = exact_stats.trueCardinality(flow_id);
                const auto est_iter = estimated_frequency.find(flow_id);
                const uint32_t est_f =
                        est_iter == estimated_frequency.end() ? 0 : est_iter->second;
                const double abs_error = std::abs(static_cast<double>(est_f) -
                                                  static_cast<double>(true_f));
                frequency_are_sum += abs_error / static_cast<double>(true_f);

                const double true_dss = computeDiffuseSupport(true_f, true_c, reward, penalty);
                const auto dss_iter = estimated_dss.find(flow_id);
                const double est_dss = dss_iter == estimated_dss.end() ? 0.0 : dss_iter->second;
                const double dss_abs_error = std::abs(est_dss - true_dss);
                dss_are_sum += true_dss == 0.0 ? 0.0 : dss_abs_error / std::abs(true_dss);
            }

            const float precision = tp + fp == 0
                                    ? 0.0f
                                    : static_cast<float>(tp) /
                                      static_cast<float>(tp + fp);
            const float recall = tp + fn == 0
                                 ? 0.0f
                                 : static_cast<float>(tp) /
                                   static_cast<float>(tp + fn);
            const float f1 = precision + recall == 0.0f
                             ? 0.0f
                             : 2.0f * precision * recall / (precision + recall);

            KrillEvalResult result;
            result.memory_kb = memory_kb;
            result.pair_witness_ratio = pair_witness_ratio;
            result.cell_num = cell_num;
            result.reward = reward;
            result.penalty = penalty;
            result.register_num = register_num;
            result.bucket_num = bucket_num;
            result.final_q = event_probability;
            result.state_change_events = static_cast<double>(state_change_events);
            result.avg_compensation = state_change_events == 0
                                      ? 0.0
                                      : compensation_sum /
                                        static_cast<double>(state_change_events);
            result.tau_f_ratio = tau_f_ratio;
            result.tau_f = tau_f;
            result.diffuse_threshold = diffuse_threshold;
            result.resident_flows = resident_count;
            result.ground_truth_count = truth_count;
            result.predicted_count = predicted;
            result.tp = tp;
            result.fp = fp;
            result.fn = fn;
            result.precision = precision;
            result.recall = recall;
            result.f1 = f1;
            result.frequency_are = truth_count == 0
                                    ? 0.0
                                    : frequency_are_sum / static_cast<double>(truth_count);
result.dss_are = truth_count == 0
                              ? 0.0
                              : dss_are_sum / static_cast<double>(truth_count);
result.fn_not_resident = truth_count - truth_resident_total;
            result.fn_resident = truth_resident_total - tp;
            result.fn_frequency_loss = fn_frequency_loss;
            result.fn_dss_loss = fn_dss_loss;
            result.update_mips = update_mips;

            std::cout << std::fixed << std::setprecision(6)
                      << " - KrillSweep: memory_kb=" << result.memory_kb
                      << ", pair_witness_ratio=" << result.pair_witness_ratio
                      << ", reward=" << result.reward
                      << ", penalty=" << result.penalty
                      << ", cell_num=" << result.cell_num
                      << ", tau_f_ratio=" << result.tau_f_ratio
                      << ", tau_f=" << result.tau_f
                      << ", tau_S=" << result.diffuse_threshold
                      << ", truth=" << result.ground_truth_count
                      << ", predicted=" << result.predicted_count
                      << ", F1=" << result.f1
                      << ", precision=" << result.precision
                      << ", recall=" << result.recall
                      << ", dss_ARE=" << result.dss_are
                      << ", freq_ARE=" << result.frequency_are
                      << ", throughput_Mips=" << result.update_mips
                      << std::endl;
        }
    }
}



