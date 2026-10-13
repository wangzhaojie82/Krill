#include "header/CocoSketch.h"

#include "header/ExactStats.h"
#include "header/MurmurHash3.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <unordered_map>

namespace {

constexpr uint32_t kBucketBits = 64 + 32;

} // namespace

CocoSketch::CocoSketch(float input_memory_kb)
        : memory_kb(input_memory_kb),
          rng(2001),
          probability(0.0, 1.0) {
    const double total_bits =
            std::max(1.0, static_cast<double>(memory_kb) * 1024.0 * 8.0);
    width = static_cast<uint32_t>(
            std::max(1.0, std::floor(total_bits /
                                     (static_cast<double>(depth) *
                                      static_cast<double>(kBucketBits)))));
    table.assign(depth, std::vector<Bucket>(width));
    seeds.reserve(depth);
    for (uint32_t i = 0; i < depth; ++i) {
        seeds.push_back(0x9e3779b9u + i * 0x85ebca6bu);
    }
}

uint64_t CocoSketch::combineKey(uint32_t flow_id, uint32_t element_id) {
    return (static_cast<uint64_t>(flow_id) << 32u) |
           static_cast<uint64_t>(element_id);
}

std::pair<uint32_t, uint32_t> CocoSketch::splitKey(uint64_t key) {
    return {static_cast<uint32_t>(key >> 32u), static_cast<uint32_t>(key)};
}

uint32_t CocoSketch::bucketIndex(uint64_t key, uint32_t row) const {
    uint32_t hash_value = 0;
    MurmurHash3_x86_32(&key, sizeof(key), seeds[row], &hash_value);
    return hash_value % width;
}

void CocoSketch::update(uint32_t flow_id, uint32_t element_id) {
    if (flow_id == 0 || element_id == 0) {
        return;
    }

    const uint64_t key = combineKey(flow_id, element_id);
    int min_row = -1;
    uint32_t min_value = std::numeric_limits<uint32_t>::max();
    int tie_count = 0;
    int first_empty_row = -1;
    std::vector<uint32_t> indices(depth);

    for (uint32_t row = 0; row < depth; ++row) {
        const uint32_t index = bucketIndex(key, row);
        indices[row] = index;
        Bucket& bucket = table[row][index];

        if (bucket.key == key) {
            ++bucket.value;
            return;
        }

        if (first_empty_row == -1 && bucket.key == 0) {
            first_empty_row = static_cast<int>(row);
        }

        if (first_empty_row == -1) {
            if (bucket.value < min_value) {
                min_value = bucket.value;
                min_row = static_cast<int>(row);
                tie_count = 1;
            } else if (bucket.value == min_value) {
                ++tie_count;
                std::uniform_int_distribution<int> tie_distribution(0, tie_count - 1);
                if (tie_distribution(rng) == 0) {
                    min_row = static_cast<int>(row);
                }
            }
        }
    }

    if (first_empty_row != -1) {
        Bucket& target = table[first_empty_row][indices[first_empty_row]];
        target.key = key;
        target.value = 1;
        return;
    }

    Bucket& target = table[min_row][indices[min_row]];
    ++target.value;
    const double replacement_probability = 1.0 / static_cast<double>(target.value);
    if (probability(rng) < replacement_probability) {
        target.key = key;
    }
}

std::vector<FlowEstimate> CocoSketch::collectEstimatesForEvaluation(
        const ExactStats& exact_stats) const {
    std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> flow_elements;
    flow_elements.reserve(table.size() * width);

    for (const auto& row : table) {
        for (const Bucket& bucket : row) {
            if (bucket.key == 0 || bucket.value == 0) {
                continue;
            }
            const auto [flow_id, element_id] = splitKey(bucket.key);
            if (flow_id == 0 || element_id == 0 ||
                exact_stats.trueFrequency(flow_id) == 0) {
                continue;
            }
            flow_elements[flow_id][element_id] += bucket.value;
        }
    }

    std::vector<FlowEstimate> estimates;
    estimates.reserve(flow_elements.size());
    for (const auto& flow_entry : flow_elements) {
        uint64_t frequency_sum = 0;
        for (const auto& element_entry : flow_entry.second) {
            frequency_sum += element_entry.second;
        }
        const uint32_t frequency = frequency_sum >
                                   static_cast<uint64_t>(std::numeric_limits<uint32_t>::max())
                                   ? std::numeric_limits<uint32_t>::max()
                                   : static_cast<uint32_t>(frequency_sum);
        const uint32_t cardinality =
                static_cast<uint32_t>(flow_entry.second.size());
        const double dss = computeDiffuseSupport(frequency,
                                                 cardinality,
                                                 reward,
                                                 penalty);
        estimates.push_back({flow_entry.first, frequency, dss});
    }
    return estimates;
}

void CocoSketch::eval(
        const std::vector<std::pair<uint32_t, uint32_t>>& data,
        const ExactStats& exact_stats,
        const std::vector<float>& tau_f_ratio_list,
        const std::vector<double>& tau_s_list) {
    if (tau_f_ratio_list.empty()) {
        std::cerr << "CocoSketch eval skipped: tau_f_ratio_list is empty."
                  << std::endl;
        return;
    }
    if (tau_s_list.empty()) {
        std::cerr << "CocoSketch eval skipped: tau_s_list is empty."
                  << std::endl;
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

    const std::vector<FlowEstimate> estimates =
            collectEstimatesForEvaluation(exact_stats);
    std::unordered_map<uint32_t, uint32_t> estimated_frequency;
    std::unordered_map<uint32_t, double> estimated_dss;
    estimated_frequency.reserve(estimates.size());
    estimated_dss.reserve(estimates.size());
    for (const FlowEstimate& estimate : estimates) {
        estimated_frequency[estimate.flow_id] = estimate.frequency;
        estimated_dss[estimate.flow_id] = estimate.dss;
    }

    for (float tau_f_ratio : tau_f_ratio_list) {
        for (double tau_s : tau_s_list) {
            const TruthSnapshot truth_snapshot = buildTruthSnapshot(exact_stats,
                                                                    tau_f_ratio,
                                                                    reward,
                                                                    penalty,
                                                                    tau_s);

            uint32_t tp = 0;
            uint32_t fp = 0;
            uint32_t predicted = 0;
            for (const FlowEstimate& estimate : estimates) {
                const bool detected = estimate.frequency >= truth_snapshot.tau_f &&
                                      estimate.dss >= tau_s;
                if (!detected) {
                    continue;
                }
                ++predicted;
                if (truth_snapshot.contains(estimate.flow_id)) {
                    ++tp;
                } else {
                    ++fp;
                }
            }

            const uint32_t truth_count = truth_snapshot.truth_count;
            const uint32_t fn = truth_count - tp;
            double frequency_are_sum = 0.0;
            double dss_are_sum = 0.0;
            for (uint32_t flow_id : truth_snapshot.diffuse_flows) {
                const uint32_t true_f = exact_stats.trueFrequency(flow_id);
                const uint32_t true_c = exact_stats.trueCardinality(flow_id);
                const auto frequency_iter = estimated_frequency.find(flow_id);
                const uint32_t est_f = frequency_iter == estimated_frequency.end()
                                       ? 0
                                       : frequency_iter->second;
                const double frequency_error =
                        std::abs(static_cast<double>(est_f) -
                                 static_cast<double>(true_f));
                frequency_are_sum += frequency_error / static_cast<double>(true_f);

                const double true_dss =
                        computeDiffuseSupport(true_f, true_c, reward, penalty);
                const auto dss_iter = estimated_dss.find(flow_id);
                const double est_dss = dss_iter == estimated_dss.end()
                                       ? 0.0
                                       : dss_iter->second;
                const double dss_error = std::abs(est_dss - true_dss);
                dss_are_sum += true_dss == 0.0
                               ? 0.0
                               : dss_error / std::abs(true_dss);
            }

            const double frequency_are = truth_count == 0
                                         ? 0.0
                                         : frequency_are_sum /
                                           static_cast<double>(truth_count);
            const double dss_are = truth_count == 0
                                   ? 0.0
                                   : dss_are_sum /
                                     static_cast<double>(truth_count);
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

            std::cout << std::fixed << std::setprecision(6)
                      << " - CocoSketch: memory_kb=" << memory_kb
                      << ", depth=" << depth
                      << ", width=" << width
                      << ", retained_flows=" << estimates.size()
                      << ", tau_f_ratio=" << tau_f_ratio
                      << ", tau_f=" << truth_snapshot.tau_f
                      << ", tau_S=" << tau_s
                      << ", truth=" << truth_count
                      << ", predicted=" << predicted
                      << ", F1=" << f1
                      << ", precision=" << precision
                      << ", recall=" << recall
                      << ", dss_ARE=" << dss_are
                      << ", freq_ARE=" << frequency_are
                      << ", throughput_Mips=" << update_mips
                      << std::endl;
        }
    }
}
