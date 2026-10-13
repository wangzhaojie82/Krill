#include "header/DUET.h"

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

// Logical storage: a 64-bit full key and a 32-bit element counter.
constexpr uint32_t kBucketBits = 64 + 32;
constexpr uint32_t kFilterRowSeed = 799957137;
constexpr uint32_t kStableRowSeed = 17157137;

uint32_t hash32(uint32_t value, uint32_t seed) {
    uint32_t result = 0;
    MurmurHash3_x86_32(&value, sizeof(value), seed, &result);
    return result;
}

} // namespace

DUET::DUET(float input_memory_kb)
        : memory_kb(input_memory_kb),
          cu_sketch(input_memory_kb * cu_ratio) {
    const double total_bits =
            std::max(1.0, static_cast<double>(memory_kb) * 1024.0 * 8.0);
    filter_width = static_cast<uint32_t>(std::max(1.0, std::floor(
            total_bits * filter_ratio / (filter_depth * kBucketBits))));
    stable_width = static_cast<uint32_t>(std::max(1.0, std::floor(
            total_bits * stable_ratio / (stable_rows * kBucketBits))));
    filter.assign(filter_depth, std::vector<Bucket>(filter_width));
    stable.assign(stable_rows, std::vector<Bucket>(stable_width));
    filter_seeds.reserve(filter_depth);
    for (uint32_t row = 0; row < filter_depth; ++row) {
        filter_seeds.push_back(0x9e3779b9u + row * 0x85ebca6bu);
    }
}

uint64_t DUET::combineKey(uint32_t flow_id, uint32_t element_id) {
    return (static_cast<uint64_t>(flow_id) << 32u) |
           static_cast<uint64_t>(element_id);
}

std::pair<uint32_t, uint32_t> DUET::splitKey(uint64_t key) {
    return {static_cast<uint32_t>(key >> 32u), static_cast<uint32_t>(key)};
}

void DUET::insertIntoFilter(uint32_t flow_id, uint32_t element_id) {
    const uint32_t row = hash32(element_id, kFilterRowSeed) % filter_depth;
    const uint32_t column = hash32(flow_id, filter_seeds[row]) % filter_width;
    const uint64_t key = combineKey(flow_id, element_id);
    Bucket& bucket = filter[row][column];

    if (bucket.element == 0) {
        bucket = {key, 1};
    } else if (bucket.element == key) {
        if (bucket.count != std::numeric_limits<uint32_t>::max()) {
            ++bucket.count;
        }
    } else {
        --bucket.count;
        if (bucket.count == 0) {
            bucket = {key, 1};
        }
    }
}

void DUET::insertIntoStable(uint32_t flow_id, uint32_t element_id, uint32_t count) {
    if (count == 0) {
        return;
    }
    const uint64_t key = combineKey(flow_id, element_id);
    const uint32_t row = hash32(flow_id, kStableRowSeed) % stable_rows;
    int empty_column = -1;
    int minimum_column = -1;
    uint32_t minimum_count = std::numeric_limits<uint32_t>::max();

    for (uint32_t column = 0; column < stable_width; ++column) {
        Bucket& bucket = stable[row][column];
        if (bucket.element == key) {
            bucket.count = static_cast<uint32_t>(std::min<uint64_t>(
                    static_cast<uint64_t>(bucket.count) + count,
                    std::numeric_limits<uint32_t>::max()));
            return;
        }
        if (bucket.element == 0 && empty_column == -1) {
            empty_column = static_cast<int>(column);
        }
        if (bucket.element != 0 &&
            (minimum_column == -1 || bucket.count < minimum_count)) {
            minimum_count = bucket.count;
            minimum_column = static_cast<int>(column);
        }
    }

    if (empty_column != -1) {
        stable[row][empty_column] = {key, count};
        return;
    }

    Bucket& minimum_bucket = stable[row][minimum_column];
    if (minimum_bucket.count > count) {
        minimum_bucket.count -= count;
    } else if (minimum_bucket.count < count) {
        minimum_bucket = {key, count - minimum_bucket.count};
    } else {
        // Equal counts cancel completely; a zero-count key earns no DSS reward.
        minimum_bucket = {};
    }
}

void DUET::update(uint32_t flow_id, uint32_t element_id) {
    if (flow_id == 0 || element_id == 0) {
        return;
    }

    const uint32_t previous_frequency = cu_sketch.query(flow_id);
    cu_sketch.update(flow_id);
    if (previous_frequency < heavy_flow_threshold) {
        insertIntoFilter(flow_id, element_id);
        if (previous_frequency + 1 == heavy_flow_threshold) {
            for (uint32_t row = 0; row < filter_depth; ++row) {
                const uint32_t column = hash32(flow_id, filter_seeds[row]) % filter_width;
                Bucket& bucket = filter[row][column];
                const auto [stored_flow, stored_element] = splitKey(bucket.element);
                if (stored_flow == flow_id) {
                    insertIntoStable(flow_id, stored_element, bucket.count);
                    bucket = {};
                }
            }
        }
    } else {
        insertIntoStable(flow_id, element_id, 1);
    }
}

std::vector<FlowEstimate> DUET::collectEstimatesForEvaluation(
        const ExactStats& exact_stats) const {
    std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint64_t>> flow_elements;
    flow_elements.reserve(static_cast<size_t>(stable_rows) * stable_width);

    // As in the original DUET query, candidate full keys come from STable.
    for (const auto& row : stable) {
        for (const Bucket& bucket : row) {
            if (bucket.element == 0 || bucket.count == 0) {
                continue;
            }
            const auto [flow_id, element_id] = splitKey(bucket.element);
            if (flow_id == 0 || element_id == 0 ||
                exact_stats.trueFrequency(flow_id) == 0) {
                continue;
            }
            flow_elements[flow_id][element_id] += bucket.count;
        }
    }

    std::vector<FlowEstimate> estimates;
    estimates.reserve(flow_elements.size());
    for (const auto& flow_entry : flow_elements) {
        double dss = 0.0;
        for (const auto& element_entry : flow_entry.second) {
            dss += static_cast<double>(reward) - static_cast<double>(penalty) *
                   (static_cast<double>(element_entry.second) - 1.0);
        }
        // CU estimates flow frequency; only retained element counts determine DSS.
        estimates.push_back({flow_entry.first, cu_sketch.query(flow_entry.first), dss});
    }
    return estimates;
}

void DUET::eval(
        const std::vector<std::pair<uint32_t, uint32_t>>& data,
        const ExactStats& exact_stats,
        const std::vector<float>& tau_f_ratio_list,
        const std::vector<double>& tau_s_list) {
    if (tau_f_ratio_list.empty()) {
        std::cerr << "DUET eval skipped: tau_f_ratio_list is empty." << std::endl;
        return;
    }
    if (tau_s_list.empty()) {
        std::cerr << "DUET eval skipped: tau_s_list is empty." << std::endl;
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

    const std::vector<FlowEstimate> estimates = collectEstimatesForEvaluation(exact_stats);
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
                if (estimate.frequency < truth_snapshot.tau_f || estimate.dss < tau_s) {
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
            double frequency_are_sum = 0.0;
            double dss_are_sum = 0.0;
            for (uint32_t flow_id : truth_snapshot.diffuse_flows) {
                const uint32_t true_f = exact_stats.trueFrequency(flow_id);
                const uint32_t true_c = exact_stats.trueCardinality(flow_id);
                const auto frequency_iter = estimated_frequency.find(flow_id);
                const uint32_t est_f = frequency_iter == estimated_frequency.end()
                                       ? 0
                                       : frequency_iter->second;
                const double frequency_error = std::abs(static_cast<double>(est_f) -
                                                         static_cast<double>(true_f));
                frequency_are_sum += frequency_error / static_cast<double>(true_f);

                const double true_dss = computeDiffuseSupport(true_f, true_c, reward, penalty);
                const auto dss_iter = estimated_dss.find(flow_id);
                const double est_dss = dss_iter == estimated_dss.end() ? 0.0 : dss_iter->second;
                const double dss_error = std::abs(est_dss - true_dss);
                dss_are_sum += true_dss == 0.0 ? 0.0 : dss_error / std::abs(true_dss);
            }

            const double frequency_are = truth_count == 0 ? 0.0 :
                    frequency_are_sum / static_cast<double>(truth_count);
            const double dss_are = truth_count == 0 ? 0.0 :
                    dss_are_sum / static_cast<double>(truth_count);
            const double precision = tp + fp == 0 ? 0.0 :
                    static_cast<double>(tp) / static_cast<double>(tp + fp);
            const double recall = truth_count == 0 ? 0.0 :
                    static_cast<double>(tp) / static_cast<double>(truth_count);
            const double f1 = precision + recall == 0.0 ? 0.0 :
                    2.0 * precision * recall / (precision + recall);

            std::cout << std::fixed << std::setprecision(6)
                      << " - DUET: memory_kb=" << memory_kb
                      << ", cu_ratio=" << cu_ratio
                      << ", filter_ratio=" << filter_ratio
                      << ", stable_ratio=" << stable_ratio
                      << ", Nth=" << heavy_flow_threshold
                      << ", filter_depth=" << filter_depth
                      << ", filter_width=" << filter_width
                      << ", stable_rows=" << stable_rows
                      << ", stable_width=" << stable_width
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
