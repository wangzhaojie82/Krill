#include "header/FLASKSketch.h"

#include "header/ExactStats.h"
#include "header/MurmurHash3.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace {

double hllAlpha(uint32_t register_num) {
    if (register_num == 16) {
        return 0.673;
    }
    if (register_num == 32) {
        return 0.697;
    }
    if (register_num == 64) {
        return 0.709;
    }
    return 0.7213 / (1.0 + 1.079 / static_cast<double>(register_num));
}

uint32_t medianValue(std::vector<uint32_t>& values) {
    if (values.empty()) {
        return 0;
    }
    const size_t middle = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + middle, values.end());
    return values[middle];
}

} // namespace

FLASKSketch::FLASKSketch(float input_memory_kb)
        : memory_kb(input_memory_kb) {
    if ((hll_register_num & (hll_register_num - 1u)) != 0u) {
        throw std::invalid_argument("hll_register_num must be a power of two.");
    }
    hll_index_bits = 0;
    uint32_t temp_register_num = hll_register_num;
    while (temp_register_num > 1) {
        ++hll_index_bits;
        temp_register_num >>= 1u;
    }
    hll_alpha = hllAlpha(hll_register_num);

    const uint32_t bucket_meta_bits = 32 + 32 + 32 + 8 + 1;
    const uint32_t hll_bits = hll_register_num * 5;
    const uint32_t bucket_bits = bucket_meta_bits + hll_bits;
    const double total_bits =
            std::max(1.0, static_cast<double>(memory_kb) * 1024.0 * 8.0);
    bucket_num = static_cast<uint32_t>(
            std::max(1.0, std::floor(total_bits /
                                      static_cast<double>(row_num * bucket_bits))));
    buckets.resize(static_cast<size_t>(row_num) * bucket_num);
    for (Bucket& bucket : buckets) {
        bucket.hll_sum = static_cast<double>(hll_register_num);
        bucket.hll_zero_num = hll_register_num;
    }
    hll_registers.assign(buckets.size() * hll_register_num, 0);
}

size_t FLASKSketch::bucketOffset(uint32_t row, uint32_t bucket) const {
    return static_cast<size_t>(row) * bucket_num + bucket;
}

size_t FLASKSketch::registerOffset(size_t bucket_offset, uint32_t register_index) const {
    return bucket_offset * hll_register_num + register_index;
}

uint32_t FLASKSketch::hash32(uint32_t value, uint32_t seed) const {
    uint32_t hash_value = 0;
    MurmurHash3_x86_32(&value, sizeof(value), seed, &hash_value);
    return hash_value;
}

uint32_t FLASKSketch::hashPair(uint32_t first, uint32_t second, uint32_t seed) const {
    const uint64_t key = (static_cast<uint64_t>(first) << 32u) | second;
    uint32_t hash_value = 0;
    MurmurHash3_x86_32(&key, sizeof(key), seed, &hash_value);
    return hash_value;
}

uint8_t FLASKSketch::leadingZeros(uint32_t hash_value) const {
    uint8_t count = 1;
    while ((hash_value & 0x80000000u) == 0 && count < 31) {
        ++count;
        hash_value <<= 1u;
    }
    return count;
}

uint32_t FLASKSketch::bucketIndex(uint32_t flow_id, uint32_t row) const {
    return hash32(flow_id, bucket_seed + row * 131u) % bucket_num;
}

void FLASKSketch::resetHLL(size_t bucket_offset) {
    const size_t begin = bucket_offset * hll_register_num;
    std::fill(hll_registers.begin() + static_cast<std::ptrdiff_t>(begin),
              hll_registers.begin() + static_cast<std::ptrdiff_t>(begin + hll_register_num),
              0);
    Bucket& bucket = buckets[bucket_offset];
    bucket.hll_sum = static_cast<double>(hll_register_num);
    bucket.hll_zero_num = hll_register_num;
}

void FLASKSketch::updateHLL(size_t bucket_offset, uint32_t register_index, uint8_t rank) {
    uint8_t& reg = hll_registers[registerOffset(bucket_offset, register_index)];
    if (rank > reg) {
        Bucket& bucket = buckets[bucket_offset];
        bucket.hll_sum -= std::ldexp(1.0, -static_cast<int>(reg));
        bucket.hll_sum += std::ldexp(1.0, -static_cast<int>(rank));
        if (reg == 0) {
            --bucket.hll_zero_num;
        }
        reg = rank;
    }
}

uint32_t FLASKSketch::queryHLL(size_t bucket_offset) const {
    const Bucket& bucket = buckets[bucket_offset];
    double estimate = hll_alpha *
                      static_cast<double>(hll_register_num) *
                      static_cast<double>(hll_register_num) / bucket.hll_sum;
    if (estimate <= 2.5 * static_cast<double>(hll_register_num) &&
        bucket.hll_zero_num > 0) {
        estimate = static_cast<double>(hll_register_num) *
                   std::log(static_cast<double>(hll_register_num) /
                            static_cast<double>(bucket.hll_zero_num));
    }
    return static_cast<uint32_t>(std::max(0.0, std::round(estimate)));
}

void FLASKSketch::reinitializeBucket(size_t bucket_offset,
                                     uint32_t flow_id,
                                     uint32_t hll_register_index,
                                     uint8_t hll_rank,
                                     uint8_t leading_zero) {
    Bucket& bucket = buckets[bucket_offset];
    bucket.key = flow_id;
    bucket.frequency = 1;
    bucket.auxiliary = 0;
    bucket.leading_zero = leading_zero;
    bucket.occupied = true;
    resetHLL(bucket_offset);
    updateHLL(bucket_offset, hll_register_index, hll_rank);
}

void FLASKSketch::update(uint32_t flow_id, uint32_t element_id) {
    const uint8_t pair_lz = leadingZeros(hashPair(flow_id, element_id, pair_seed));
    const uint32_t hll_hash = hash32(element_id, hll_seed);
    const uint32_t hll_register_index = hll_hash >> (32u - hll_index_bits);
    uint32_t hll_remain_hash = hll_hash << hll_index_bits;
    uint8_t hll_rank = 1;
    while ((hll_remain_hash & 0x80000000u) == 0 && hll_rank < 31) {
        ++hll_rank;
        hll_remain_hash <<= 1u;
    }

    for (uint32_t row = 0; row < row_num; ++row) {
        const uint32_t index = bucketIndex(flow_id, row);
        const size_t offset = bucketOffset(row, index);
        Bucket& bucket = buckets[offset];

        if (!bucket.occupied) {
            reinitializeBucket(offset, flow_id, hll_register_index, hll_rank, pair_lz);
            continue;
        }

        if (bucket.key == flow_id) {
            bucket.leading_zero = std::max(bucket.leading_zero, pair_lz);
            ++bucket.frequency;
            ++bucket.auxiliary;
            updateHLL(offset, hll_register_index, hll_rank);
            continue;
        }

        ++bucket.frequency;
        --bucket.auxiliary;
        const uint32_t spread = std::max(1u, queryHLL(offset));
        const double resident_frequency =
                std::max(0.0, (static_cast<double>(bucket.frequency) +
                               static_cast<double>(bucket.auxiliary)) / 2.0);
        const double repeat_ratio =
                resident_frequency <= 0.0
                ? std::numeric_limits<double>::infinity()
                : resident_frequency / static_cast<double>(spread);
        const double old_score = repeat_ratio <= 0.0
                                 ? std::numeric_limits<double>::infinity()
                                 : alpha * static_cast<double>(bucket.leading_zero) / repeat_ratio;
        if (static_cast<double>(pair_lz) > old_score) {
            bucket.key = flow_id;
            bucket.leading_zero = pair_lz;
        }
    }
}

std::pair<uint32_t, uint32_t> FLASKSketch::query(uint32_t flow_id) const {
    std::vector<uint32_t> frequency_estimates;
    frequency_estimates.reserve(row_num);
    uint32_t cardinality_estimate = std::numeric_limits<uint32_t>::max();

    for (uint32_t row = 0; row < row_num; ++row) {
        const uint32_t index = bucketIndex(flow_id, row);
        const size_t offset = bucketOffset(row, index);
        const Bucket& bucket = buckets[offset];
        cardinality_estimate = std::min(cardinality_estimate, queryHLL(offset));

        double estimate = 0.0;
        if (bucket.occupied && bucket.key == flow_id) {
            estimate = (static_cast<double>(bucket.frequency) +
                        static_cast<double>(bucket.auxiliary)) / 2.0;
        } else {
            estimate = (static_cast<double>(bucket.frequency) -
                        static_cast<double>(bucket.auxiliary)) / 2.0;
        }
        frequency_estimates.push_back(
                static_cast<uint32_t>(std::max(0.0, std::round(estimate))));
    }

    if (cardinality_estimate == std::numeric_limits<uint32_t>::max()) {
        cardinality_estimate = 0;
    }
    return {medianValue(frequency_estimates), cardinality_estimate};
}

std::vector<FlowEstimate> FLASKSketch::collectEstimatesForEvaluation(
        const ExactStats& exact_stats) const {
    std::unordered_set<uint32_t> retained_flow_ids;
    retained_flow_ids.reserve(buckets.size());
    for (const Bucket& bucket : buckets) {
        if (bucket.occupied && bucket.key != 0) {
            retained_flow_ids.insert(bucket.key);
        }
    }

    std::vector<FlowEstimate> estimates;
    estimates.reserve(retained_flow_ids.size());
    for (uint32_t flow_id : retained_flow_ids) {
        if (exact_stats.trueFrequency(flow_id) == 0) {
            continue;
        }
        const auto estimate = query(flow_id);
        estimates.push_back({flow_id,
                             estimate.first,
                             computeDiffuseSupport(estimate.first,
                                                   estimate.second,
                                                   reward,
                                                   penalty)});
    }
    return estimates;
}

void FLASKSketch::eval(
        const std::vector<std::pair<uint32_t, uint32_t>>& data,
        const ExactStats& exact_stats,
        const std::vector<float>& tau_f_ratio_list,
        const std::vector<double>& tau_s_list) {

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

    const std::vector<FlowEstimate> candidate_estimates =
            collectEstimatesForEvaluation(exact_stats);
    std::unordered_map<uint32_t, FlowEstimate> estimates;
    estimates.reserve(candidate_estimates.size());
    for (const FlowEstimate& estimate : candidate_estimates) {
        estimates[estimate.flow_id] = estimate;
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
            for (const auto& item : estimates) {
                const uint32_t flow_id = item.first;
                const uint32_t est_f = item.second.frequency;
                const double est_s = item.second.dss;
                const bool detected = est_f >= truth_snapshot.tau_f && est_s >= tau_s;
                if (!detected) {
                    continue;
                }
                ++predicted;
                if (truth_snapshot.contains(flow_id)) {
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
                const auto estimate_iter = estimates.find(flow_id);
                const uint32_t est_f = estimate_iter == estimates.end()
                                       ? 0
                                       : estimate_iter->second.frequency;
                const double abs_error = std::abs(static_cast<double>(est_f) -
                                                  static_cast<double>(true_f));
                frequency_are_sum += abs_error / static_cast<double>(true_f);

                const double true_dss = computeDiffuseSupport(true_f, true_c, reward, penalty);
                const double est_dss = estimate_iter == estimates.end()
                                       ? 0.0
                                       : estimate_iter->second.dss;
                const double dss_abs_error = std::abs(est_dss - true_dss);
                dss_are_sum += true_dss == 0.0 ? 0.0 : dss_abs_error / std::abs(true_dss);
            }
            const double frequency_are = truth_count == 0
                                         ? 0.0
                                         : frequency_are_sum / static_cast<double>(truth_count);
            const double dss_are = truth_count == 0
                                   ? 0.0
                                   : dss_are_sum / static_cast<double>(truth_count);
            const float precision = tp + fp == 0
                                    ? 0.0f
                                    : static_cast<float>(tp) / static_cast<float>(tp + fp);
            const float recall = tp + fn == 0
                                 ? 0.0f
                                 : static_cast<float>(tp) / static_cast<float>(tp + fn);
            const float f1 = precision + recall == 0.0f
                             ? 0.0f
                             : 2.0f * precision * recall / (precision + recall);

            std::cout << std::fixed << std::setprecision(6)
                      << " - FLASK: memory_kb=" << memory_kb
                      << ", hll_regs=" << hll_register_num
                      << ", alpha=" << alpha
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


