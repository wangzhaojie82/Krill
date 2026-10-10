#include "header/BloomFilterBaseline.h"

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

namespace {

double normalizedHash01(uint32_t hash_value) {
    return static_cast<double>(hash_value) /
           static_cast<double>(std::numeric_limits<uint32_t>::max());
}

} // namespace

BloomFilterBaseline::BloomFilter::BloomFilter(uint64_t bit_num, uint32_t hash_num)
        : bit_num(std::max<uint64_t>(64, bit_num)),
          hash_num(std::max<uint32_t>(1, hash_num)) {
    words.assign((this->bit_num + 63u) / 64u, 0);
}

bool BloomFilterBaseline::BloomFilter::setBit(uint64_t index) {
    uint64_t& word = words[index >> 6u];
    const uint64_t mask = 1ull << (index & 63u);
    if ((word & mask) != 0) {
        return false;
    }
    word |= mask;
    ++set_bit_num;
    return true;
}

bool BloomFilterBaseline::BloomFilter::queryAndInsert(uint32_t flow_id,
                                                        uint32_t element_id) {
    const uint32_t key[2] = {flow_id, element_id};
    uint32_t h1 = 0;
    uint32_t h2 = 0;
    MurmurHash3_x86_32(key, static_cast<int>(sizeof(key)), seed1, &h1);
    MurmurHash3_x86_32(key, static_cast<int>(sizeof(key)), seed2, &h2);
    if (h2 == 0) {
        h2 = 0x9e3779b9u;
    }

    bool seen = true;
    uint64_t index_value = static_cast<uint64_t>(h1);
    const uint64_t step = static_cast<uint64_t>(h2);
    for (uint32_t i = 0; i < hash_num; ++i) {
        const uint64_t index = index_value % bit_num;
        if (setBit(index)) {
            seen = false;
        }
        index_value += step + static_cast<uint64_t>(2u * i + 1u);
    }
    return seen;
}

double BloomFilterBaseline::BloomFilter::fillRatio() const {
    return bit_num == 0
           ? 0.0
           : static_cast<double>(set_bit_num) / static_cast<double>(bit_num);
}

double BloomFilterBaseline::BloomFilter::estimatedFalsePositiveRate() const {
    if (bit_num == 0 || hash_num == 0) {
        return 0.0;
    }
    return std::pow(fillRatio(), static_cast<double>(hash_num));
}

BloomFilterBaseline::BloomFilterBaseline(float input_memory_kb)
        : memory_kb(input_memory_kb),
          pair_seen(64, bloom_hash_num) {
    if (cell_num == 0) {
        throw std::invalid_argument("cell_num must be positive.");
    }
    if (bloom_ratio <= 0.0f || bloom_ratio >= 1.0f) {
        throw std::invalid_argument("bloom_ratio must be in (0, 1).");
    }

    const uint64_t total_bits = static_cast<uint64_t>(
            std::max<double>(1.0,
                             static_cast<double>(memory_kb) * 1024.0 * 8.0));
    uint64_t bloom_bits = static_cast<uint64_t>(
            static_cast<double>(total_bits) * static_cast<double>(bloom_ratio));
    bloom_bits = std::max<uint64_t>(64, bloom_bits);
    if (bloom_bits >= total_bits) {
        bloom_bits = total_bits / 2;
    }

    const uint64_t sketch_bits = total_bits - bloom_bits;
    const uint64_t cell_bits = static_cast<uint64_t>(flow_id_bits) +
                               static_cast<uint64_t>(dss_bits) +
                               static_cast<uint64_t>(freq_bits);
    const uint64_t bucket_bits =
            std::max<uint64_t>(1, static_cast<uint64_t>(cell_num) * cell_bits);
    bucket_num = static_cast<uint32_t>(std::max<uint64_t>(1, sketch_bits / bucket_bits));

    pair_seen = BloomFilter(bloom_bits, bloom_hash_num);
    cells.assign(static_cast<size_t>(bucket_num) * static_cast<size_t>(cell_num),
                 Cell{});
}

int64_t BloomFilterBaseline::exponentialProbabilisticDecay(
        uint32_t flow_id,
        uint32_t element_id,
        const Cell& min_cell,
        int64_t reward_value) {
    constexpr double decay_base = 1.08;
    const double val =
            std::floor(std::log2(static_cast<double>(std::max<int64_t>(0, min_cell.dss)) + 1.0));
    const double probability = std::pow(decay_base, -val);
    const uint32_t key[3] = {flow_id, element_id, min_cell.flow_id};
    uint32_t hash_value = 0;
    MurmurHash3_x86_32(key, static_cast<int>(sizeof(key)), 6001, &hash_value);
    return normalizedHash01(hash_value) <= probability ? reward_value : 0;
}

void BloomFilterBaseline::update(uint32_t flow_id, uint32_t element_id) {
    constexpr uint32_t bucket_seed = 1001;
    const int64_t reward_value = static_cast<int64_t>(reward);
    const int64_t penalty_value = static_cast<int64_t>(penalty);
    const bool seen = pair_seen.queryAndInsert(flow_id, element_id);
    const int64_t delta = seen ? -penalty_value : reward_value;

    uint32_t bucket_hash = 0;
    MurmurHash3_x86_32(&flow_id, sizeof(flow_id), bucket_seed, &bucket_hash);
    Cell* bucket = cells.data() +
                   static_cast<size_t>(bucket_hash % bucket_num) *
                   static_cast<size_t>(cell_num);

    Cell* empty_cell = nullptr;
    Cell* min_cell = nullptr;
    for (uint32_t i = 0; i < cell_num; ++i) {
        Cell& cell = bucket[i];
        if (cell.flow_id == flow_id) {
            cell.dss += delta;
            ++cell.freq;
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
        empty_cell->flow_id = flow_id;
        empty_cell->dss = delta;
        empty_cell->freq = 1;
        return;
    }

    if (delta < 0 || min_cell == nullptr) {
        return;
    }

    min_cell->dss -= exponentialProbabilisticDecay(flow_id,
                                                   element_id,
                                                   *min_cell,
                                                   reward_value);
    if (min_cell->dss <= 0) {
        min_cell->flow_id = flow_id;
        min_cell->dss = reward_value;
        min_cell->freq = 1;
    }
}

std::vector<FlowEstimate> BloomFilterBaseline::collectEstimatesForEvaluation(
        const ExactStats& exact_stats) {
    std::vector<FlowEstimate> estimates;
    estimates.reserve(cells.size());
    for (const Cell& cell : cells) {
        if (cell.flow_id == 0 || exact_stats.trueFrequency(cell.flow_id) == 0) {
            continue;
        }
        estimates.push_back({cell.flow_id,
                             cell.freq,
                             static_cast<double>(cell.dss)});
    }
    return estimates;
}

void BloomFilterBaseline::eval(
        const std::vector<std::pair<uint32_t, uint32_t>>& data,
        const ExactStats& exact_stats,
        const std::vector<float>& tau_f_ratio_list,
        const std::vector<double>& tau_s_list) {
    if (tau_f_ratio_list.empty()) {
        std::cerr << "eval skipped: tau_f_ratio_list is empty." << std::endl;
        return;
    }
    if (tau_s_list.empty()) {
        std::cerr << "eval skipped: tau_s_list is empty." << std::endl;
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
        int64_t est_s = 0;
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
        estimated_dss[cell.flow_id] = static_cast<double>(cell.dss);
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
            for (const auto& cell : resident_stats) {
                const bool truth = truth_snapshot.contains(cell.flow_id);
                const bool detected = cell.est_f >= tau_f &&
                                      static_cast<double>(cell.est_s) >= tau_s;
                if (detected) {
                    ++predicted;
                    if (truth) {
                        ++tp;
                    } else {
                        ++fp;
                    }
                }
            }

            const uint32_t truth_count = truth_snapshot.truth_count;
            const uint32_t fn = truth_count - tp;
            double frequency_are_sum = 0.0;
            double dss_are_sum = 0.0;
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
                const double est_dss =
                        dss_iter == estimated_dss.end() ? 0.0 : dss_iter->second;
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
                      << " - BloomFilterBaseline: memory_kb=" << memory_kb
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


