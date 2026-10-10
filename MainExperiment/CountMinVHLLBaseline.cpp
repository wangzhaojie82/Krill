#include "header/CountMinVHLLBaseline.h"

#include "header/ExactStats.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <unordered_set>

CountMinVHLLBaseline::CountMinVHLLBaseline(float input_memory_kb)
        : memory_kb(input_memory_kb),
          count_min(std::max(1.0f, input_memory_kb * count_min_ratio)),
          cardinality_estimator(std::max(1.0f, input_memory_kb * vhll_ratio)) {
    const double heap_bits =
            std::max(1.0, static_cast<double>(memory_kb) *
                          static_cast<double>(heap_ratio) * 1024.0 * 8.0);
    heap_capacity = static_cast<size_t>(
            std::max(1.0, std::floor(heap_bits / static_cast<double>(heap_entry_bits))));
    heap.reserve(heap_capacity);
    heap_index.reserve(heap_capacity);
}

bool CountMinVHLLBaseline::candidateLess(size_t lhs, size_t rhs) const {
    if (heap[lhs].est_f != heap[rhs].est_f) {
        return heap[lhs].est_f < heap[rhs].est_f;
    }
    return heap[lhs].flow_id < heap[rhs].flow_id;
}

void CountMinVHLLBaseline::swapCandidates(size_t lhs, size_t rhs) {
    std::swap(heap[lhs], heap[rhs]);
    heap_index[heap[lhs].flow_id] = lhs;
    heap_index[heap[rhs].flow_id] = rhs;
}

void CountMinVHLLBaseline::siftUp(size_t index) {
    while (index > 0) {
        const size_t parent = (index - 1) / 2;
        if (!candidateLess(index, parent)) {
            break;
        }
        swapCandidates(index, parent);
        index = parent;
    }
}

void CountMinVHLLBaseline::siftDown(size_t index) {
    while (true) {
        const size_t left = index * 2 + 1;
        const size_t right = left + 1;
        size_t smallest = index;
        if (left < heap.size() && candidateLess(left, smallest)) {
            smallest = left;
        }
        if (right < heap.size() && candidateLess(right, smallest)) {
            smallest = right;
        }
        if (smallest == index) {
            break;
        }
        swapCandidates(index, smallest);
        index = smallest;
    }
}

void CountMinVHLLBaseline::updateCandidate(uint32_t flow_id, uint32_t est_f) {
    const auto iter = heap_index.find(flow_id);
    if (iter != heap_index.end()) {
        const size_t index = iter->second;
        heap[index].est_f = est_f;
        siftDown(index);
        return;
    }

    if (heap.size() < heap_capacity) {
        heap.push_back({flow_id, est_f});
        const size_t index = heap.size() - 1;
        heap_index[flow_id] = index;
        siftUp(index);
        return;
    }

    if (!heap.empty() && est_f > heap.front().est_f) {
        heap_index.erase(heap.front().flow_id);
        heap.front() = {flow_id, est_f};
        heap_index[flow_id] = 0;
        siftDown(0);
    }
}

void CountMinVHLLBaseline::update(uint32_t flow_id, uint32_t element_id) {
    count_min.update(flow_id);
    cardinality_estimator.update(flow_id, element_id);
    updateCandidate(flow_id, count_min.query(flow_id));
}

std::vector<FlowEstimate> CountMinVHLLBaseline::collectEstimatesForEvaluation(
        const ExactStats& exact_stats) {
    cardinality_estimator.prepareQuery();

    std::vector<FlowEstimate> estimates;
    estimates.reserve(heap.size());
    for (const Candidate& candidate : heap) {
        if (candidate.flow_id == 0 || exact_stats.trueFrequency(candidate.flow_id) == 0) {
            continue;
        }
        const uint32_t est_f = count_min.query(candidate.flow_id);
        const uint32_t est_c = cardinality_estimator.query(candidate.flow_id);
        estimates.push_back({candidate.flow_id,
                             est_f,
                             computeDiffuseSupport(est_f, est_c, reward, penalty)});
    }
    return estimates;
}

void CountMinVHLLBaseline::eval(
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

    cardinality_estimator.prepareQuery();

    struct CandidateStat {
        uint32_t flow_id = 0;
        uint32_t est_f = 0;
        uint32_t est_c = 0;
        double est_s = 0.0;
    };

    std::vector<CandidateStat> candidate_stats;
    candidate_stats.reserve(heap.size());
    std::unordered_map<uint32_t, uint32_t> estimated_frequency;
    estimated_frequency.reserve(heap.size());
    std::unordered_map<uint32_t, double> estimated_dss;
    estimated_dss.reserve(heap.size());

    for (const auto& candidate : heap) {
        const uint32_t est_f = count_min.query(candidate.flow_id);
        const uint32_t est_c = cardinality_estimator.query(candidate.flow_id);
        const double est_s = computeDiffuseSupport(est_f, est_c, reward, penalty);
        candidate_stats.push_back({candidate.flow_id, est_f, est_c, est_s});
        estimated_frequency[candidate.flow_id] = est_f;
        estimated_dss[candidate.flow_id] = est_s;
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
            for (const auto& candidate : candidate_stats) {
                const bool detected = candidate.est_f >= truth_snapshot.tau_f &&
                                      candidate.est_s >= tau_s;
                if (!detected) {
                    continue;
                }
                ++predicted;
                if (truth_snapshot.contains(candidate.flow_id)) {
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
                      << " - CountMinVHLLBaseline: memory_kb=" << memory_kb
                      << ", cm_ratio=" << count_min_ratio
                      << ", heap_ratio=" << heap_ratio
                      << ", vhll_ratio=" << vhll_ratio
                      << ", reward=" << reward
                      << ", penalty=" << penalty
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


