#include "header/ErrorBreakdownExperiment.h"

#include "header/BloomFilterBaseline.h"
#include "header/CountMinVHLLBaseline.h"
#include "header/DatasetLoader.h"
#include "header/Krill.h"
#include "header/ExactStats.h"
#include "header/FLASKSketch.h"
#include "header/CocoSketch.h"
#include "header/DUET.h"
#include "header/FlowEstimate.h"
#include "header/PISketch.h"

#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

constexpr uint32_t kReward = 12;
constexpr uint32_t kPenalty = 1;

struct DetectionErrorBreakdown {
    std::string method;
    float memory_kb = 0.0f;
    uint32_t tau_f = 0;
    double tau_s = 0.0;
    uint64_t retained = 0;
    uint64_t negative_flows = 0;
    uint64_t truth = 0;
    uint64_t predicted = 0;
    uint64_t tp = 0;
    uint64_t fp = 0;
    uint64_t fn = 0;
    uint64_t frequency_fp = 0;
    uint64_t dss_fp = 0;
    uint64_t not_retained = 0;
    uint64_t frequency_miss = 0;
    uint64_t dss_miss = 0;
    double precision = 0.0;
    double recall = 0.0;
    double f1 = 0.0;
};

DetectionErrorBreakdown computeBreakdown(
        const std::string& method,
        float memory_kb,
        const std::vector<FlowEstimate>& estimates,
        const ExactStats& exact_stats,
        float tau_f_ratio,
        double tau_s) {
    const TruthSnapshot truth_snapshot = buildTruthSnapshot(exact_stats,
                                                            tau_f_ratio,
                                                            kReward,
                                                            kPenalty,
                                                            tau_s);

    std::unordered_map<uint32_t, FlowEstimate> estimate_by_flow;
    estimate_by_flow.reserve(estimates.size());
    for (const FlowEstimate& estimate : estimates) {
        if (estimate.flow_id != 0) {
            estimate_by_flow[estimate.flow_id] = estimate;
        }
    }

    DetectionErrorBreakdown result;
    result.method = method;
    result.memory_kb = memory_kb;
    result.tau_f = truth_snapshot.tau_f;
    result.tau_s = tau_s;
    result.retained = estimate_by_flow.size();
    result.truth = truth_snapshot.truth_count;
    const uint64_t total_flows =
            static_cast<uint64_t>(exact_stats.true_flow_frequency.size());
    result.negative_flows = total_flows >= result.truth
                            ? total_flows - result.truth
                            : 0;

    for (const auto& item : estimate_by_flow) {
        const FlowEstimate& estimate = item.second;
        if (estimate.frequency < truth_snapshot.tau_f || estimate.dss < tau_s) {
            continue;
        }

        ++result.predicted;
        if (truth_snapshot.contains(estimate.flow_id)) {
            ++result.tp;
            continue;
        }

        ++result.fp;
        const uint32_t true_f = exact_stats.trueFrequency(estimate.flow_id);
        if (true_f < truth_snapshot.tau_f) {
            ++result.frequency_fp;
        } else {
            ++result.dss_fp;
        }
    }

    for (uint32_t flow_id : truth_snapshot.diffuse_flows) {
        const auto estimate_iter = estimate_by_flow.find(flow_id);
        if (estimate_iter == estimate_by_flow.end()) {
            ++result.not_retained;
            continue;
        }

        const FlowEstimate& estimate = estimate_iter->second;
        if (estimate.frequency < truth_snapshot.tau_f) {
            ++result.frequency_miss;
        } else if (estimate.dss < tau_s) {
            ++result.dss_miss;
        }
    }

    result.fn = result.not_retained + result.frequency_miss + result.dss_miss;
    result.precision = result.predicted == 0
                       ? 0.0
                       : static_cast<double>(result.tp) /
                         static_cast<double>(result.predicted);
    result.recall = result.truth == 0
                    ? 0.0
                    : static_cast<double>(result.tp) /
                      static_cast<double>(result.truth);
    result.f1 = result.precision + result.recall == 0.0
                ? 0.0
                : 2.0 * result.precision * result.recall /
                  (result.precision + result.recall);
    return result;
}

void printBreakdown(const DetectionErrorBreakdown& result) {
    const auto printRate = [](const char* name,
                              uint64_t numerator,
                              uint64_t denominator) {
        std::cout << ", " << name << '=';
        if (denominator == 0) {
            std::cout << "N/A";
            return;
        }
        std::cout << static_cast<double>(numerator) /
                    static_cast<double>(denominator);
    };

    std::cout << std::fixed << std::setprecision(6)
              << "Error Breakdown: method=" << result.method
              << ", memory_kb=" << result.memory_kb
              << ", tau_f=" << result.tau_f
              << ", tau_S=" << result.tau_s
              << ", retained=" << result.retained
              << ", truth=" << result.truth
              << ", non_target_flows=" << result.negative_flows
              << ", predicted=" << result.predicted
              << ", precision=" << result.precision << ", recall=" << result.recall << ", F1="
              << result.f1
              << ", Frequency_FP=" << result.frequency_fp
              << ", DSS_FP=" << result.dss_fp
              << ", Not_Retained=" << result.not_retained
              << ", Frequency_Miss=" << result.frequency_miss
              << ", DSS_Miss=" << result.dss_miss
              ;

    printRate("FDR", result.fp, result.predicted);
    printRate("Frequency_FDR", result.frequency_fp, result.predicted);
    printRate("DSS_FDR", result.dss_fp, result.predicted);
    printRate("FPR", result.fp, result.negative_flows);
    printRate("Frequency_FPR", result.frequency_fp, result.negative_flows);
    printRate("DSS_FPR", result.dss_fp, result.negative_flows);
    printRate("FNR", result.fn, result.truth);
    printRate("Not_Retained_FNR", result.not_retained, result.truth);
    printRate("Frequency_Miss_FNR", result.frequency_miss, result.truth);
    printRate("DSS_Miss_FNR", result.dss_miss, result.truth);
    printRate("Frequency_FP_Share", result.frequency_fp, result.fp);
    printRate("DSS_FP_Share", result.dss_fp, result.fp);
    printRate("Not_Retained_FN_Share", result.not_retained, result.fn);
    printRate("Frequency_Miss_FN_Share", result.frequency_miss, result.fn);
    printRate("DSS_Miss_FN_Share", result.dss_miss, result.fn);
    std::cout << std::endl;
}

template<typename Sketch>
void updateSketch(Sketch& sketch, const PairData& data) {
    for (const auto& item : data) {
        sketch.update(item.first, item.second);
    }
}

} // namespace

void runDetectionErrorBreakdown(
        const LoadedDataset& loaded,
        const ErrorBreakdownExperimentConfig& config) {
    if (config.memory_kb_list.empty()) {
        std::cerr << "Error breakdown skipped: memory_kb_list is empty." << std::endl;
        return;
    }

    std::cout << "\n=== Detection Error Breakdown ===" << std::endl;
    for (float memory_kb : config.memory_kb_list) {
        {
            Krill sketch(memory_kb);
            updateSketch(sketch, loaded.data);
            printBreakdown(computeBreakdown("DS",
                                            memory_kb,
                                            sketch.collectEstimatesForEvaluation(
                                                    loaded.exact_stats),
                                            loaded.exact_stats,
                                            config.tau_f_ratio,
                                            config.tau_s));
        }
        {
            BloomFilterBaseline sketch(memory_kb);
            updateSketch(sketch, loaded.data);
            printBreakdown(computeBreakdown("BloomFilterBaseline",
                                            memory_kb,
                                            sketch.collectEstimatesForEvaluation(
                                                    loaded.exact_stats),
                                            loaded.exact_stats,
                                            config.tau_f_ratio,
                                            config.tau_s));
        }
        {
            PISketch sketch(memory_kb);
            updateSketch(sketch, loaded.data);
            printBreakdown(computeBreakdown("PIS",
                                            memory_kb,
                                            sketch.collectEstimatesForEvaluation(
                                                    loaded.exact_stats),
                                            loaded.exact_stats,
                                            config.tau_f_ratio,
                                            config.tau_s));
        }
        {
            CountMinVHLLBaseline sketch(memory_kb);
            updateSketch(sketch, loaded.data);
            printBreakdown(computeBreakdown("CM-vHLL",
                                            memory_kb,
                                            sketch.collectEstimatesForEvaluation(
                                                    loaded.exact_stats),
                                            loaded.exact_stats,
                                            config.tau_f_ratio,
                                            config.tau_s));
        }
        {
            FLASKSketch sketch(memory_kb);
            updateSketch(sketch, loaded.data);
            printBreakdown(computeBreakdown("FLASK",
                                            memory_kb,
                                            sketch.collectEstimatesForEvaluation(
                                                    loaded.exact_stats),
                                            loaded.exact_stats,
                                            config.tau_f_ratio,
                                            config.tau_s));
        }
        {
            CocoSketch sketch(memory_kb);
            updateSketch(sketch, loaded.data);
            printBreakdown(computeBreakdown("CocoSketch",
                                            memory_kb,
                                            sketch.collectEstimatesForEvaluation(
                                                    loaded.exact_stats),
                                            loaded.exact_stats,
                                            config.tau_f_ratio,
                                            config.tau_s));
        }
        {
            DUET sketch(memory_kb);
            updateSketch(sketch, loaded.data);
            printBreakdown(computeBreakdown("DUET",
                                            memory_kb,
                                            sketch.collectEstimatesForEvaluation(
                                                    loaded.exact_stats),
                                            loaded.exact_stats,
                                            config.tau_f_ratio,
                                            config.tau_s));
        }
    }
}


