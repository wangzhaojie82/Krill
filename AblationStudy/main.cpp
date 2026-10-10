#include "header/AblationSketch.h"
#include "header/CompensationAblation.h"
#include "header/DatasetLoader.h"
#include "header/ReplacementAblation.h"

#include <chrono>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace {

    void printResult(const AblationResult &result) {
        std::cout << std::fixed << std::setprecision(6)
                  << "variant=" << result.variant
                  << ", total_kb=" << result.total_memory_kb
                  << ", pair_kb=";
        if (result.pair_memory_limited) {
            std::cout << result.pair_memory_kb;
        } else {
            std::cout << "unlimited";
        }
        std::cout << ", residence_kb=";
        if (result.residence_memory_limited) {
            std::cout << result.residence_memory_kb;
        } else {
            std::cout << "unlimited";
        }
        std::cout
                << ", cells_per_bucket=" << result.cell_num
                << ", tau_f=" << result.tau_f
                << ", tau_S=" << result.tau_s
                << ", resident=" << result.resident_flows
                << ", truth=" << result.ground_truth_flows
                << ", predicted=" << result.predicted_flows
                << ", TP/FP/FN=" << result.tp << '/' << result.fp << '/' << result.fn
                << ", precision=" << result.precision
                << ", recall=" << result.recall
                << ", F1=" << result.f1
                << ", dss_ARE=" << result.dss_are
                << ", freq_ARE=" << result.frequency_are
                << ", Mips=" << result.update_mips
                << std::endl;
    }

    void runVariant(const AblationConfig &config,
                    const std::vector<std::pair<uint32_t, uint32_t>> &data,
                    const ExactStats &exact_stats,
                    const std::vector<float> &tau_f_ratio_list,
                    const std::vector<double> &tau_s_list,
                    const std::string &study_name) {
        std::cout << "\n[" << study_name << "] "
                  << AblationSketch::pairMechanismName(config.pair_mechanism)
                  << "+"
                  << AblationSketch::residenceMechanismName(config.residence_mechanism)
                  << ", total_kb=" << config.total_memory_kb
                  << ", pair_kb=";
        if (config.pair_mechanism == PairMechanism::Exact) {
            std::cout << "unlimited";
        } else {
            std::cout << config.total_memory_kb * config.pair_memory_ratio;
        }
        std::cout << ", residence_kb=";
        if (config.residence_mechanism == ResidenceMechanism::Ideal) {
            std::cout << "unlimited";
        } else {
            std::cout << config.total_memory_kb * config.residence_memory_ratio;
        }
        std::cout
                  << ", cells_per_bucket=" << config.cell_num << std::endl;
        const AblationSketch sketch(config);
        for (const AblationResult &result: sketch.run(data,
                                                      exact_stats,
                                                      tau_f_ratio_list,
                                                      tau_s_list)) {
            printResult(result);
        }
    }

    void runReplacementVariant(
            const ReplacementAblationConfig &config,
            const std::vector<std::pair<uint32_t, uint32_t>> &data,
            const ExactStats &exact_stats,
            const std::vector<float> &tau_f_ratio_list,
            const std::vector<double> &tau_s_list,
            ReplacementMechanism mechanism) {
        std::cout << "\n[Replacement strategy ablation] "
                  << ReplacementAblation::mechanismName(mechanism)
                  << ", total_kb=" << config.total_memory_kb
                  << ", pair_kb="
                  << config.total_memory_kb * config.pair_memory_ratio
                  << ", residence_kb="
                  << config.total_memory_kb * config.residence_memory_ratio
                  << ", cells_per_bucket=" << config.cell_num << std::endl;

        const ReplacementAblation experiment(config);
        for (const AblationResult &result : experiment.run(
                data, exact_stats, tau_f_ratio_list, tau_s_list, mechanism)) {
            printResult(result);
        }
    }

    void runCompensationVariant(
            const CompensationAblationConfig &config,
            const std::vector<std::pair<uint32_t, uint32_t>> &data,
            const ExactStats &exact_stats,
            const std::vector<float> &tau_f_ratio_list,
            const std::vector<double> &tau_s_list,
            CompensationMechanism compensation,
            ResidenceMechanism residence) {
        std::cout << "\n[Compensation ablation] "
                  << CompensationAblation::compensationName(compensation)
                  << "+" << CompensationAblation::residenceName(residence)
                  << ", total_kb=" << config.total_memory_kb
                  << ", pair_kb="
                  << config.total_memory_kb * config.pair_memory_ratio
                  << ", residence_kb=";
        if (residence == ResidenceMechanism::Ideal) {
            std::cout << "unlimited";
        } else {
            std::cout << config.total_memory_kb * config.residence_memory_ratio;
        }
        std::cout << ", cells_per_bucket=" << config.cell_num << std::endl;

        const CompensationAblation experiment(config);
        for (const AblationResult &result : experiment.run(
                data, exact_stats, tau_f_ratio_list, tau_s_list,
                compensation, residence)) {
            printResult(result);
        }
    }

    AblationConfig makeConfig(PairMechanism pair_mechanism,
                              ResidenceMechanism residence_mechanism,
                              double total_memory_kb,
                              double pair_memory_ratio,
                              double residence_memory_ratio,
                              uint32_t cell_num,
                              uint32_t reward,
                              uint32_t penalty,
                              double decay_base) {
        AblationConfig config;
        config.pair_mechanism = pair_mechanism;
        config.residence_mechanism = residence_mechanism;
        config.total_memory_kb = total_memory_kb;
        config.pair_memory_ratio = pair_memory_ratio;
        config.residence_memory_ratio = residence_memory_ratio;
        config.cell_num = cell_num;
        config.reward = reward;
        config.penalty = penalty;
        config.decay_base = decay_base;
        return config;
    }

} // namespace

int main(int argc, char **argv) {
    try {
        std::cout << "Krill Ablation Study starts ..." << std::endl;
        const auto start_time = std::chrono::steady_clock::now();

        const bool smoke_test = argc > 1 && std::string(argv[1]) == "--smoke";
        LoadedDataset loaded;
        if (smoke_test) {
            for (uint32_t flow_id = 1; flow_id <= 100; ++flow_id) {
                for (uint32_t element_id = 1; element_id <= 100; ++element_id) {
                    loaded.data.emplace_back(flow_id, element_id);
                    loaded.exact_stats.update(flow_id, element_id);
                }
                for (uint32_t repeat = 0; repeat < 100; ++repeat) {
                    loaded.data.emplace_back(flow_id, 1);
                    loaded.exact_stats.update(flow_id, 1);
                }
            }
            for (uint32_t repeat = 0; repeat < 2000; ++repeat) {
                loaded.data.emplace_back(1001, 1);
                loaded.exact_stats.update(1001, 1);
            }
            std::cout << "Synthetic smoke dataset loaded. data="
                      << loaded.data.size() << std::endl;
        } else {
            loaded = loadDatasetECommerceEvents();
            // loaded = loadDatasetTenrec();
            // loaded = loadDataSetStackOverflow();
        }

        // Every experiment uses one total-memory budget. Finite PairWitness/Bloom Filter
        // and FlowResidence receive 30% and 70%, respectively.
        const std::vector<double> total_memory_kb_list =
                smoke_test ? std::vector<double>{200}
                           : std::vector<double>{50, 100, 200, 400, 800, 1200};
        const double pair_memory_ratio = 0.3;
        const double residence_memory_ratio = 0.7;

        const std::vector<float> tau_f_ratio_list = {0.0001};
        const std::vector<double> tau_s_list = {1024};

        const uint32_t cell_num = 4;
        const uint32_t reward = 12;
        const uint32_t penalty = 1;
        const double decay_base = 1.08;

        const bool run_component_ablation = true;
        const bool run_challenge_strategy_ablation = false;
        const bool run_compensation_ablation = false;

        const std::vector<PairMechanism> pair_mechanisms = {
                PairMechanism::PairWitness,
                PairMechanism::BloomFilter,
                PairMechanism::Exact}; // Ideal setting, Exact checks for new/repeated, no memory limitations.

        if (run_component_ablation) {
            for (double total_memory_kb: total_memory_kb_list) {

                for (PairMechanism mechanism: pair_mechanisms) {
                    runVariant(makeConfig(mechanism,
                                          ResidenceMechanism::Ideal,
                                          total_memory_kb,
                                          pair_memory_ratio,
                                          residence_memory_ratio,
                                          cell_num,
                                          reward,
                                          penalty,
                                          decay_base),
                               loaded.data,
                               loaded.exact_stats,
                               tau_f_ratio_list,
                               tau_s_list,
                               "PairWitness ablation / ideal residence");
                    runVariant(makeConfig(mechanism,
                                          ResidenceMechanism::FlowResidence,
                                          total_memory_kb,
                                          pair_memory_ratio,
                                          residence_memory_ratio,
                                          cell_num,
                                          reward,
                                          penalty,
                                          decay_base),
                               loaded.data,
                               loaded.exact_stats,
                               tau_f_ratio_list,
                               tau_s_list,
                               "PairWitness ablation / finite residence");
                }

            }
        }

        if (run_challenge_strategy_ablation) {
            ReplacementAblationConfig replacement_config;
            replacement_config.pair_memory_ratio = pair_memory_ratio;
            replacement_config.residence_memory_ratio = residence_memory_ratio;
            replacement_config.cell_num = cell_num;
            replacement_config.reward = reward;
            replacement_config.penalty = penalty;
            replacement_config.decay_base = decay_base;

            const std::vector<ReplacementMechanism> replacement_mechanisms = {
                    ReplacementMechanism::ProbabilisticDecay,
                    ReplacementMechanism::DirectDecay,
                    ReplacementMechanism::DirectReplace};

            for (double total_memory_kb : total_memory_kb_list) {
                replacement_config.total_memory_kb = total_memory_kb;
                for (ReplacementMechanism mechanism : replacement_mechanisms) {
                    runReplacementVariant(replacement_config,
                                           loaded.data,
                                           loaded.exact_stats,
                                           tau_f_ratio_list,
                                           tau_s_list,
                                           mechanism);
                }
            }
        }

        if (run_compensation_ablation) {
            CompensationAblationConfig compensation_config;
            compensation_config.pair_memory_ratio = pair_memory_ratio;
            compensation_config.residence_memory_ratio = residence_memory_ratio;
            compensation_config.cell_num = cell_num;
            compensation_config.reward = reward;
            compensation_config.penalty = penalty;
            compensation_config.decay_base = decay_base;

            const std::vector<CompensationMechanism> compensation_mechanisms = {
                    CompensationMechanism::Compensated,
                    CompensationMechanism::Uncompensated};
            const std::vector<ResidenceMechanism> residence_mechanisms = {
                    ResidenceMechanism::Ideal,
                    ResidenceMechanism::FlowResidence
            };

            for (double total_memory_kb : total_memory_kb_list) {
                compensation_config.total_memory_kb = total_memory_kb;
                for (ResidenceMechanism residence : residence_mechanisms) {
                    for (CompensationMechanism compensation : compensation_mechanisms) {
                        runCompensationVariant(compensation_config,
                                                loaded.data,
                                                loaded.exact_stats,
                                                tau_f_ratio_list,
                                                tau_s_list,
                                                compensation,
                                                residence);
                    }
                }
            }
        }

        const auto end_time = std::chrono::steady_clock::now();
        const auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
                end_time - start_time).count();
        std::cout << "\nAblation experiments end. Elapsed time: "
                  << elapsed / 3600 << "h "
                  << (elapsed % 3600) / 60 << "m "
                  << elapsed % 60 << "s" << std::endl;
    } catch (const std::exception &exception) {
        std::cerr << "Ablation experiment failed: " << exception.what() << std::endl;
        return 1;
    }
    return 0;
}


