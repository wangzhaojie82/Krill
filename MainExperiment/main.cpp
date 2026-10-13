#include "header/Krill.h"
#include "header/DatasetLoader.h"
#include "header/BloomFilterBaseline.h"
#include "header/PISketch.h"
#include "header/CountMinVHLLBaseline.h"
#include "header/FLASKSketch.h"
#include "header/CocoSketch.h"
#include "header/DUET.h"
#include "header/StreamOrderExperiment.h"
#include "header/ErrorBreakdownExperiment.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <cstdint>
#include <vector>

int main() {
    std::cout << "Krill experiment starts ..." << std::endl;
    const auto start_time = std::chrono::steady_clock::now();

    auto loaded = loadDatasetECommerceEvents();
//    auto loaded = loadDatasetTenrec();
//    auto loaded = loadDataSetStackOverflow();

    const std::vector<float> memory_kb_list = {100, 200, 300, 400, 500};

    const std::vector<float> pair_witness_ratio_list = {0.3};
    const std::vector<uint32_t> cell_num_list = {4};
    const std::vector<uint32_t> reward_list = {12};
    const std::vector<uint32_t> penalty_list = {1};
    const std::vector<float> tau_f_ratio_list = {0.0001};
    const std::vector<double> tau_s_list = {1024};

    // Experimental control switch
    const bool run_baseline_compare = true; // overall compare
    const bool run_parameter_sweep = false; // param sensitivity
    const bool run_order_robustness = false; // robustness to stream order
    const bool run_error_breakdown = false; // error breakdown

    if (run_baseline_compare) {
        for (float memory_kb : memory_kb_list) {

                    std::cout << "\nBaseline Comp: memory_kb=" << memory_kb << std::endl;

                    Krill sketch(memory_kb);
                    sketch.eval(loaded.data, loaded.exact_stats, tau_f_ratio_list, tau_s_list);

                    BloomFilterBaseline bloom_filter_baseline_sketch(memory_kb);
                    bloom_filter_baseline_sketch.eval(loaded.data, loaded.exact_stats, tau_f_ratio_list, tau_s_list);

                    PISketch piSketch(memory_kb);
                    piSketch.eval(loaded.data, loaded.exact_stats, tau_f_ratio_list,tau_s_list);

                    CountMinVHLLBaseline cmVhllBaseline(memory_kb);
                    cmVhllBaseline.eval(loaded.data, loaded.exact_stats, tau_f_ratio_list, tau_s_list);

                    FLASKSketch flaskSketch(memory_kb);
                    flaskSketch.eval(loaded.data, loaded.exact_stats, tau_f_ratio_list, tau_s_list);

                    CocoSketch cocoSketch(memory_kb);
                    cocoSketch.eval(loaded.data, loaded.exact_stats, tau_f_ratio_list, tau_s_list);

                    DUET duet(memory_kb);
                    duet.eval(loaded.data, loaded.exact_stats, tau_f_ratio_list, tau_s_list);

        }
    }

    if (run_parameter_sweep) {
        uint64_t run_index = 0;
        const uint64_t total_run_count =
                static_cast<uint64_t>(memory_kb_list.size()) *
                static_cast<uint64_t>(pair_witness_ratio_list.size()) *
                static_cast<uint64_t>(cell_num_list.size()) *
                static_cast<uint64_t>(reward_list.size()) *
                static_cast<uint64_t>(penalty_list.size());

        for (float memory_kb : memory_kb_list) {
            for (float pair_witness_ratio : pair_witness_ratio_list) {
                for (uint32_t reward : reward_list) {
                    for (uint32_t penalty : penalty_list) {
                        for (uint32_t cell_num : cell_num_list) {

                            ++run_index;
                            std::cout << "[" << run_index << "/" << total_run_count
                                      << "] memory_kb=" << memory_kb
                                      << ", pair_witness_ratio=" << pair_witness_ratio
                                      << ", reward=" << reward
                                      << ", penalty=" << penalty
                                      << ", cell_num=" << cell_num
                                      << std::endl;

                            Krill sketch(memory_kb,
                                                pair_witness_ratio,
                                                cell_num,
                                                reward,
                                                penalty);
                            sketch.evalForSweep(loaded.data,
                                                loaded.exact_stats,
                                                tau_f_ratio_list,
                                                tau_s_list);
                        }
                    }
                }
            }
        }
    }

    if (run_order_robustness) {
        StreamOrderExperimentConfig config;
        config.dataset_name = "ECommerceEvents";
        config.memory_kb_list = {100, 200, 300, 400, 500};
        config.pair_witness_ratio = 0.3f;
        config.cell_num = 4;
        config.reward = 12;
        config.penalty = 1;
        config.tau_f_ratio = 0.0001f;
        config.tau_s = 1024.0;
        config.random_seed = 2001;
        config.block_shuffle_seed = 1809;
        config.block_size = 1000000;

        runStreamOrderRobustness(loaded, config);
    }

    if (run_error_breakdown) {
        ErrorBreakdownExperimentConfig config;
        config.memory_kb_list = {200};
        config.tau_f_ratio = 0.0001f;
        config.tau_s = 1024.0;

        runDetectionErrorBreakdown(loaded, config);
    }

    const auto default_end_time = std::chrono::steady_clock::now();
    const auto default_elapsed_time =
            std::chrono::duration_cast<std::chrono::seconds>(default_end_time - start_time).count();
    std::cout << "\nExperiment ends. Elapsed time: "
              << default_elapsed_time / 3600 << "h "
              << (default_elapsed_time % 3600) / 60 << "m "
              << default_elapsed_time % 60 << "s" << std::endl;
    return 0;
}



