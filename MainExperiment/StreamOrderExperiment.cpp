#include "header/StreamOrderExperiment.h"

#include "header/Krill.h"
#include "header/ExactStats.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>

namespace {

using PairData = std::vector<std::pair<uint32_t, uint32_t>>;

void shuffleWithinBlocks(PairData& data,
                         std::size_t block_size,
                         uint32_t seed) {
    if (block_size == 0) {
        throw std::invalid_argument("block_size must be positive.");
    }
    std::mt19937 generator(seed);
    for (std::size_t begin = 0; begin < data.size(); begin += block_size) {
        const std::size_t end = std::min(begin + block_size, data.size());
        std::shuffle(data.begin() + static_cast<std::ptrdiff_t>(begin),
                     data.begin() + static_cast<std::ptrdiff_t>(end),
                     generator);
    }
}

} // namespace

void runStreamOrderRobustness(const LoadedDataset& loaded,
                              const StreamOrderExperimentConfig& config) {
    if (loaded.data.empty()) {
        throw std::invalid_argument("The loaded stream is empty.");
    }
    if (config.memory_kb_list.empty()) {
        throw std::invalid_argument("memory_kb_list must not be empty.");
    }
    if (config.block_size == 0) {
        throw std::invalid_argument("block_size must be positive.");
    }

    const auto runOne = [&](const std::string& order,
                            uint32_t seed,
                            std::size_t block_size,
                            const PairData& ordered_data,
                            float memory_kb) {
        Krill sketch(memory_kb,
                             config.pair_witness_ratio,
                             config.cell_num,
                             config.reward,
                             config.penalty);
        const KrillEvalResult result =
                sketch.evaluateForOrder(ordered_data,
                                        loaded.exact_stats,
                                        config.tau_f_ratio,
                                        config.tau_s);
        std::cout << std::fixed << std::setprecision(6)
                  << "Order=" << order
                  << ", seed=" << seed
                  << ", memory_kb=" << memory_kb
                  << ", precision=" << result.precision
                  << ", recall=" << result.recall
                  << ", F1=" << result.f1
                  << ", freq_ARE=" << result.frequency_are
                  << ", dss_ARE=" << result.dss_are
                  << ", throughput_Mips=" << result.update_mips
                  << std::endl;
    };

    std::cout << "Stream-order robustness starts: dataset="
              << config.dataset_name
              << ", items=" << loaded.data.size()
              << ", tau_f_ratio=" << config.tau_f_ratio
              << ", tau_S=" << config.tau_s
              << std::endl;

    for (float memory_kb : config.memory_kb_list) {
        runOne("Original", 0, 0, loaded.data, memory_kb);

        PairData reversed = loaded.data;
        std::reverse(reversed.begin(), reversed.end());
        runOne("Reverse", 0, 0, reversed, memory_kb);

        PairData random_order = loaded.data;
        std::mt19937 random_generator(config.random_seed);
        std::shuffle(random_order.begin(), random_order.end(), random_generator);
        runOne("Random", config.random_seed, 0, random_order, memory_kb);

        PairData block_shuffle = loaded.data;
        shuffleWithinBlocks(block_shuffle,
                            config.block_size,
                            config.block_shuffle_seed);
        runOne("BlockShuffle",
               config.block_shuffle_seed,
               config.block_size,
               block_shuffle,
               memory_kb);
    }
}


