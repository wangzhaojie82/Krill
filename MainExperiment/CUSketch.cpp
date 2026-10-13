#include "header/CUSketch.h"

#include "header/MurmurHash3.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

CUSketch::CUSketch(float memory_kb) {
    const double total_bits =
            std::max(1.0, static_cast<double>(memory_kb) * 1024.0 * 8.0);
    // Logical storage: four rows of 32-bit counters.
    width = static_cast<uint32_t>(std::max(1.0,
            std::floor(total_bits / (depth * 32.0))));
    counters.assign(depth, std::vector<uint32_t>(width, 0));
}

void CUSketch::update(uint32_t flow_id, uint32_t weight) {
    struct CounterIndex {
        uint32_t row = 0;
        uint32_t column = 0;
    };

    std::array<CounterIndex, depth> minimum_indices{};
    uint32_t minimum_count = 0;
    uint32_t minimum = std::numeric_limits<uint32_t>::max();

    for (uint32_t row = 0; row < depth; ++row) {
        uint32_t hash_value = 0;
        MurmurHash3_x86_32(&flow_id, sizeof(flow_id), row, &hash_value);
        const uint32_t column = hash_value % width;
        const uint32_t value = counters[row][column];
        if (value < minimum) {
            minimum = value;
            minimum_count = 0;
        }
        if (value == minimum) {
            minimum_indices[minimum_count++] = {row, column};
        }
    }

    for (uint32_t i = 0; i < minimum_count; ++i) {
        const CounterIndex index = minimum_indices[i];
        const uint64_t updated = static_cast<uint64_t>(
                counters[index.row][index.column]) + weight;
        counters[index.row][index.column] = static_cast<uint32_t>(std::min<uint64_t>(
                updated, std::numeric_limits<uint32_t>::max()));
    }
}

uint32_t CUSketch::query(uint32_t flow_id) const {
    uint32_t minimum = std::numeric_limits<uint32_t>::max();
    for (uint32_t row = 0; row < depth; ++row) {
        uint32_t hash_value = 0;
        MurmurHash3_x86_32(&flow_id, sizeof(flow_id), row, &hash_value);
        const uint32_t column = hash_value % width;
        minimum = std::min(minimum, counters[row][column]);
    }
    return minimum;
}
