
#include <algorithm>
#include "header/CountMin.h"

CountMin::CountMin(float memory_kb) {

    this->counter_bits = 32;
    this->depth = 4;

    uint32_t total_counters = static_cast<uint32_t>(
            std::max(1.0f, memory_kb) * 1024.0f * 8.0f /
            static_cast<float>(this->counter_bits));

    width = std::max(1, static_cast<int>(std::round(static_cast<double>(total_counters) /
                                                    static_cast<double>(depth))));

    counters = std::vector<std::vector<uint32_t>>(depth, std::vector<uint32_t>(width, 0));
}



void CountMin::update(const uint32_t flow_label, uint32_t weight) {
    uint32_t hash_value = 0;
    for (int i = 0; i < depth; i++) {
        MurmurHash3_x86_32(&flow_label, 4, i, &hash_value);
        uint32_t j = hash_value % width;
        const uint64_t new_value = static_cast<uint64_t>(counters[i][j]) +
                                   static_cast<uint64_t>(weight);
        counters[i][j] = new_value > UINT32_MAX
                         ? UINT32_MAX
                         : static_cast<uint32_t>(new_value);
    }
}


uint32_t CountMin::query(const uint32_t flow_label) {

    int min_value = 0x7FFFFFFF;
    uint32_t hash_value = 0;
    for (int i = 0; i < depth; i++) {
        MurmurHash3_x86_32(&flow_label, 4, i, &hash_value);
        uint32_t j = hash_value % width;
        int val = counters[i][j];
        if (val < min_value) {
            min_value = val;
        }
    }
    return min_value;
}

