#ifndef CU_SKETCH_H
#define CU_SKETCH_H

#include <cstdint>
#include <vector>

class CUSketch {
private:
    static constexpr uint32_t depth = 4;
    uint32_t width = 1;
    std::vector<std::vector<uint32_t>> counters;

public:
    explicit CUSketch(float memory_kb);

    void update(uint32_t flow_id, uint32_t weight = 1);
    uint32_t query(uint32_t flow_id) const;
};

#endif
