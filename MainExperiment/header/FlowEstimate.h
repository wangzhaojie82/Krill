#ifndef FLOW_ESTIMATE_H
#define FLOW_ESTIMATE_H

#include <cstdint>

struct FlowEstimate {
    uint32_t flow_id = 0;
    uint32_t frequency = 0;
    double dss = 0.0;
};

#endif

