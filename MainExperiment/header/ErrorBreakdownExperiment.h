#ifndef ERROR_BREAKDOWN_EXPERIMENT_H
#define ERROR_BREAKDOWN_EXPERIMENT_H

#include <vector>

struct LoadedDataset;

struct ErrorBreakdownExperimentConfig {
    std::vector<float> memory_kb_list = {100, 200, 300, 400, 500};
    float tau_f_ratio = 0.0001f;
    double tau_s = 1024.0;
};

void runDetectionErrorBreakdown(
        const LoadedDataset& loaded,
        const ErrorBreakdownExperimentConfig& config);

#endif

