# Krill Main Experiment

This repository contains the standalone experiment code for **Krill** and the baselines used in the experimental evaluation.

## Project Structure

```text
MainExperiment/
|-- CMakeLists.txt
|-- README.md
|-- main.cpp
|-- DatasetLoader.cpp
|-- Krill.cpp
|-- BloomFilterBaseline.cpp
|-- PISketch.cpp
|-- CountMinVHLLBaseline.cpp
|-- FLASKSketch.cpp
|-- StreamOrderExperiment.cpp
|-- ErrorBreakdownExperiment.cpp
|-- CountMin.cpp
|-- vHLL.cpp
|-- MurmurHash3.cpp
`-- header/
    |-- Krill.h
    |-- BloomFilterBaseline.h
    |-- PISketch.h
    |-- CountMinVHLLBaseline.h
    |-- FLASKSketch.h
    |-- StreamOrderExperiment.h
    |-- ErrorBreakdownExperiment.h
    |-- DatasetLoader.h
    |-- ExactStats.h
    |-- FlowEstimate.h
    |-- CountMin.h
    |-- vHLL.h
    `-- MurmurHash3.h
```

## Environment

Recommended environment:

- CMake 3.26 or later
- C++17 compiler
- Windows with MSVC, MinGW/GCC, or Clang
- Sufficient memory and disk bandwidth for the selected dataset

## Build

From `Krill/MainExperiment`, configure and build with CMake. The release optimization flags are already defined in `CMakeLists.txt`.

```powershell
cmake -S . -B cmake-build-release
cmake --build cmake-build-release -j
```

Run the executable:

```powershell
.\cmake-build-release\Krill.exe
```

## Datasets

Dataset files are not included. Download the datasets separately and update the file paths in `DatasetLoader.cpp` if needed.

Select one dataset loader in `main.cpp`:

```cpp
auto loaded = loadDatasetECommerceEvents();
//    auto loaded = loadDatasetTenrec();
//    auto loaded = loadDataSetStackOverflow();
```

The available dataset entries for reproduction are:

1. E-Commerce Events History in Cosmetics. 

> ECommerce Events History contains behavior data for 5 months (Oct 2019 – Feb 2020) from a medium cosmetics online store. Each row in the file represents an event. All events are related to products and users. Data collected by Open CDP project. 
> 
> Download link placeholder:  https://www.kaggle.com/datasets/mkechinov/ecommerce-events-history-in-cosmetics-shop/

2. Tenrec recommendation dataset

> Tenrec is a large-scale multipurpose benchmark dataset for recommender systems where data was collected from two feeds (articles and videos) recommendation platforms.
> 
> Download link: https://tenrec0.github.io/
> 
> Citation：
> 
> Yuan G, Yuan F, Li Y, et al. Tenrec: A large-scale multipurpose benchmark dataset for recommender systems[J]. Advances in Neural Information Processing Systems, 2022, 35: 11480-11493.

3. SNAP Stack Overflow temporal network

> This is a temporal network of interactions on the stack exchange web site Stack Overflow. 
> 
> Download link: https://snap.stanford.edu/data/sx-stackoverflow.html
> 
> Citation: 
> 
> Ashwin Paranjape, Austin R. Benson, and Jure Leskovec. "Motifs in Temporal Networks." In Proceedings of the Tenth ACM International Conference on Web Search and Data Mining, 2017.

Each loader converts raw records into (flow, element) pairs and updates exact statistics for ground truth.

## Experiment Configuration

The experiment switches are in `main.cpp`:

```cpp
const bool run_baseline_compare = true;
const bool run_parameter_sweep = false;
const bool run_order_robustness = false;
const bool run_error_breakdown = false;
```

Use one experiment switch at a time for clean logs.

### Baseline Comparison

Set:

```cpp
const bool run_baseline_compare = true;
```

This runs Krill and the baselines over `memory_kb_list` and prints metrics to the console.

### Parameter Sensitivity

Set:

```cpp
const bool run_parameter_sweep = true;
```

Edit the parameter lists in `main.cpp` as needed. Results are printed to the console with the scanned parameters and the following metrics.

### Stream Order Robustness

Set:

```cpp
const bool run_order_robustness = true;
```

Configuration is in the `StreamOrderExperimentConfig` block in `main.cpp`.

### Error Breakdown

Set:

```cpp
const bool run_error_breakdown = false;
```

This prints detection error categories for Krill and baseline to the console. Adjust parameters in the  `ErrorBreakdownExperimentConfig` block.

## Reproducing Results

A typical workflow is:

1. Download the selected dataset.
2. Update the corresponding file paths in `DatasetLoader.cpp`.
3. Select the dataset loader in `main.cpp`.
4. Select the experiment switch in `main.cpp`.
5. Build in Release mode with CMake.
6. Run `Krill.exe` and get the console output.
