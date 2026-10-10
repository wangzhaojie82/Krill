# Krill Ablation Study

This repository contains the standalone ablation-study code for **Krill**.

## Project Structure

```text
AblationStudy/
|-- CMakeLists.txt
|-- README.md
|-- main.cpp
|-- DatasetLoader.cpp
|-- AblationSketch.cpp
|-- ReplacementAblation.cpp
|-- CompensationAblation.cpp
|-- MurmurHash3.cpp
`-- header/
    |-- AblationSketch.h
    |-- ReplacementAblation.h
    |-- CompensationAblation.h
    |-- DatasetLoader.h
    |-- ExactStats.h
    `-- MurmurHash3.h
```

## Build

From `\Krill\AblationStudy`:

```powershell
cmake -S . -B cmake-build-release -DCMAKE_BUILD_TYPE=Release
cmake --build cmake-build-release --config Release -j
```

Run the executable:

```powershell
.\cmake-build-release\KrillAblationStudy.exe
```

If CMake uses a Visual Studio generator, the executable may be under:

```powershell
.\cmake-build-release\Release\KrillAblationStudy.exe
```

## Datasets

Dataset files are not included. Download the datasets separately and update the file paths in `DatasetLoader.cpp` if needed.

Select one dataset loader in `main.cpp`:

```cpp
loaded = loadDatasetECommerceEvents();
// loaded = loadDatasetTenrec();
// loaded = loadDataSetStackOverflow();
```

The dataset sources are the same as in `../MainExperiment/README.md`.

## Experiment Configuration

The experiment switches are in `main.cpp`:

```cpp
const bool run_component_ablation = true;
const bool run_challenge_strategy_ablation = false;
const bool run_compensation_ablation = false;
```

Use one experiment switch at a time for clean logs.

Main parameters are also configured in `main.cpp`. Results are printed to the console with the ablation variant, memory settings, F1, precision, recall, ARE, and throughput.

## Reproducing Results

A typical workflow is:

1. Download the selected dataset.
2. Update the corresponding file paths in `DatasetLoader.cpp`.
3. Select the dataset loader in `main.cpp`.
4. Select the ablation switch in `main.cpp`.
5. Build in Release mode with CMake.
6. Run `KrillAblationStudy.exe` and save the console output.

A small synthetic smoke test is available:

```powershell
.\cmake-build-release\KrillAblationStudy.exe --smoke
```
