cmake_minimum_required(VERSION 3.26)
project(Krill LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -funroll-loops -ffast-math -DNDEBUG")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")

add_executable(Krill
        main.cpp
        DatasetLoader.cpp
        Krill.cpp
        BloomFilterBaseline.cpp
        PISketch.cpp
        CountMinVHLLBaseline.cpp
        FLASKSketch.cpp
        CocoSketch.cpp
        DUET.cpp
        CUSketch.cpp
        StreamOrderExperiment.cpp
        ErrorBreakdownExperiment.cpp
        CountMin.cpp
        MurmurHash3.cpp
        header/DatasetLoader.h
        header/Krill.h
        header/BloomFilterBaseline.h
        header/PISketch.h
        header/CountMinVHLLBaseline.h
        header/FLASKSketch.h
        header/CocoSketch.h
        header/DUET.h
        header/CUSketch.h
        header/StreamOrderExperiment.h
        header/ErrorBreakdownExperiment.h
        header/FlowEstimate.h
        header/CountMin.h
        header/ExactStats.h
        header/MurmurHash3.h
        header/vHLL.h
        vHLL.cpp)

target_include_directories(Krill PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/header)

target_compile_definitions(Krill PRIVATE
        PROJECT_ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
