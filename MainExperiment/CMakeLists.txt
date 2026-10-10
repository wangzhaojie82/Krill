cmake_minimum_required(VERSION 3.26)
project(Krill LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

add_executable(Krill
        main.cpp
        DatasetLoader.cpp
        Krill.cpp
        BloomFilterBaseline.cpp
        PISketch.cpp
        CountMinVHLLBaseline.cpp
        FLASKSketch.cpp
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
        PROJECT_ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
        $<$<CONFIG:Release>:NDEBUG>)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(Krill PRIVATE
            $<$<CONFIG:Release>:-O3>
            $<$<CONFIG:Release>:-march=native>
            $<$<CONFIG:Release>:-funroll-loops>
            $<$<CONFIG:Release>:-ffast-math>)
elseif(MSVC)
    target_compile_options(Krill PRIVATE
            /utf-8
            $<$<CONFIG:Release>:/O2>
            $<$<CONFIG:Release>:/Ob2>
            $<$<CONFIG:Release>:/Oi>
            $<$<CONFIG:Release>:/Ot>
            $<$<CONFIG:Release>:/fp:fast>)
endif()
