cmake_minimum_required(VERSION 3.26)
project(KrillAblationStudy LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -funroll-loops -ffast-math -DNDEBUG")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")

add_executable(KrillAblationStudy
        main.cpp
        DatasetLoader.cpp
        AblationSketch.cpp
        ReplacementAblation.cpp
        CompensationAblation.cpp
        MurmurHash3.cpp
        header/DatasetLoader.h
        header/AblationSketch.h
        header/ReplacementAblation.h
        header/CompensationAblation.h
        header/ExactStats.h
        header/MurmurHash3.h)

target_include_directories(KrillAblationStudy PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/header)

target_compile_definitions(KrillAblationStudy PRIVATE
        PROJECT_ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}")