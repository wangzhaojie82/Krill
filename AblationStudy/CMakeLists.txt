cmake_minimum_required(VERSION 3.26)
project(KrillAblationStudy LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

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
        PROJECT_ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
        $<$<CONFIG:Release>:NDEBUG>)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(KrillAblationStudy PRIVATE
            $<$<CONFIG:Release>:-O3>
            $<$<CONFIG:Release>:-march=native>
            $<$<CONFIG:Release>:-funroll-loops>
            $<$<CONFIG:Release>:-ffast-math>)
elseif(MSVC)
    target_compile_options(KrillAblationStudy PRIVATE
            /utf-8
            $<$<CONFIG:Release>:/O2>
            $<$<CONFIG:Release>:/Ob2>
            $<$<CONFIG:Release>:/Oi>
            $<$<CONFIG:Release>:/Ot>
            $<$<CONFIG:Release>:/fp:fast>)
endif()
