cmake_minimum_required(VERSION 3.20)
project(cmp355_research LANGUAGES CXX)

list(APPEND CMAKE_PREFIX_PATH
    "${CMAKE_CURRENT_SOURCE_DIR}/../seqan3/cmake")

find_package(seqan3 3.0 REQUIRED)

add_executable(main
    main.cpp
    fm_index/fm_index.cpp
    single_seed/single_seed_index.cpp
    single_seed/single_seed_sorted.cpp
    single_seed/single_seed_vector.cpp
    paired_seed/paired_seed_index.cpp
    paired_seed/paired_seed_sorted.cpp
    paired_seed/paired_seed_vector.cpp
)
target_link_libraries(main PRIVATE seqan3::seqan3)
target_include_directories(main SYSTEM PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/../seqan3/include"
)