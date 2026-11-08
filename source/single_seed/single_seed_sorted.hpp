#pragma once
#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>
#include <utility>
#include <stdexcept>

struct SingleSeedSortedIndex {
    std::vector<std::pair<uint64_t, size_t>> positions;
    size_t k;
    SingleSeedSortedIndex(size_t _k) : k(_k) {}
};

SingleSeedSortedIndex build_single_seed_sorted_index(
    const std::string& reference,
    size_t k
);

std::vector<std::vector<size_t>> search_single_seed_sorted_index(
    const SingleSeedSortedIndex& index,
    const std::string& reference,
    const std::vector<std::string>& queries
);

std::vector<std::vector<size_t>> solve_with_single_seed_sorted_index(
    const std::string& reference,
    const std::vector<std::string>& queries,
    size_t k
);
