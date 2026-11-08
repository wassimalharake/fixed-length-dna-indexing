#pragma once
#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>
#include <unordered_map>
#include <stdexcept>

struct SingleSeedIndex {
    std::unordered_map<uint64_t, std::vector<size_t>> positions;
    size_t k;
    SingleSeedIndex(size_t _k) : k(_k) {}
};

SingleSeedIndex build_single_seed_index(
    const std::string& reference,
    size_t k
);

std::vector<std::vector<size_t>> search_single_seed_index(
    const SingleSeedIndex& index,
    const std::string& reference,
    const std::vector<std::string>& queries
);

std::vector<std::vector<size_t>> solve_with_single_seed_index(
    const std::string& reference,
    const std::vector<std::string>& queries,
    size_t k
);