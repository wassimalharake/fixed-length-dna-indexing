#pragma once
#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>
#include <stdexcept>

struct SingleSeedVectorIndex {
    std::vector<uint32_t> offsets;
    std::vector<uint32_t> positions;
    size_t k;
    SingleSeedVectorIndex(size_t _k) : k(_k) {}
};

SingleSeedVectorIndex build_single_seed_vector_index(
    const std::string& reference,
    size_t k
);

std::vector<std::vector<size_t>> search_single_seed_vector_index(
    const SingleSeedVectorIndex& index,
    const std::string& reference,
    const std::vector<std::string>& queries
);

std::vector<std::vector<size_t>> solve_with_single_seed_vector_index(
    const std::string& reference,
    const std::vector<std::string>& queries,
    size_t k
);
