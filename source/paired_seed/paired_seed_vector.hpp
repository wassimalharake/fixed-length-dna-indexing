#pragma once
#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>
#include <stdexcept>

struct PairedSeedVectorIndex {
    std::vector<uint32_t> offsets;
    std::vector<uint32_t> positions;
    size_t m, k;
    PairedSeedVectorIndex(size_t _m, size_t _k) : m(_m), k(_k) {}
};

PairedSeedVectorIndex build_paired_seed_vector_index(
    const std::string& reference,
    size_t m, size_t k
);

std::vector<std::vector<size_t>> search_paired_seed_vector_index(
    const PairedSeedVectorIndex& index,
    const std::string& reference,
    const std::vector<std::string>& queries
);

std::vector<std::vector<size_t>> solve_with_paired_seed_vector_index(
    const std::string& reference,
    const std::vector<std::string>& queries,
    size_t m, size_t k
);
