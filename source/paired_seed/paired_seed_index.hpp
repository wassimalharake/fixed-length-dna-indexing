#pragma once
#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>
#include <unordered_map>
#include <stdexcept>


struct PairedSeedIndex {
    std::unordered_map<uint64_t, std::vector<size_t>> positions;
    size_t m, k;
    PairedSeedIndex(size_t _m, size_t _k) : m(_m), k(_k) {}
};

PairedSeedIndex build_paired_seed_index(
    const std::string& reference,
    size_t m, size_t k
);

std::vector<std::vector<size_t>> search_paired_seed_index(
    const PairedSeedIndex& index,
    const std::string& reference,
    const std::vector<std::string>& queries
);

std::vector<std::vector<size_t>> solve_with_paired_seed_index(
    const std::string& reference,
    const std::vector<std::string>& queries,
    size_t m, size_t k
);