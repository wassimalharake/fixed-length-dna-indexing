#pragma once
#include <cstddef>
#include <string>
#include <vector>
#include <algorithm>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/fm_index/all.hpp>

using FMIndex = seqan3::fm_index<seqan3::dna4, seqan3::text_layout::single>;

FMIndex build_fm_index(const std::string& reference);

std::vector<std::vector<size_t>> search_fm_index(
    const FMIndex& index,
    const std::vector<std::string>& queries
);

std::vector<std::vector<size_t>> solve_with_fm_index(
    const std::string& reference,
    const std::vector<std::string>& queries
);