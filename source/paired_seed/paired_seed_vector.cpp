#include "paired_seed_vector.hpp"
#include <limits>
using namespace std;

static uint64_t encode_base(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: throw invalid_argument("Invalid DNA base");
    }
}

static uint64_t encode_seed_pair(const string& s, size_t m, size_t k) {
    uint64_t seed = 0;
    for (size_t i = 0; i < k; i++)
        seed <<= 2, seed += encode_base(s[i]);
    for (size_t i = 0; i < k; i++)
        seed <<= 2, seed += encode_base(s[m - k + i]);
    return seed;
}

PairedSeedVectorIndex build_paired_seed_vector_index(
        const string& reference,
        size_t m, size_t k) {

    if (k == 0 or k > 5 or 2 * k > m)
        throw invalid_argument("Invalid m or k (vector supports k <= 5)");
    if (reference.size() > numeric_limits<uint32_t>::max())
        throw invalid_argument("Reference too long for uint32_t positions");

    PairedSeedVectorIndex index(m, k);
    size_t keys = size_t(1) << (4 * k);
    index.offsets.resize(keys + 1, 0);
    if (m > reference.size())
        return index;

    size_t count = reference.size() - m + 1;
    index.positions.resize(count);

    uint64_t rolling_mask = (uint64_t(1) << (2 * (k - 1))) - 1;
    rolling_mask = (rolling_mask << 2 * k) | rolling_mask;

    // First pass: count positions for each key.
    uint64_t key = encode_seed_pair(reference, m, k);
    index.offsets[key + 1]++;
    for (size_t i = 1; i < count; i++) {
        key = (key & rolling_mask) << 2;
        key |= encode_base(reference[i + k - 1]) << (2 * k);
        key |= encode_base(reference[i + m - 1]);
        index.offsets[key + 1]++;
    }

    for (size_t i = 1; i < index.offsets.size(); i++)
        index.offsets[i] += index.offsets[i - 1];

    // Second pass: fill the positions in increasing order for each key.
    key = encode_seed_pair(reference, m, k);
    index.positions[index.offsets[key]++] = 0;
    for (size_t i = 1; i < count; i++) {
        key = (key & rolling_mask) << 2;
        key |= encode_base(reference[i + k - 1]) << (2 * k);
        key |= encode_base(reference[i + m - 1]);
        index.positions[index.offsets[key]++] = static_cast<uint32_t>(i);
    }

    // Restore the starting offsets after they were used as write cursors.
    for (size_t i = index.offsets.size() - 1; i > 0; i--)
        index.offsets[i] = index.offsets[i - 1];
    index.offsets[0] = 0;
    return index;
}

vector<vector<size_t>> search_paired_seed_vector_index(
        const PairedSeedVectorIndex& index,
        const string& reference,
        const vector<string>& queries) {

    vector<vector<size_t>> matches(queries.size());
    for (size_t i = 0; i < queries.size(); i++) {
        if (queries[i].size() != index.m)
            throw invalid_argument("Invalid query length");

        uint64_t key = encode_seed_pair(queries[i], index.m, index.k);
        for (size_t j = index.offsets[key]; j < index.offsets[key + 1]; j++) {
            size_t pos = index.positions[j];
            if (reference.compare(pos + index.k, index.m - 2 * index.k,
                    queries[i], index.k, index.m - 2 * index.k) == 0)
                matches[i].push_back(pos);
        }
    }

    return matches;
}

vector<vector<size_t>> solve_with_paired_seed_vector_index(
        const string& reference,
        const vector<string>& queries,
        size_t m, size_t k) {

    PairedSeedVectorIndex index = build_paired_seed_vector_index(reference, m, k);
    vector<vector<size_t>> matches = search_paired_seed_vector_index(index, reference, queries);
    return matches;
}
