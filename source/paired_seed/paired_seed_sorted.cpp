#include "paired_seed_sorted.hpp"
#include <algorithm>
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

PairedSeedSortedIndex build_paired_seed_sorted_index(
        const string& reference,
        size_t m, size_t k) {

    if (k == 0 or k > 16 or 2 * k > m)
        throw invalid_argument("Invalid m or k");

    PairedSeedSortedIndex index(m, k);
    if (m > reference.size())
        return index;

    index.positions.reserve(reference.size() - m + 1);
    uint64_t key = encode_seed_pair(reference, m, k);
    index.positions.push_back({key, 0});

    uint64_t rolling_mask = (uint64_t(1) << (2 * (k - 1))) - 1;
    rolling_mask = (rolling_mask << 2 * k) | rolling_mask;

    for (size_t i = 1; i + m <= reference.size(); i++) {
        key = (key & rolling_mask) << 2;
        key |= encode_base(reference[i + k - 1]) << 2 * k;
        key |= encode_base(reference[i + m - 1]);
        index.positions.push_back({key, i});
    }
    sort(index.positions.begin(), index.positions.end());
    return index;
}

vector<vector<size_t>> search_paired_seed_sorted_index(
        const PairedSeedSortedIndex& index,
        const string& reference,
        const vector<string>& queries) {

    vector<vector<size_t>> matches(queries.size());
    for (size_t i = 0; i < queries.size(); i++) {
        if (queries[i].size() != index.m)
            throw invalid_argument("Invalid query length");

        uint64_t key = encode_seed_pair(queries[i], index.m, index.k);
        auto it = lower_bound(index.positions.begin(), index.positions.end(),
                make_pair(key, size_t(0)));

        for (; it != index.positions.end() && it -> first == key; ++it) {
            size_t pos = it -> second;
            if (reference.compare(pos + index.k, index.m - 2 * index.k,
                    queries[i], index.k, index.m - 2 * index.k) == 0)
                matches[i].push_back(pos);
        }
    }

    return matches;
}

vector<vector<size_t>> solve_with_paired_seed_sorted_index(
        const string& reference,
        const vector<string>& queries,
        size_t m, size_t k) {

    PairedSeedSortedIndex index = build_paired_seed_sorted_index(reference, m, k);
    vector<vector<size_t>> matches = search_paired_seed_sorted_index(index, reference, queries);
    return matches;
}
