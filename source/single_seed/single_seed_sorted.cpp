#include "single_seed_sorted.hpp"
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

static uint64_t encode_seed_single(const string& s, size_t k) {
    if (s.size() < k)
        throw invalid_argument("String shorter than seed length");
    uint64_t seed = 0;
    for (size_t i = 0; i < k; i++)
        seed <<= 2, seed += encode_base(s[i]);
    return seed;
}

SingleSeedSortedIndex build_single_seed_sorted_index(
        const string& reference, size_t k) {

    if (k == 0 or k > 32 or k > reference.size())
        throw invalid_argument("Invalid k");

    SingleSeedSortedIndex index(k);
    index.positions.reserve(reference.size() - k + 1);

    uint64_t key = encode_seed_single(reference, k);
    index.positions.push_back({key, 0});

    uint64_t rolling_mask = (uint64_t(1) << (2 * (k - 1))) - 1;
    for (size_t i = 1; i + k <= reference.size(); i++) {
        key = (key & rolling_mask) << 2;
        key |= encode_base(reference[i + k - 1]);
        index.positions.push_back({key, i});
    }
    sort(index.positions.begin(), index.positions.end());
    return index;
}

vector<vector<size_t>> search_single_seed_sorted_index(
        const SingleSeedSortedIndex& index,
        const string& reference,
        const vector<string>& queries) {

    vector<vector<size_t>> matches(queries.size());
    for (size_t i = 0; i < queries.size(); i++) {
        uint64_t key = encode_seed_single(queries[i], index.k);
        auto it = lower_bound(index.positions.begin(), index.positions.end(),
                make_pair(key, size_t(0)));

        for (; it != index.positions.end() && it -> first == key; ++it) {
            size_t pos = it -> second;
            if (reference.compare(pos + index.k, queries[i].size() - index.k,
                    queries[i], index.k, queries[i].size() - index.k) == 0)
                matches[i].push_back(pos);
        }
    }

    return matches;
}

vector<vector<size_t>> solve_with_single_seed_sorted_index(
        const string& reference,
        const vector<string>& queries,
        size_t k) {

    SingleSeedSortedIndex index = build_single_seed_sorted_index(reference, k);
    vector<vector<size_t>> matches = search_single_seed_sorted_index(index, reference, queries);
    return matches;
}