#include "single_seed_index.hpp"
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

uint64_t encode_seed_single(const string& s, size_t k) {
    if (s.size() < k)
        throw invalid_argument("String shorter than seed length");
    uint64_t seed = 0;
    for (size_t i = 0; i < k; i++)
        seed <<= 2, seed += encode_base(s[i]);
    return seed;
}

SingleSeedIndex build_single_seed_index(
        const string& reference, size_t k) {

    if (k == 0 or k > 32 or k > reference.size())
        throw invalid_argument("Invalid k");

    SingleSeedIndex index(k);

    uint64_t key = encode_seed_single(reference, k);
    index.positions[key].push_back(0);
    
    uint64_t rolling_mask = (uint64_t(1) << (2 * (k - 1))) - 1;
    for (size_t i = 1; i + k <= reference.size(); i++) {
        key = (key & rolling_mask) << 2;
        key |= encode_base(reference[i + k - 1]);
        index.positions[key].push_back(i);
    }
    return index;
}

vector<vector<size_t>> search_single_seed_index(
        const SingleSeedIndex& index,
        const string& reference,
        const vector<string>& queries) {

    vector<vector<size_t>> matches(queries.size());
    for (size_t i = 0; i < queries.size(); i++) {
        uint64_t seed = encode_seed_single(queries[i], index.k);
        auto it = index.positions.find(seed);
        if (it == index.positions.end())
            continue;

        for (size_t pos : it -> second) {
            if (reference.compare(pos + index.k, queries[i].size() - index.k,
                    queries[i], index.k, queries[i].size() - index.k) == 0)
                matches[i].push_back(pos);
        }
    }

    return matches;
}

vector<vector<size_t>> solve_with_single_seed_index(
        const string& reference,
        const vector<string>& queries,
        size_t k) {

    SingleSeedIndex index = build_single_seed_index(reference, k);
    vector<vector<size_t>> matches = search_single_seed_index(index, reference, queries);
    return matches;
}