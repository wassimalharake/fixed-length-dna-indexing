#include "single_seed_vector.hpp"
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

static uint64_t encode_seed_single(const string& s, size_t k) {
    if (s.size() < k)
        throw invalid_argument("String shorter than seed length");
    uint64_t seed = 0;
    for (size_t i = 0; i < k; i++)
        seed <<= 2, seed += encode_base(s[i]);
    return seed;
}

SingleSeedVectorIndex build_single_seed_vector_index(
        const string& reference, size_t k) {

    if (k == 0 or k > 10 or k > reference.size())
        throw invalid_argument("Invalid k (vector supports k <= 10)");
    if (reference.size() > UINT32_MAX)
        throw invalid_argument("Reference too long for uint32_t positions");

    SingleSeedVectorIndex index(k);
    size_t keys = size_t(1) << (2 * k);
    index.offsets.resize(keys + 1, 0);

    size_t count = reference.size() - k + 1;
    index.positions.resize(count);

    uint64_t rolling_mask = (uint64_t(1) << (2 * (k - 1))) - 1;
    uint64_t key = encode_seed_single(reference, k);
    index.offsets[key + 1]++;
    for (size_t i = 1; i < count; i++) {
        key = (key & rolling_mask) << 2;
        key |= encode_base(reference[i + k - 1]);
        index.offsets[key + 1]++;
    }

    for (size_t i = 1; i < index.offsets.size(); i++)
        index.offsets[i] += index.offsets[i - 1];

    key = encode_seed_single(reference, k);
    index.positions[index.offsets[key]++] = 0;
    for (size_t i = 1; i < count; i++) {
        key = (key & rolling_mask) << 2;
        key |= encode_base(reference[i + k - 1]);
        index.positions[index.offsets[key]++] = static_cast<uint32_t>(i);
    }

    for (size_t i = index.offsets.size() - 1; i > 0; i--)
        index.offsets[i] = index.offsets[i - 1];
    index.offsets[0] = 0;
    return index;
}

vector<vector<size_t>> search_single_seed_vector_index(
        const SingleSeedVectorIndex& index,
        const string& reference,
        const vector<string>& queries) {

    vector<vector<size_t>> matches(queries.size());
    for (size_t i = 0; i < queries.size(); i++) {
        uint64_t key = encode_seed_single(queries[i], index.k);
        for (size_t j = index.offsets[key]; j < index.offsets[key + 1]; j++) {
            size_t pos = index.positions[j];
            if (reference.compare(pos + index.k, queries[i].size() - index.k,
                    queries[i], index.k, queries[i].size() - index.k) == 0)
                matches[i].push_back(pos);
        }
    }

    return matches;
}

vector<vector<size_t>> solve_with_single_seed_vector_index(
        const string& reference,
        const vector<string>& queries,
        size_t k) {

    SingleSeedVectorIndex index = build_single_seed_vector_index(reference, k);
    vector<vector<size_t>> matches = search_single_seed_vector_index(index, reference, queries);
    return matches;
}
