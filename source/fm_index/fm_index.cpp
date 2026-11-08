#include "fm_index.hpp"
using namespace std;
using namespace seqan3;

dna4_vector to_dna4_vector(const string& s) {
    dna4_vector v;
    v.reserve(s.size());
    for (const char& c : s) {
        dna4 base;
        assign_char_to(c, base);
        v.push_back(base);
    }
    return v;
}

FMIndex build_fm_index(const string& reference) {
    dna4_vector genome = dna4_vector(to_dna4_vector(reference));
    return fm_index{genome};
}

vector<vector<size_t>> search_fm_index(
        const FMIndex& index,
        const vector<string>& queries) {

    vector<vector<size_t>> matches(queries.size());
    for (size_t i = 0; i < queries.size(); i++) {
        dna4_vector pattern(to_dna4_vector(queries[i]));

        auto cursor = index.cursor();
        if (!cursor.extend_right(pattern))
            continue;

        for (auto& [ref_id, pos] : cursor.locate())
            matches[i].push_back(pos);
        sort(matches[i].begin(), matches[i].end());
    }
    return matches;
}

vector<vector<size_t>> solve_with_fm_index(
        const string& reference,
        const vector<string>& queries) {

    fm_index index = build_fm_index(reference);
    vector<vector<size_t>> matches = search_fm_index(index, queries);
    return matches;
}