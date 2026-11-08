#include <chrono>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include "fm_index/fm_index.hpp"

#include "single_seed/single_seed_index.hpp"
#include "single_seed/single_seed_sorted.hpp"
#include "single_seed/single_seed_vector.hpp"

#include "paired_seed/paired_seed_index.hpp"
#include "paired_seed/paired_seed_sorted.hpp"
#include "paired_seed/paired_seed_vector.hpp"

using namespace std;

size_t get_memory_kb() {
    ifstream file("/proc/self/status");
    string line;

    while (getline(file, line)) {
        if (line.starts_with("VmRSS:")) {
            istringstream stream(line.substr(6));
            size_t kb;
            stream >> kb;
            return kb;
        }
    }

    throw runtime_error("Cannot read memory usage");
}

template<typename Build, typename Search>
auto benchmark(const string& name, Build build, Search search) {
    using Clock = chrono::steady_clock;
    size_t memory_before = get_memory_kb();
    auto start = Clock::now();
    auto index = build();
    auto built = Clock::now();
    size_t memory_after = get_memory_kb();
    auto matches = search(index);
    auto finished = Clock::now();
    double build_ms = chrono::duration<double, milli>(built - start).count();
    double search_ms = chrono::duration<double, milli>(finished - built).count();
    int64_t index_memory_kb = static_cast<int64_t>(memory_after) - static_cast<int64_t>(memory_before);

    cerr << fixed << setprecision(3);
    cerr << name << ":\n";
    cerr << "  Construction: " << build_ms << " ms\n";
    cerr << "  Searching:    " << search_ms << " ms\n";
    cerr << "  Total:        " << build_ms + search_ms << " ms\n";
    cerr << "  Est. index RAM: " << index_memory_kb / 1024.0 << " MiB\n";
    return matches;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: ./build/main <fmi|ssi|psi|sss|pss|ssv|psv>\n";
        return 1;
    }
    string algo = argv[1];

    int m, q, k;
    string reference;
    cin >> m >> q >> k >> reference;

    vector<string> queries(q);
    for (string& p : queries) cin >> p;

    if (algo == "fmi")
        auto fmi = benchmark("FM-index",
            [&]() { return build_fm_index(reference); },
            [&](auto& index) { return search_fm_index(index, queries); }
        );
    else if (algo == "ssi")
        auto ssi = benchmark("Single-seed-hashmap",
            [&]() { return build_single_seed_index(reference, k); },
            [&](auto& index) { return search_single_seed_index(index, reference, queries); }
        );
    else if (algo == "psi")
        auto psi = benchmark("Paired-seed-hashmap",
            [&]() { return build_paired_seed_index(reference, m, k); },
            [&](auto& index) { return search_paired_seed_index(index, reference, queries); }
        );
    else if (algo == "sss")
        auto sss = benchmark("Single-seed-sorted",
            [&]() { return build_single_seed_sorted_index(reference, k); },
            [&](auto& index) { return search_single_seed_sorted_index(index, reference, queries); }
        );
    else if (algo == "pss")
        auto pss = benchmark("Paired-seed-sorted",
            [&]() { return build_paired_seed_sorted_index(reference, m, k); },
            [&](auto& index) { return search_paired_seed_sorted_index(index, reference, queries); }
        );
    else if (algo == "ssv")
        auto ssv = benchmark("Single-seed-static",
            [&]() { return build_single_seed_vector_index(reference, k); },
            [&](auto& index) { return search_single_seed_vector_index(index, reference, queries); }
        );
    else if (algo == "psv")
        auto psv = benchmark("Paired-seed-static",
            [&]() { return build_paired_seed_vector_index(reference, m, k); },
            [&](auto& index) { return search_paired_seed_vector_index(index, reference, queries); }
        );
    else
        cout << "UNKOWN ALGO" << endl;
}