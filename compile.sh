#!/bin/bash
set -e

cmake -S source -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON > /dev/null

cmake --build build -j2