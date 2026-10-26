#!/usr/bin/env bash
set -euo pipefail

echo "=== 1. Verifying POA Logic ==="
python3 scripts/verify_logic.py

if command -v cargo >/dev/null 2>&1; then
    echo -e "\n=== 2. Running Rust POA CLI ==="
    cargo run --release -- -i data/sample_poa.fa -c consensus.fa -g graph.gfa -d graph.dot
    echo "Generated consensus:"
    cat consensus.fa
    rm -f consensus.fa graph.gfa graph.dot
else
    echo -e "\n[Note: cargo toolchain not yet installed; algorithmic verification ran via Python suite]"
fi

echo -e "\nExample run completed successfully!"
