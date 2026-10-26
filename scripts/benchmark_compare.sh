#!/usr/bin/env bash
set -euo pipefail

echo "=== Benchmark: POA Sequence-to-DAG Scaling ==="

python3 -c '
import time, random
from verify_logic import SimplePoa

bases = ["A", "C", "G", "T"]
for length in [100, 200, 500]:
    ancestor = "".join(random.choice(bases) for _ in range(length))
    reads = []
    for _ in range(5):
        r = list(ancestor)
        for i in range(0, len(r), 30):
            r[i] = random.choice([b for b in bases if b != r[i]])
        reads.append("".join(r))

    t0 = time.perf_counter()
    p = SimplePoa()
    p.init_seq(reads[0])
    for r in reads[1:]:
        p.align_and_fuse(r)
    cons = p.consensus()
    elapsed = time.perf_counter() - t0

    print(f"5 reads of {length:4d} bp: Graph Nodes = {len(p.bases):4d}, Consensus Time = {elapsed*1000:6.2f} ms")
'

echo "Benchmark completed."
