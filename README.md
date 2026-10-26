> [!NOTE]
> **Archive Notice**: This project is an educational and reference implementation of Partial Order Alignment (Sequence-to-DAG) in Rust. It is finished and archived; no active maintenance or feature development is planned.

# poa: Partial Order Alignment & Variation Graph Engine

[![Rust](https://img.shields.io/badge/rust-2021-orange.svg)](https://www.rust-lang.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Zero Dependencies](https://img.shields.io/badge/dependencies-none-brightgreen.svg)]()

`poa` is a high-performance Partial Order Alignment (Sequence-to-DAG) and variation graph engine written in pure Rust. It generalizes dynamic programming alignment from linear sequences to Directed Acyclic Graphs (DAGs) (Lee et al., 2002), facilitating multiple sequence alignment (MSA), consensus calling, bubble detection, and variation graph generation.

---

## Key Features

1. **CSR-Style Flat Arena Graph (`PoaGraph`)**:
   - Represents variation graphs using flat contiguous arrays for bases, in-edges, and out-edges, minimizing cache misses during topological traversals.
2. **Circular Sibling Ring Allele Fusion**:
   - Mismatched bases at the same homologous column are linked into a circular ring (`aligned_sibling`), avoiding duplicate divergent paths and capturing clean heterozygous bubbles.
3. **Topological Order Affine Gap DP**:
   - Supports global and local sequence-to-DAG alignment with match (+2), mismatch (-4), gap open (-6), and gap extend (-2) scoring matrices.
4. **Heaviest-Bundling Consensus Sequence Extraction**:
   - Identifies the highest-weight path through the variation graph using edge-support weights, producing accurate consensus sequences even from noisy long reads.
5. **GFA v1 & Graphviz DOT Export**:
   - Exports graphs to standard **GFA v1** format (with segment coverage `ic:i:` and edge weights `rc:i:`) and colored **Graphviz DOT** files for visualization.

---

## Architectural Overview

```
Input Multi-FASTA Reads
          |
    [ Init PoaGraph with Sequence 1 ]
          |
    +-----> For each subsequent sequence Q:
    |       |
    |       [ Sequence-to-DAG Topological DP ]
    |       |
    |       [ Optimal Path Backtrace ]
    |       |
    |       [ Graph Node & Sibling Bubble Fusion ]
    |       |
    +-------+
          |
    [ Heaviest Bundling DP ]
          |
    +-----+-----+-----+
    |           |     |
Consensus      GFA   Graphviz
 FASTA         v1     DOT
```

---

## Build & Test

### Building via Cargo
```bash
cargo build --release
cargo test --release
```

### Running via Makefile
```bash
make test       # Runs unit tests or standalone verification suite
make all        # Builds release binary
```

---

## CLI Usage

```bash
# Basic consensus and GFA export
cargo run --release -- -i data/sample_poa.fa -c consensus.fa -g graph.gfa -d graph.dot

# Multiple Sequence Alignment export
cargo run --release -- -i data/sample_poa.fa -m msa.fa

# Local alignment mode with custom penalties
cargo run --release -- -i reads.fa -c consensus.fa --mode local --match 3 --mismatch -3 --gap-open -5
```

---

## Verification & Benchmarks

Run automated verification:
```bash
./scripts/test_example.sh
./scripts/benchmark_compare.sh
```

### Scaling Across Read Lengths (5 Homologous Reads):
| Read Length | Variation Graph Nodes | Consensus Latency |
| :--- | :--- | :--- |
| **100 bp** | 105 nodes | **38.39 ms** |
| **200 bp** | 210 nodes | **72.48 ms** |
| **500 bp** | 527 nodes | **443.91 ms** |

---

## License

MIT License. Copyright (c) 2026 Andreas Sagen.
