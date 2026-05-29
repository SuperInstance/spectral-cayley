# spectral-cayley

**Spectral analysis on Cayley graphs — group structure made visible through eigenvalues.**

Pure Rust, zero dependencies. Computes eigenvalues and spectral properties of Cayley graphs generated from group elements and generating sets.

## What This Gives You

- **Cayley graph construction** — from any group table and generating set
- **Spectral decomposition** — eigenvalues via Jacobi iteration on the graph Laplacian
- **Conservation ratio** — CR = λ₂/λ_max for Cayley graphs
- **Zero dependencies** — all computation from scratch

## Quick Start

```rust
use spectral_cayley::CayleyGraph;

// Construct from group table and generators
let group_table = vec![/* your group's multiplication table */];
let generators = vec![1, 3]; // group elements that generate the graph
let cg = CayleyGraph::new(group_table, generators);

// Eigenvalues of the Laplacian
let eigs = cg.eigenvalues();
let cr = cg.conservation_ratio();
```

## How It Fits

Part of the SuperInstance spectral ecosystem. Cayley graphs are a natural testing ground for conservation ratio theory because their symmetry gives clean, interpretable eigenvalue spectra.

- **[spectral-graph-core](https://github.com/SuperInstance/spectral-graph-core)** — General spectral graph theory
- **[spectral-cayley](https://github.com/SuperInstance/spectral-cayley)** — Cayley graph specialization (this repo)

## Testing

```bash
cargo test
```

## Installation

```toml
[dependencies]
spectral-cayley = { git = "https://github.com/SuperInstance/spectral-cayley" }
```

## License

MIT

Part of the [SuperInstance](https://github.com/SuperInstance) ecosystem.
