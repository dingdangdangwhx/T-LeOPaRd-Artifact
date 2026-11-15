# T-LeOPaRd Artifact

Artifact for **T-LeOPaRd: Efficient Threshold Oblivious PRFs from Lattices**.
It contains parameter-selection scripts and C implementations with and without
NIZKs, using 2-out-of-3 replicated sharing at correctness level kappa = 32.

| Setting | Partial modulus | Plain modulus |
|---|---|---|
| With NIZKs | 2^72 - 107 | 2^72 - 107 |
| Semi-honest | 2^56 - 27 | 2^64 - 59 |

## Requirements

Linux x86-64, AVX2, AES-NI, AVX-512F/BW/DQ/VL, GCC/G++, Make, CMake,
Python 3, GMP and MPFR. On Ubuntu:

```sh
sudo apt install build-essential cmake python3 libgmp-dev libmpfr-dev util-linux procps
```

The proof backends and HEXL dependencies are included. Builds use
`-march=native`; compile on the machine used for measurement.

## Run

```sh
bash run.sh
```

This builds from source and runs all four configurations, with one warmup and
three measured runs per configuration. Correctness is checked before reporting
the mean runtime and communication of each algorithm, plus Request/BlindEval
breakdowns. Raw logs, hardware information, `summary.csv`, `summary.json`, and
`summary.md` are written to `results/<timestamp>/`.

```sh
RUNS=5 OMP_NUM_THREADS=8 bash run.sh
```

BlindEval is reported per server, with each server verifying the request
independently; server proofs are generated sequentially. Finalize verifies the
two replies in parallel. Request communication counts one message. Preprocessing
covers ten instances, and server preprocessing covers all three servers without
reuse amortization. KB means 1024 bytes. The with-NIZK Eval timer measures PRF
arithmetic; the semi-honest implementation includes secret-seeded evaluation noise.
Communication is serialized-file size with NIZKs and packed-field accounting
without NIZKs, excluding transport overhead.

## Parameters

SageMath is needed only to recompute parameters; benchmark builds use the
included proof-parameter headers.

```sh
make parameters
make proof-parameters
```

`make parameters` writes the kappa = 32 protocol candidates to
`results/parameters.json`. To reproduce all parameter-table rows:

```sh
sage -python parameters/compute.py --kappa 16 32 64 --output results/parameters.json
```

These are protocol parameter screens, not a certificate for the composed proof
system. The benchmark uses a fixed reproducible client fixture and fresh server
coins. This is research code, not a production cryptographic service.

## Layout

- `parameters/`: protocol parameter selection and NIZK2 parameter generation.
- `src/protocol/`: with-NIZK tools and the semi-honest benchmark.
- `src/lazer/`, `third_party/`: proof and arithmetic dependencies.

Run `make clean` to remove build products. Dependency licenses are retained.
`OMP_NUM_THREADS=8 make check` additionally tests server-proof rejection and
single/aggregate commitment extraction for both with-NIZK profiles.
