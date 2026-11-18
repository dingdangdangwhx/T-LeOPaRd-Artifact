# T-LeOPaRd Artifact

Artifact for **T-LeOPaRd: Efficient Threshold Oblivious PRFs from Lattices**.
It contains parameter-selection scripts and C implementations with and without
NIZKs, using 2-out-of-3 replicated sharing at correctness level kappa = 32.

The anonymous [full version of the paper](paper/full-version.pdf) is included.

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
two replies in parallel within the total thread budget (two groups of four
by default). PreProc measures one physical server publishing both component
vectors for one epoch; it needs no client input. Request samples a fresh mask
and sends its commitment online. Both settings include the public rounding
offset and the secret-seeded direct evaluation noise. KB means 1024 bytes.
For the two proof-enabled profiles the main Request row uses a fixed
synthetic client request (seed offset 9, one rejection round); its time
includes the client process, not seed search or deployment-state I/O.
The report also gives the fresh-coin end-to-end Request time and rejection
counts. Other proof-enabled online timings include CLI and state handling;
semi-honest times are in-process.
Communication is serialized-file size with NIZKs and packed-field accounting
without NIZKs, excluding transport overhead.

To repeat only the default one-round client measurement:

```sh
python3 tests/request-round.py partial72
python3 tests/request-round.py plain72
```

The three published client timings per profile are in [the reference measurements](paper/client-request.csv).
Both profiles replay offset 9; add `--start 0 --limit 64` to repeat the search.
The fixture fixes the public inputs, mask and opening, then verifies the
byte-identical replays with the normal verifier. Its request timer covers
the client CLI rather than the deployment driver. These are conditional
one-round timings, not fresh-coin averages. Only the fixture binary links
test coins; it does not accept user inputs. `run.sh` measures it alongside
the ordinary epoch protocol; it never sets the live prover to a fixed seed.

## Parameters

SageMath is needed only to recompute parameters; benchmark builds use the
included proof-parameter headers.

```sh
make parameters
make proof-parameters
make client-parameters
```

`make parameters` writes the kappa = 32 protocol candidates to
`results/parameters.json`. To reproduce all parameter-table rows:

```sh
sage -python parameters/compute.py --kappa 16 32 64 --output results/parameters.json
```

The q72 client uses a binary-opening commitment (modulus 65437, rank 14),
masked Bin1 range checks, and a verifier-enforced outer squared budget
1653435896617783. The conditional inner MSIS screen uses extraction slack 3,
rank 24 and randomness dimension 74. `make client-parameters` checks the
commitment, lifting/no-wrap inequalities, and ideal challenge filter.
These screens are not a certificate for the composed proof system: the complete
recursive extraction, simulation, and concrete Fiat--Shamir loss analysis remain
separate proof obligations.

The server proof retains uniform key shares in its unbounded message block.
Commitment errors are affine residuals rather than separate witnesses, with
the same joint opening L2 bound. Generated headers include the corresponding
LNP22 parameters; sparse affine reconstruction and verification are optimized.

The epoch driver uses fresh client and server coins without selecting rejection
rounds. It enforces one publication per server and epoch, and limits responses
to 1024 (Partial) or 4096 (Plain) uses per server in the proof-enabled profiles.
The proof-free profiles use one use per server per epoch. The warmup checks
all three quorums; negative tests additionally exercise the proof-enabled
epoch boundary and rejection paths. This is research code, not a production service.
The revised q72 proof/preprocessing format is not compatible with earlier files;
regenerate them together. The semi-honest parameter sets are unchanged.

## Layout

- `parameters/`: protocol parameter selection and NIZK2 parameter generation.
- `src/protocol/`: with-NIZK tools and the semi-honest benchmark.
- `src/lazer/`, `third_party/`: proof and arithmetic dependencies.

Run `make clean` to remove build products. Dependency licenses are retained.
`OMP_NUM_THREADS=8 make check` additionally tests both proof backends, changed
norm budgets, ARP boundaries, packed constraints, cached preprocessing, and
single/aggregate server commitment extraction for both with-NIZK profiles.
