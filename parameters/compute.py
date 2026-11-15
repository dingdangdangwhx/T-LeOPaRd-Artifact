"""Byte-aligned protocol parameter candidates; run with sage -python.

This screens honest correctness and, in malicious mode, the existing exact-L2
server relation. It does not generate proof-system parameters or certify the
composed client proof. Client ranks retain the manuscript's ternary-reference
MLWE model; validation for the implemented binary client remains separate.
"""

import argparse
import contextlib
import io
import json
import os
import sys
from functools import lru_cache
from pathlib import Path

from sage.all import Integer, RealField, ceil, log, pi, sqrt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
os.environ.setdefault("LWE_ESTIMATOR_DIR", str(ROOT / "third_party"))

from t_leopard.hardness import MLWEDimensions, TARGET_RHF, lattice_estimator
from t_leopard.proof_bounds import noise_budgets
from t_leopard.threshold import (
    ReplicatedSharing, correctness_modulus_bound, threshold_widths,
)
from t_leopard.widths import max_queries_per_input

RF = RealField(256)
SHARING = ReplicatedSharing(3, 2)
DIMENSIONS = MLWEDimensions()
OFFSETS = {24: 3, 32: 99, 40: 195, 48: 59, 56: 27, 64: 59, 72: 107}


@lru_cache(None)
def modulus(bits):
    offset = Integer(OFFSETS.get(bits, 3))
    while not (Integer(2) ** bits - offset).is_prime():
        offset += 8
    return Integer(2) ** bits - offset


def at_bits(mode, kappa, sessions, bits, security_model):
    q, d, h, p = modulus(bits), Integer(64), Integer(1), Integer(4)
    ell = Integer(DIMENSIONS.rank(TARGET_RHF[128], d, Integer(bits), Integer(1), True))
    m = Integer(1)
    for _ in range(32):
        s0, minimum_s, minimum_s1, _ = threshold_widths(
            SHARING, max_queries_per_input(mode, Integer(kappa)), d, m + ell, h,
            sessions=Integer(sessions),
        )
        # The reference estimator rounds a computed Gaussian stddev upwards.
        next_m = Integer(DIMENSIONS.rank(
            TARGET_RHF[128], d, Integer(bits), ceil(RF(s0) / sqrt(2 * RF(pi))), False,
        ))
        if next_m == m:
            break
        m = next_m
    else:
        raise RuntimeError("rank/width iteration did not converge")
    unit = RF("1.55") * sqrt(2 * RF(pi))
    es_log = int(ceil(log(RF(minimum_s) / unit, 2)))
    s = unit * Integer(2) ** es_log
    esp_log = int(ceil(log(RF(minimum_s1) * s / RF(minimum_s) / unit, 2)))
    s1 = unit * Integer(2) ** esp_log
    sigma, sigma1 = s / sqrt(2 * RF(pi)), s1 / sqrt(2 * RF(pi))
    rms = sqrt((SHARING.width + 1) * sigma1**2 + SHARING.width * (m + ell) * d * sigma**2)
    rms_q = p * h * d * Integer(2) ** kappa * rms
    conservative_q = correctness_modulus_bound(
        Integer(kappa), p, d, m, ell, h, s, s1, Integer(1), SHARING,
    )
    bounds = noise_budgets(s, s1, m + ell, h, d)
    bs = 2 * bounds.preprocessing
    accepted = SHARING.width * (sqrt(RF((m + ell) * d)) * bs + bounds.online)
    uniqueness_q = SHARING.quorums * p * h * d * Integer(2) ** kappa * (2 * accepted + 1)
    required_q = max(conservative_q, uniqueness_q) if security_model == "malicious" else rms_q
    return {
        "model": security_model, "mode": mode, "kappa": kappa, "sessions": sessions,
        "logq": bits, "q": int(q), "q_offset": int(Integer(2) ** bits - q),
        "m": int(m), "ell": int(ell), "p": int(p), "d": int(d), "h": int(h),
        "s0": float(s0), "s": float(s), "s1": float(s1),
        "sigma_es": float(sigma), "sigma_online": float(sigma1),
        "es_sampler_log": es_log, "online_sampler_log": esp_log,
        "B_s": int(bs), "B_1": int(bounds.online),
        "digit_base": int(bounds.digit_base),
        "digit_low": int(bounds.digit_low), "digit_high": int(bounds.digit_high),
        "correctness_required_logq": float(log(conservative_q, 2)),
        "rms_required_logq": float(log(rms_q, 2)),
        "uniqueness_required_logq": float(log(uniqueness_q, 2)),
        "required_logq": float(log(required_q, 2)),
        "log2_rms_failure": float(log(rms_q / q, 2) - kappa),
        "log2_uniqueness": float(log(uniqueness_q / q, 2) - kappa),
        "passes": bool(q >= required_q),
    }


@lru_cache(None)
def commitment_rank(bits, m, initial_rank):
    est = lattice_estimator()
    q = modulus(bits)
    for rank in range(initial_rank, initial_rank + 64):
        opening = sqrt(RF(4 * (2 * rank + 2 * m) * 64))
        error = sqrt(RF(2 * rank * 64 + 1)) * opening
        scale = q.isqrt()
        extraction = 2 * error < scale and 2 * (scale + 1) * error < q
        sis_delta = RF(2) ** (log(2 * opening, 2)**2 / (4 * rank * 64 * log(q, 2)))
        costs = []
        for model in (est.BKZ.enum, est.BKZ.sieve):
            with contextlib.redirect_stdout(io.StringIO()):
                result = est.estimate_lwe(
                    rank * 64, est.alphaf(est.sigmaf(RF("0.5")), q), q,
                    secret_distribution=(0, 1), reduction_cost_model=model,
                )
            costs.extend(float(log(value["rop"], 2)) for value in result.values() if "delta_0" in value)
        if extraction and sis_delta <= RF("1.0044") and min(costs) >= 128:
            return rank, float(error), float(sis_delta), min(costs)
    raise RuntimeError("commitment rank search failed")


def search(mode, kappa, sessions, security_model):
    bits = max(24, 8 * int(ceil(RF(kappa + 16) / 8)))
    history = []
    for _ in range(24):
        row = at_bits(mode, kappa, sessions, bits, security_model)
        history.append({"logq": bits, "required_logq": row["required_logq"]})
        next_bits = max(24, 8 * int(ceil(RF(row["required_logq"]) / 8)))
        if next_bits == bits and row["passes"]:
            break
        bits = next_bits
    else:
        raise RuntimeError("modulus iteration failed")
    row["previous_byte_aligned_passes"] = at_bits(mode, kappa, sessions, bits - 8, security_model)["passes"] if bits > 24 else False
    assert not row["previous_byte_aligned_passes"]
    rank, error, sis_delta, lwe_cost = commitment_rank(bits, row["m"], row["ell"])
    row.update({"n_s": rank, "commitment_extraction_error": error,
                "commitment_sis_delta": sis_delta, "commitment_min_legacy_lwe_bits": lwe_cost,
                "search_history": history})
    rows = rank + 2 * row["m"]
    row["sizes_KiB"] = {
        "pk_payload_with_seed_hash": (3 * rows * 64 * bits // 8 + 32) / 1024,
        "projected_component_files_without_seed_hash": (3 * rows * 64 * max(8, bits // 8) + 60) / 1024,
        "server_preprocessing_per_slot": 2 * (row["m"] + row["ell"]) * 64 * bits / 8 / 1024,
        "server_online_per_bundle": 64 * bits / 8 / 1024,
        "client_online_cx": row["m"] * 64 * bits / 8 / 1024,
    }
    row["proof_parameters_regenerated"] = False
    row["backend_modulus_currently_supported"] = bits in OFFSETS
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--kappa", type=int, choices=(16, 32, 64), nargs="+", default=(16, 32, 64))
    parser.add_argument("--models", choices=("malicious", "semi-honest-one-time", "semi-honest-matched-reuse"),
                        nargs="+", default=("malicious", "semi-honest-one-time"))
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for model in args.models:
        for mode in ("POPRF", "OPRF"):
            for kappa in args.kappa:
                sessions = 1 if model == "semi-honest-one-time" else (
                    2**10 if mode == "POPRF" or kappa == 16 else 2**12 if kappa == 32 else 2**32
                )
                print(f"screening {model} {mode} kappa={kappa} K={sessions}", file=sys.stderr, flush=True)
                rows.append(search(mode, kappa, sessions, model))
                # Preserve completed cases if a later estimator call is slow.
                args.output.write_text(json.dumps({
                    "scope": "Conditional protocol candidates, not new measured profiles or composed-proof certification. Reference client ranks require separate binary-client validation. Semi-honest rows require uniform-offset RMS correctness analysis. Sizes retain commitment-based offline layout, include one 32-byte seed hash in PK, and omit serialization headers and NIZKs.",
                    "rows": rows,
                }, indent=2) + "\n")


if __name__ == "__main__":
    main()
