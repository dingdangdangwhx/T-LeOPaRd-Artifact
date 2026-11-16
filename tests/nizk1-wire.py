#!/usr/bin/env python3
"""Exercise serialized q72 client proofs; benchmark seeds are test fixtures."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import struct
import subprocess
import tempfile


def fields(text):
    return dict(line.split("=", 1) for line in text.splitlines() if "=" in line)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--client", type=Path, default=Path("src/protocol/threshold-client"))
    parser.add_argument("--seeds", default="0,1,2,3,4,5,6,7,12")
    parser.add_argument("--patterns", action="store_true")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    client = str(args.client.resolve())
    args.output.mkdir(parents=True, exist_ok=True)
    tag = "504f5052462d7461672d303030303100"
    cases = [(int(s), "5ac31980277104bed261") for s in args.seeds.split(",")]
    if args.patterns:
        cases += [(12, "00" * 10), (12, "ff" * 10), (12, "aa" * 10)]
    rows = []
    for case, (seed, bits) in enumerate(cases):
        env = dict(os.environ, OMP_NUM_THREADS=str(args.threads), LEOPARD_PUBLIC_THREADS=str(args.threads),
                   LEOPARD_BENCH_SEED_OFFSET=str(seed), LEOPARD_REQUEST_CACHED_PREPROC="0")
        with tempfile.TemporaryDirectory(prefix="nizk1-security-") as tmp:
            root = Path(tmp)
            pp, preq, prestate, req, state = [root / name for name in
                                           ("pp", "preq", "prestate", "req", "state")]

            def run(label, *arguments, accept=True):
                result = subprocess.run([client, *map(str, arguments)], env=env,
                                        capture_output=True, text=True, timeout=180)
                (args.output / f"case{case}-{label}.txt").write_text(result.stdout + result.stderr)
                data = fields(result.stdout)
                accepted = result.returncode == 0 and data.get("ok") == "1"
                if accepted != accept:
                    raise RuntimeError(f"case {case}, {label}: unexpected acceptance={accepted}")
                return data

            run("setup", "setup", pp)
            run("preproc", "preproc-client", pp, 1, preq, prestate)
            request = run("request", "request", pp, prestate, 3, 2, "0,2", tag, bits, req, state)
            verify = run("verify", "verify-request", pp, tag, req)
            original = req.read_bytes()
            changed = bytearray(original)
            budget = struct.unpack_from("<Q", changed, 84)[0]
            struct.pack_into("<Q", changed, 84, budget + 1)
            bad = root / "changed-budget"
            bad.write_bytes(changed)
            run("reject-budget", "verify-request", pp, tag, bad, accept=False)
            changed = bytearray(original)
            changed[92] ^= 1
            bad = root / "changed-ar-seed"
            bad.write_bytes(changed)
            run("reject-seed", "verify-request", pp, tag, bad, accept=False)
            row = {"seed": seed, "input": bits, "request": float(request["time_sec"]),
                   "verify": float(verify["time_sec"]),
                   "rounds": int(request["req_inner_rejection_rounds"]),
                   "proof_bytes": int(request["req_proof_bytes"]),
                   "request_bytes": int(request["req_bytes"]),
                   "request_sha256": hashlib.sha256(original).hexdigest()}
            rows.append(row)
            print(json.dumps(row), flush=True)
    (args.output / "summary.json").write_text(json.dumps(rows, indent=2) + "\n")
    print(f"PASS: {len(rows)} honest requests and {2 * len(rows)} tampering checks")


if __name__ == "__main__":
    main()
