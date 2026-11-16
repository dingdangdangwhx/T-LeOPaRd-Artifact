#!/usr/bin/env python3
"""Check that cached preprocessing removes only duplicate public bytes."""
import hashlib
import json
import os
from pathlib import Path
import struct
import subprocess
import tempfile

HERE = Path(__file__).resolve().parents[1]
CLIENT = HERE / "src/protocol/threshold-client"
OUT = HERE / "results/cache-test"
OUT.mkdir(parents=True, exist_ok=True)
ENV = dict(os.environ, OMP_NUM_THREADS="8", LEOPARD_PUBLIC_THREADS="8",
           LEOPARD_BENCH_SEED_OFFSET="2", LEOPARD_REQUEST_CACHED_PREPROC="0")
TAG = "504f5052462d7461672d303030303100"
checks = 0


def run(label, *args, accept=True):
    global checks
    result = subprocess.run([str(CLIENT), *map(str, args)], env=ENV,
                            capture_output=True, text=True, timeout=180)
    (OUT / f"{label}.txt").write_text(result.stdout + result.stderr)
    fields = dict(line.split("=", 1) for line in result.stdout.splitlines() if "=" in line)
    assert result.returncode >= 0, f"signal/crash in {label}"
    assert (result.returncode == 0 and fields.get("ok") == "1") == accept, label
    checks += 1
    return fields


with tempfile.TemporaryDirectory(prefix="cached-preprocessing-") as tmp:
    root = Path(tmp)
    pp, preq, prestate, stateless, cached, state = [root / name for name in
                                                 ("pp", "preq", "prestate", "stateless", "cached", "state")]
    run("setup", "setup", pp)
    run("preprocess", "preproc-client", pp, 1, preq, prestate)
    args = ("request", pp, prestate, 3, 2, "0,2", TAG, "5ac31980277104bed261")
    old = run("request-stateless", *args, stateless, state)
    run("verify-stateless", "verify-request", pp, TAG, stateless)
    run("verify-stateless-with-cache", "verify-request", pp, TAG, stateless, preq)
    ENV["LEOPARD_REQUEST_CACHED_PREPROC"] = "1"
    new = run("request-cached", *args, cached, state)
    run("verify-cached", "verify-request", pp, TAG, cached, preq)
    run("reject-no-cache", "verify-request", pp, TAG, cached, accept=False)
    before, after = stateless.read_bytes(), cached.read_bytes()
    removed = int(new["req_cr_omitted_bytes"])
    assert len(before) - len(after) == removed
    assert before[4:108] == after[4:108]
    assert before[108 + removed:] == after[108:]
    assert old["req_proof_bytes"] == new["req_proof_bytes"]
    mutations = {"slot": (28, "<Q", 1), "cr-rows": (36, "<Q", 102),
                 "budget": (84, "<Q", struct.unpack_from("<Q", after, 84)[0] + 1)}
    for name, (off, fmt, value) in mutations.items():
        bad = bytearray(after)
        struct.pack_into(fmt, bad, off, value)
        path = root / name
        path.write_bytes(bad)
        run("reject-" + name, "verify-request", pp, TAG, path, preq, accept=False)
    for name, off in (("seed", 92), ("proof", -1)):
        bad = bytearray(after)
        bad[off] ^= 1
        path = root / name
        path.write_bytes(bad)
        run("reject-" + name, "verify-request", pp, TAG, path, preq, accept=False)
    badcache = bytearray(preq.read_bytes())
    badcache[-2] ^= 1
    wrong = root / "wrong-cache"
    wrong.write_bytes(badcache)
    run("reject-wrong-cache", "verify-request", pp, TAG, cached, wrong, accept=False)
    run("reject-wrong-cache-stateless", "verify-request", pp, TAG, stateless, wrong, accept=False)
    summary = dict(checks=checks, request_stateless_bytes=len(before),
                   request_cached_bytes=len(after), omitted_bytes=removed,
                   proof_bytes=int(new["req_proof_bytes"]),
                   stateless_sha256=hashlib.sha256(before).hexdigest(),
                   identical_transcript_and_proof=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print("Cached request regression PASS")
