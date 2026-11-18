"""Search/replay the first accepting round of one public synthetic epoch request."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
TAG = "504f5052462d7461672d303030303100"
EPOCH_SEED = "6cb0da21981144bd9a457b80469d539e"
ONE_ROUND_OFFSET = 9


def options():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profile", choices=("partial72", "plain72"))
    parser.add_argument("--start", type=int, default=ONE_ROUND_OFFSET)
    parser.add_argument("--limit", type=int, default=1)
    parser.add_argument("--runs", type=int, default=3)
    parser.add_argument("--lock-fd", type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if not 1 <= args.runs <= 1000:
        parser.error("--runs must be in [1,1000]")
    if args.start < 0 or args.limit < 1 or args.start + args.limit > 2**32:
        parser.error("candidate range is invalid")
    return args


def main(args):
    os.environ.setdefault("OMP_NUM_THREADS", "8")
    name = "request-round-" + args.profile
    subprocess.run(["make", "-C", ROOT / "src/protocol", "PROFILE=" + args.profile,
                    "threshold-client", "../../build/" + name], check=True, stdout=subprocess.DEVNULL)
    with tempfile.TemporaryDirectory() as temp:
        work = Path(temp)
        def invoke(offset):
            start = time.perf_counter()
            result = subprocess.run([ROOT / "build" / name, str(offset)], cwd=work,
                                    capture_output=True, text=True, check=True)
            elapsed = time.perf_counter() - start
            v = dict(line.split("=", 1) for line in result.stdout.splitlines()
                     if "=" in line and line.partition("=")[0].isidentifier())
            if v["ok"] != "1" or v["req_cached_preprocessing"] != "0" or v["req_inner_batch_rounds"] != "1":
                raise RuntimeError("unexpected fixture mode")
            raw = (work / "fixture.req").read_bytes()
            n = int(v["req_cr_wire_bytes"]) + int(v["req_public_cx_bytes"])
            if len(raw) != int(v["req_bytes"]) or len(raw) < 108+n or not n:
                raise RuntimeError("unexpected q72 request format")
            cr = int(v["req_cr_wire_bytes"])
            breakdown = {"epoch ID": 8, "Cx": int(v["req_public_cx_bytes"]),
                         "commitment to R": cr, "other public": int(v["req_public_other_bytes"]) - cr,
                         "outer proof": int(v["req_proof_bytes"])}
            if sum(breakdown.values()) != 8+len(raw):
                raise RuntimeError("request byte accounting mismatch")
            return dict(byte_breakdown=breakdown, rounds=int(v["req_inner_rejection_rounds"]), wall_s=elapsed,
                        inner_s=float(v["req_inner_prove_sec"]), request_core_s=float(v["time_sec"]),
                        request_bytes=8+len(raw), request_hash=hashlib.sha256(raw).hexdigest(),
                        statement_hash=hashlib.sha256(raw[108:108+n]).hexdigest())

        started = time.perf_counter()
        statement = None
        history = []
        for offset in range(args.start, args.start + args.limit):
            row = invoke(offset)
            if statement is not None and statement != row["statement_hash"]:
                raise RuntimeError("search changed the mask commitment or blinded matrix")
            statement = row["statement_hash"]
            history.append(row["rounds"])
            if row["rounds"] == 1:
                break
        else:
            raise SystemExit(f"no first-round acceptance in {len(history)} trials: {history}")
        search_s = time.perf_counter() - started
        replay = []
        for _ in range(args.runs):
            r = invoke(offset)
            if r["rounds"] != 1 or r["request_hash"] != row["request_hash"]:
                raise RuntimeError("selected request did not reproduce")
            subprocess.run([ROOT / "src/protocol/threshold-client", "epoch-verify",
                            "fixture.pp", TAG, "fixture.req", EPOCH_SEED], cwd=work,
                           stdout=subprocess.DEVNULL, check=True)
            replay.append(r)
        print(json.dumps(dict(profile=args.profile, fixture_only=True, offset=offset,
                              search_seconds=search_s, search_rounds=history,
                              ordinary_verifications=len(replay), replay=replay,
                              means={k: statistics.mean(r[k] for r in replay) for k in
                                     ("wall_s", "inner_s", "request_core_s", "request_bytes")}), indent=2))


if __name__ == "__main__":
    args = options()
    path = ROOT / ".artifact.lock"
    with (path.open("a") if args.lock_fd is None else os.fdopen(os.dup(args.lock_fd), "a")) as lock:
        if not os.path.samestat(os.fstat(lock.fileno()), path.stat()):
            raise SystemExit("unexpected artifact lock")
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise SystemExit("another artifact benchmark is running")
        main(args)
