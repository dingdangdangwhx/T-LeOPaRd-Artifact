"""Summarize the four artifact benchmarks using only the Python standard library."""
import csv
import json
import math
import statistics
import sys
from pathlib import Path

ALGORITHMS = ("Setup", "KeyGen", "PreProc",
              "Request", "BlindEval", "Finalize", "Eval")
PROFILES = {
    "partial72": ("With NIZKs", "Partial", 72, 42, 45),
    "plain72": ("With NIZKs", "Plain", 72, 42, 45),
    "semi_partial56": ("Semi-honest", "Partial", 56, 32, 35),
    "semi_plain64": ("Semi-honest", "Plain", 64, 37, 40),
}


def malicious_run(path):
    row = json.loads(path.read_text())
    if row["correctness"] != "PASS" or set(row["seconds"]) != set(ALGORITHMS) or set(row["bytes"]) != set(ALGORITHMS):
        raise ValueError(f"incomplete measurement: {path.name}")
    det = row["detail"]
    req = det["request"]
    row["request_seconds"] = {key: float(req["req_" + key + "_sec"]) for key in (
        "protocol_transcript", "inner_instance", "inner_prepare", "inner_prove",
        "outer_witness", "outer_statement", "outer_prove", "serialize")}
    row["blind_seconds"] = {key: statistics.mean(det[key]) for key in ("verify_request", "nizk2_prove")}
    cr = int(req["req_cr_wire_bytes"])
    row["request_bytes"] = {"epoch ID": 8, "Cx": int(req["req_public_cx_bytes"]),
                            "commitment to R": cr, "other public": int(req["req_public_other_bytes"]) - cr,
                            "outer proof": int(req["req_proof_bytes"])}
    row["reply_bytes"] = {"public": det["reply_public"], "proof": det["reply_proofs"]}
    for key, algorithm in (("request_bytes", "Request"), ("reply_bytes", "BlindEval")):
        if sum(row[key].values()) != row["bytes"][algorithm]:
            raise ValueError(f"byte accounting mismatch: {path.name}, {algorithm}")
    return row


def semi_runs(path, expected):
    rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    for row in rows:
        if row["correctness"] != "PASS" or row["quorums"] != 3 or row["proofs_executed"] != 0:
            raise ValueError(f"proof-free correctness failed: {path.name}")
        if (row["logq"], row["m"], row["ell"]) != expected[2:]:
            raise ValueError(f"unexpected profile: {path.name}")
        if len(row["seconds"]) != len(ALGORITHMS) or len(row["bytes"]) != len(ALGORITHMS):
            raise ValueError("incomplete proof-free measurements")
        if any(not math.isfinite(x) or x < 0 for x in row["seconds"] + row["bytes"]):
            raise ValueError("invalid proof-free measurement")
    results = []
    for row in rows:
        message_bytes = dict(zip(ALGORITHMS, row["bytes"]))
        results.append({"correctness": "PASS", "seconds": dict(zip(ALGORITHMS, row["seconds"])),
                        "bytes": message_bytes,
                        "request_seconds": dict(zip(("fast core", "selector matrices", "load mask", "chain", "Cx"),
                                                     row["request_breakdown"]))})
    return results


def aggregate(rows):
    means = {}
    for group in ("seconds", "bytes", "request_seconds", "blind_seconds", "request_bytes", "reply_bytes"):
        if group in rows[0]:
            means[group] = {key: statistics.mean(row[group][key] for row in rows) for key in rows[0][group]}
    return means


def add_fixture(profile, data, count):
    rows = data["replay"]
    if (data["profile"] != profile["id"] or data["fixture_only"] is not True
            or data["offset"] != 9 or data["ordinary_verifications"] != count
            or not rows or len(rows) != count or len({r["request_hash"] for r in rows}) != 1):
        raise ValueError("invalid one-round fixture record")
    for r in rows:
        if (r["rounds"] != 1 or sum(r["byte_breakdown"].values()) != r["request_bytes"]
                or not 0 < r["inner_s"] <= r["wall_s"] < math.inf):
            raise ValueError("invalid one-round fixture measurement")
    stats = profile["mean"]
    profile["fresh_request"] = {
        "seconds": stats["seconds"]["Request"], "bytes": stats["bytes"]["Request"],
        "request_seconds": stats["request_seconds"], "request_bytes": stats["request_bytes"],
        "rounds": [int(r["detail"]["request"]["req_inner_rejection_rounds"]) for r in profile["runs"]]}
    profile["one_round_fixture"] = data
    stats["seconds"]["Request"] = statistics.mean(r["wall_s"] for r in rows)
    stats["bytes"]["Request"] = rows[0]["request_bytes"]
    stats["request_seconds"] = {"inner prove": statistics.mean(r["inner_s"] for r in rows)}
    stats["request_bytes"] = rows[0]["byte_breakdown"]


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: python3 scripts/summarize.py RESULT_DIRECTORY RUNS")
    root, count = Path(sys.argv[1]), int(sys.argv[2])
    if count < 1:
        raise ValueError("RUNS must be positive")
    profiles = {}
    for name, info in PROFILES.items():
        if name.startswith("semi_"):
            rows = semi_runs(root / f"{name}.jsonl", info)
        else:
            rows = [malicious_run(root / f"{name}_run{i}.json") for i in range(1, count + 1)]
        if len(rows) != count:
            raise ValueError(f"wrong run count: {name}")
        profiles[name] = {"setting": info[0], "profile": info[1], "logq": info[2],
                          "m": info[3], "ell": info[4], "correctness": "PASS",
                          "mean": aggregate(rows), "runs": rows, "id": name}
        if not name.startswith("semi_"):
            add_fixture(profiles[name], json.loads((root / f"{name}_one_round.json").read_text()), count)

    scope = ("BlindEval: per-server mean, independent request verification, sequential server proofs. "
             "Finalize: client-side parallel verification within the total thread budget. "
             "Request: one message including the fresh commitment. PreProc: one server, one epoch, unamortized. "
             "With NIZKs, the headline Request row is a synthetic one-round client CLI replay at offset 9, "
             "excluding seed search and the deployment driver; fresh-coin end-to-end Request is also reported. "
             "Other algorithms use the epoch benchmark. KB: 1024 bytes; no transport framing.")
    summary = {"runs": count, "warmups": 1, "scope": scope, "profiles": profiles}
    (root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    with (root / "summary.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("setting", "profile", "logq", "algorithm", "mean_seconds", "mean_bytes", "mean_KB", "correctness", "runs", "scope"))
        for value in profiles.values():
            for algorithm in ALGORITHMS:
                writer.writerow((value["setting"], value["profile"], value["logq"], algorithm,
                                 value["mean"]["seconds"][algorithm], value["mean"]["bytes"][algorithm],
                                 value["mean"]["bytes"][algorithm] / 1024, "PASS", count,
                                 "fixed-one-round-client" if algorithm == "Request" and "fresh_request" in value else "epoch"))
            if "fresh_request" in value:
                r = value["fresh_request"]
                writer.writerow((value["setting"], value["profile"], value["logq"], "Request",
                                 r["seconds"], r["bytes"], r["bytes"]/1024, "PASS", count, "fresh-coins-end-to-end"))

    lines = ["# T-LeOPaRd Results", "", f"One warmup and {count} measured runs per configuration. Correctness: PASS.",
             "", scope, "", "Hardware and compiler: `system.txt`. Raw logs are retained.", ""]
    for value in profiles.values():
        mean = value["mean"]
        lines += [f"## {value['setting']} / {value['profile']} (q{value['logq']})", "",
                  "| Algorithm | Time (ms) | Communication (KB) |", "|---|---:|---:|"]
        for algorithm in ALGORITHMS:
            lines.append(f"| F.{algorithm} | {1000 * mean['seconds'][algorithm]:.3f} | {mean['bytes'][algorithm]/1024:.2f} |")
        if "fresh_request" in value:
            r = value["fresh_request"]
            lines += ["", f"Request above: fixed synthetic request, offset 9, one accepting round. "
                      f"Fresh-coin end-to-end comparison: {1000*r['seconds']:.3f} ms, "
                      f"{r['bytes']/1024:.2f} KB; rejection rounds: {r['rounds']}."]
        for key, title, factor, unit in (("request_seconds", "Request timing", 1000, "ms"),
                                         ("blind_seconds", "BlindEval timing per server", 1000, "ms"),
                                         ("request_bytes", "Request communication", 1/1024, "KB"),
                                         ("reply_bytes", "Reply communication per server", 1/1024, "KB")):
            if key not in mean:
                continue
            lines += ["", f"### {title}", "", f"| Component | {unit} |", "|---|---:|"]
            for component, number in mean[key].items():
                lines.append(f"| {component} | {factor * number:.3f} |")
        lines.append("")
    lines += ["Nested preparation/proving timers are included in their parent timers, not additional costs.", ""]
    (root / "summary.md").write_text("\n".join(lines))
    print("\n".join(lines[:10]))
    for name, value in profiles.items():
        mean = value["mean"]
        print(f"{name}: Request {mean['seconds']['Request']:.6f} s / {mean['bytes']['Request']/1024:.2f} KB; "
              f"BlindEval/server {mean['seconds']['BlindEval']:.6f} s / {mean['bytes']['BlindEval']/1024:.2f} KB; PASS")


if __name__ == "__main__":
    main()
