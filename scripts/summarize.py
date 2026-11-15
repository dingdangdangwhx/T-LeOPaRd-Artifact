"""Summarize the four artifact benchmarks using only the Python standard library."""
import csv
import json
import math
import re
import statistics
import sys
from pathlib import Path

ALGORITHMS = ("Setup", "KeyGen", "PreProcClient", "PreProcServer",
              "Request", "BlindEval", "Finalize", "Eval")
PROFILES = {
    "partial72": ("With NIZKs", "Partial", 72, 42, 45),
    "plain72": ("With NIZKs", "Plain", 72, 42, 45),
    "semi_partial56": ("Semi-honest", "Partial", 56, 32, 35),
    "semi_plain64": ("Semi-honest", "Plain", 64, 37, 40),
}


def field(text, name):
    match = re.search(r"^\s*" + re.escape(name) + r"\s*=\s*([0-9.eE+-]+)", text, re.M)
    if match is None:
        raise ValueError(f"missing measurement: {name}")
    value = float(match[1])
    if not math.isfinite(value) or value < 0:
        raise ValueError(f"invalid measurement: {name}")
    return value


def section(text, name, end):
    return text.split(name + ":\n", 1)[1].split(end + ":\n", 1)[0]


def malicious_run(path):
    text = path.read_text()
    if not re.search(r"correctness\s*= PASS\b", text):
        raise ValueError(f"correctness failed: {path.name}")
    if not re.search(r"threshold mode\s*= 2-of-3", text):
        raise ValueError("expected 2-of-3 sharing")
    if not re.search(r"reply bundles\s*= 2 bundles covering 3", text):
        raise ValueError("expected two aggregate proofs covering all components")
    timing = section(text, "Timing", "Finalize Breakdown")
    req = section(text, "Request Breakdown", "BlindEval Breakdown")
    blind = section(text, "BlindEval Breakdown", "Communication")
    comm = section(text, "Communication", "Outputs")
    result = {"correctness": "PASS", "seconds": {}, "bytes": {}}
    byte_names = dict(zip(ALGORITHMS, ("pp bytes", "pk bytes", "preq bytes",
                                      "prep bytes", "req bytes", "rep bytes", None, None)))
    for name in ALGORITHMS:
        result["seconds"][name] = field(timing, "F." + name) / (2 if name == "BlindEval" else 1)
        message_bytes = field(comm, byte_names[name]) if byte_names[name] else 0
        result["bytes"][name] = message_bytes / (2 if name == "BlindEval" else 1)
    result["request_seconds"] = {key: field(req, key) for key in (
        "protocol total", "NIZK total", "inner instance", "inner prepare",
        "inner prove", "outer witness", "outer statement", "outer prove", "serialize proof")}
    result["blind_seconds"] = {key: field(blind, key) for key in (
        "verify-request", "public ctx", "FS + ar-seed", "deserialize req", "outer statement", "outer verify")}
    for key in ("protocol total", "NIZK2 core total", "build instance", "NIZK2 prover", "export reply"):
        result["blind_seconds"][key] = field(blind, key) / 2
    result["request_bytes"] = {key: field(comm, key) for key in (
        "Cx bytes", "other public", "req proof bytes")}
    result["reply_bytes"] = {key: field(comm, key) / 2 for key in (
        "rep public bytes", "rep proof bytes")}
    if abs(sum(result["request_bytes"].values()) - result["bytes"]["Request"]) > 1e-6:
        raise ValueError("request communication components do not sum to the total")
    if abs(sum(result["reply_bytes"].values()) - result["bytes"]["BlindEval"]) > 1e-6:
        raise ValueError("reply communication components do not sum to the total")
    result["inner_rejection_rounds"] = field(req, "inner reject rnds")
    result["inner_batch_rounds"] = field(req, "inner batch rnds")
    return result


def semi_runs(path, expected):
    rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    for row in rows:
        if row["correctness"] != "PASS" or row["quorums"] != 3 or row["proofs_executed"] != 0:
            raise ValueError(f"proof-free correctness failed: {path.name}")
        if (row["logq"], row["m"], row["ell"]) != expected[2:]:
            raise ValueError(f"unexpected profile: {path.name}")
        if len(row["seconds"]) != 8 or len(row["bytes"]) != 8:
            raise ValueError("incomplete proof-free measurements")
        if any(not math.isfinite(x) or x < 0 for x in row["seconds"] + row["bytes"]):
            raise ValueError("invalid proof-free measurement")
    return [{"correctness": "PASS", "seconds": dict(zip(ALGORITHMS, row["seconds"])),
             "bytes": dict(zip(ALGORITHMS, row["bytes"])),
             "request_seconds": dict(zip(("fast core", "selector matrices", "load slot", "chain", "Cx"),
                                          row["request_breakdown"]))} for row in rows]


def aggregate(rows):
    means = {}
    for group in ("seconds", "bytes", "request_seconds", "blind_seconds", "request_bytes", "reply_bytes"):
        if group in rows[0]:
            means[group] = {key: statistics.mean(row[group][key] for row in rows) for key in rows[0][group]}
    return means


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
            rows = [malicious_run(root / f"{name}_run{i}.log") for i in range(1, count + 1)]
        if len(rows) != count:
            raise ValueError(f"wrong run count: {name}")
        profiles[name] = {"setting": info[0], "profile": info[1], "logq": info[2],
                          "m": info[3], "ell": info[4], "correctness": "PASS",
                          "mean": aggregate(rows), "runs": rows}

    scope = ("BlindEval: per-server mean, independent request verification, sequential server proofs. "
             "Finalize: client-side parallel verification. Request: one message. "
             "PreProcClient: 10 instances. PreProcServer: 3 servers x 10 instances. "
             "With-NIZK Eval: PRF arithmetic only. KB: 1024 bytes; communication excludes transport overhead.")
    summary = {"runs": count, "warmups": 1, "scope": scope, "profiles": profiles}
    (root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    with (root / "summary.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("setting", "profile", "logq", "algorithm", "mean_seconds", "mean_bytes", "mean_KB", "correctness", "runs"))
        for value in profiles.values():
            for algorithm in ALGORITHMS:
                writer.writerow((value["setting"], value["profile"], value["logq"], algorithm,
                                 value["mean"]["seconds"][algorithm], value["mean"]["bytes"][algorithm],
                                 value["mean"]["bytes"][algorithm] / 1024, "PASS", count))

    lines = ["# T-LeOPaRd Results", "", f"One warmup and {count} measured runs per configuration. Correctness: PASS.",
             "", scope, "", "Hardware and compiler: `system.txt`. Raw logs are retained.", ""]
    for value in profiles.values():
        mean = value["mean"]
        lines += [f"## {value['setting']} / {value['profile']} (q{value['logq']})", "",
                  "| Algorithm | Time (ms) | Communication (KB) |", "|---|---:|---:|"]
        for algorithm in ALGORITHMS:
            lines.append(f"| F.{algorithm} | {1000 * mean['seconds'][algorithm]:.3f} | {mean['bytes'][algorithm]/1024:.2f} |")
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
