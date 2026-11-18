"""Section 4 protocol driver for the 2-of-3 artifact. C kernels are separate processes."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import statistics
import subprocess
import tempfile
import time

BIN = Path(__file__).resolve().parents[1] / "src/protocol"
CLIENT, SERVER = BIN / "threshold-client", BIN / "threshold-server"
LIMIT = {"partial72": 1024, "plain72": 4096}
TAG = "504f5052462d7461672d303030303100"


def call(tool, *args):
    p = subprocess.run([tool, *map(str, args)], capture_output=True, text=True)
    fields = dict(line.split("=", 1) for line in p.stdout.splitlines()
                  if "=" in line and line.partition("=")[0].isidentifier())
    if p.returncode < 0:
        raise RuntimeError(f"{args[0]} terminated by signal {-p.returncode}")
    if p.returncode or fields.get("ok") != "1":
        raise ValueError(f"{args[0]} rejected: {p.stderr[-500:]}")
    return fields


def shake(domain, data, n):
    return hashlib.shake_128(domain + data).digest(n)


class Protocol:
    def __init__(self, root, profile):
        self.root = root
        self.profile = profile
        self.pp, self.pk, self.sk = [root / name for name in ("pp", "pk", "sk")]
        self.master = root / "master"
        self.secret = root / "eval-secret"
        if int(call(SERVER, "epoch-limit")["limit"]) != LIMIT[profile]:
            raise ValueError("compiled server does not match requested profile")
        self.setup = call(CLIENT, "epoch-setup", self.pp)
        start = time.perf_counter()
        self.keys = call(SERVER, "keygen", self.pp, 3, 2, self.pk, self.sk, self.master)
        self.secret.write_bytes(os.urandom(16))
        seed_hash = shake(b"T-LeOPaRd-ROs", self.secret.read_bytes(), 32)
        (self.pk / "seed-hash").write_bytes(seed_hash)
        parts = {}
        for path in self.pk.glob("pk_*.bin"):
            wire = path.read_bytes()
            mask = int.from_bytes(wire[4:12], "little")
            if wire != parts.setdefault(mask, wire):
                raise ValueError("replicated component commitment mismatch")
        if set(parts) != {3, 5, 6}:
            raise ValueError("missing key components")
        self.pkhash = shake(b"T-LeOPaRd-PK-v1", b"".join(parts[i] for i in sorted(parts)) + seed_hash, 32)
        self.key_time = time.perf_counter() - start
        self.db = sqlite3.connect(root / "epochs.sqlite", isolation_level=None)
        self.db.execute("CREATE TABLE epochs(id BLOB, party INTEGER, uses INTEGER, PRIMARY KEY(id,party))")

    def seed(self, epoch):
        if len(epoch) != 8:
            raise ValueError("invalid epoch identifier")
        return shake(b"T-LeOPaRd-epoch-v1", self.pkhash + epoch, 16).hex()

    def directories(self, epoch):
        root = self.root / self.seed(epoch)
        return root / "published", root / "private"

    def preproc(self, epoch, party):
        pub, private = self.directories(epoch)
        if party not in range(3):
            raise ValueError("invalid party")
        # Reserve before sampling. A failed or interrupted publication must never
        # allow a second independent sample for the same matrix and share.
        self.db.execute("INSERT INTO epochs VALUES(?,?,-1)", (epoch, party))
        pub.parent.mkdir(exist_ok=True)
        result = call(SERVER, "epoch-preproc", self.pp, self.sk, 3, 2,
                      self.seed(epoch), pub, private, party)
        packet = epoch + b"".join((pub / f"prep_{party}_{i}.bin").read_bytes() for i in range(2))
        if len(packet) != int(result["prep_total_bytes"]) + len(epoch):
            raise ValueError("incomplete epoch publication")
        (pub / f"server-{party}.packet").write_bytes(packet)
        self.db.execute("UPDATE epochs SET uses=0 WHERE id=? AND party=?", (epoch, party))
        return len(packet)

    def request(self, epoch, quorum, tag, x, client):
        pub, _ = self.directories(epoch)
        if len(set(quorum)) != 2 or any(i not in range(3) for i in quorum):
            raise ValueError("invalid quorum")
        saved = client / "publication"
        saved.mkdir()
        for party in quorum:
            packet = (pub / f"server-{party}.packet").read_bytes()
            if packet[:8] != epoch or (len(packet) - 8) % 2:
                raise ValueError("inconsistent epoch publication")
            # Both component vectors use the same fixed-width wire format.
            size = (len(packet) - 8) // 2
            for i in range(2):
                (saved / f"prep_{party}_{i}.bin").write_bytes(packet[8+i*size:8+(i+1)*size])
        result = call(CLIENT, "epoch-request", self.pp, self.seed(epoch), 3, 2,
                      ",".join(map(str, quorum)), tag, x, client / "req", client / "st")
        (client / "context").write_text(json.dumps([epoch.hex(), tag, self.pkhash.hex()]))
        return epoch + (client / "req").read_bytes(), result

    def blind(self, party, message, tag):
        epoch = message[:8]
        _, private = self.directories(epoch)
        self.db.execute("BEGIN IMMEDIATE")
        try:
            row = self.db.execute("SELECT uses FROM epochs WHERE id=? AND party=?", (epoch, party)).fetchone()
            if row is None or not 0 <= row[0] < LIMIT[self.profile]:
                raise ValueError("unknown or exhausted epoch")
            with tempfile.TemporaryDirectory(dir=self.root) as tmp:
                tmp = Path(tmp)
                req = tmp / "incoming"
                req.write_bytes(message[8:])
                verified = call(CLIENT, "epoch-verify", self.pp, tag, req, self.seed(epoch))
                proved = call(SERVER, "epoch-prove-reply", self.pp, self.pk, self.sk,
                              private, tag, req, tmp / "out", party)
                reply = (tmp / "out" / f"repb_{party}_0.bin").read_bytes()
            self.db.execute("UPDATE epochs SET uses=uses+1 WHERE id=? AND party=?", (epoch, party))
            self.db.commit()
            return reply, verified, proved
        except BaseException:
            self.db.rollback()
            raise

    def finalize(self, epoch, tag, client):
        pub = client / "publication"
        if json.loads((client / "context").read_text()) != [epoch.hex(), tag, self.pkhash.hex()]:
            raise ValueError("request context mismatch")
        call(SERVER, "verify-replies", self.pp, self.pk, pub, tag, client / "req", client)
        return call(CLIENT, "epoch-finalize", client / "st", pub, client, self.pkhash.hex())["y_hex"]

    def evaluate(self, tag, x):
        return call(CLIENT, "epoch-eval", self.pp, tag, x, self.master,
                    self.pkhash.hex(), self.secret)["y_hex"]


def timed(fn, *args):
    start = time.perf_counter()
    value = fn(*args)
    return value, time.perf_counter() - start


def run(profile, regress=False):
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        p = Protocol(root, profile)
        epoch = (1).to_bytes(8, "big")
        offline = [timed(p.preproc, epoch, i) for i in range(3)]
        seconds = {"Setup": float(p.setup["time_sec"]), "KeyGen": p.key_time,
                   "PreProc": statistics.mean(t for _, t in offline)}
        sizes = {"Setup": int(p.setup["pp_bytes"]), "KeyGen": int(p.keys["pk_total_bytes"]) + 32,
                 "PreProc": statistics.mean(b for b, _ in offline)}
        previous = None
        def rejects(fn, *args):
            try:
                fn(*args)
            except (ValueError, sqlite3.IntegrityError):
                return
            raise AssertionError("invalid operation accepted")
        for index, quorum in enumerate(((0, 2), (0, 1), (1, 2)) if regress else ((0, 2),)):
            client = root / f"client-{index}"
            client.mkdir()
            x = "5ac31980277104bed261"
            (msg, rq), rt = timed(p.request, epoch, quorum, TAG, x, client)
            cr = msg[116:116 + int(rq["req_cr_wire_bytes"])]
            if not cr or cr == previous:
                raise AssertionError("client commitment missing or reused")
            previous = cr
            if regress and index == 0:
                rejects(p.preproc, epoch, 0)
                rejects(p.blind, 0, bytes(8) + msg[8:], TAG)
                rejects(p.blind, 0, msg, "00" * 16)
                rejects(p.blind, 0, msg[:-1], TAG)
                other = (2).to_bytes(8, "big")
                p.preproc(other, 0)
                rejects(p.blind, 0, other + msg[8:], TAG)
                assert p.db.execute("SELECT uses FROM epochs WHERE party=0").fetchone() == (0,)
            replies, bt, vs, ps, public, proofs = [], [], [], [], [], []
            for party in quorum:
                (reply, v, s), elapsed = timed(p.blind, party, msg, TAG)
                (client / f"repb_{party}_0.bin").write_bytes(reply)
                public.append(int(s["rep_public_total_bytes"]))
                proofs.append(int(s["rep_proof_total_bytes"]))
                if public[-1] + proofs[-1] != len(reply):
                    raise ValueError("reply byte accounting mismatch")
                replies.append(len(reply))
                bt.append(elapsed)
                vs.append(float(v["time_sec"]))
                ps.append(float(s["blindeval_prover_total_sec"]))
            final, ft = timed(p.finalize, epoch, TAG, client)
            if regress and index == 0:
                path = client / f"repb_{quorum[0]}_0.bin"
                intact = path.read_bytes()
                path.write_bytes(intact[:-1])
                rejects(p.finalize, epoch, TAG, client)
                path.write_bytes(intact)
            direct, et = timed(p.evaluate, TAG, x)
            if final != direct or regress and direct != p.evaluate(TAG, x):
                raise AssertionError("Finalize/Eval mismatch")
            if index == 0:
                seconds.update(Request=rt, BlindEval=statistics.mean(bt), Finalize=ft, Eval=et)
                sizes.update(Request=len(msg), BlindEval=statistics.mean(replies), Finalize=0, Eval=0)
                detail = dict(request=rq, verify_request=vs, nizk2_prove=ps,
                              reply_public=statistics.mean(public), reply_proofs=statistics.mean(proofs))
        if regress:
            # Move only the test ledger to its last admissible use.
            p.db.execute("UPDATE epochs SET uses=? WHERE party=1", (LIMIT[profile]-1,))
            p.blind(1, msg, TAG)
            rejects(p.blind, 1, msg, TAG)
        p.db.close()
        return dict(seconds=seconds, bytes=sizes, detail=detail, correctness="PASS")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("profile", choices=LIMIT)
    ap.add_argument("--test", action="store_true")
    args = ap.parse_args()
    print(json.dumps(run(args.profile, args.test)))
