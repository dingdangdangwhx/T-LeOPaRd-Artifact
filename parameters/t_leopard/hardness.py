"""Hardness of MSIS and MLWE, measured by the root Hermite factor (RHF).

Both problems are considered hard enough when the best attack needs an RHF of at
most the target, 1.0045 for about 128 bits and 1.0060 for about 90 bits
(``sec:performance``).

* MSIS uses the closed-form estimate of the LNP22 scripts.
* MLWE uses the 2022 version of the lattice estimator of Albrecht, Player and
  Scott. It is slow, so ranks are looked up in a table of precomputed
  dimensions first (``data/mlwe_dimensions.csv``), and every estimator call is
  memoized for the rest of the run.
"""

import csv
import hashlib
import importlib.util
import logging
import os
from functools import lru_cache
from pathlib import Path

from sage.all import RR, Integer, ceil, log, round, sqrt

LOG = logging.getLogger("t_leopard")

SCRIPTS_DIR = Path(__file__).resolve().parent.parent
DIMENSION_TABLE = SCRIPTS_DIR / "data" / "mlwe_dimensions.csv"
DEFAULT_ESTIMATOR_DIR = SCRIPTS_DIR.parent / "third_party"

#: SHA-256 of estimator.py at commit c50ab18d53ae of bitbucket.org/malb/lwe-estimator,
#: the version bundled with the artifact.
ESTIMATOR_SHA256 = "fcb4846093f1a2117ac8703b757b47c1121227cfa223199da2b10045044d6568"

#: Root Hermite factor targets per security level.
TARGET_RHF = {128: RR("1.0045"), 90: RR("1.0060")}


# --------------------------------------------------------------------------- MSIS

def msis_root_hermite_factor(bound, rank, d, log_q):
    """RHF needed to find an MSIS solution of Euclidean norm at most ``bound``.

    For a rank-``rank`` module over a degree-``d`` ring and modulus 2^log_q,
    log2(RHF) = log2(bound)^2 / (4 rank d log_q), following [GamNgu08, Mic08] as in
    the LNP22 scripts. For rank 0 the result is infinite, i.e. not secure.
    """
    log_bound = log(bound, 2)
    return RR(Integer(2) ** (log_bound**2 / (4 * rank * d * log_q)))


def smallest_secure_msis_rank(target_rhf, d, log_q, solution_norm):
    """Smallest rank n >= 0 whose MSIS instance needs an RHF of at most ``target_rhf``.

    ``solution_norm(n)`` gives the norm bound of the relevant solution, which may
    depend on n.
    """
    n = Integer(0)
    while msis_root_hermite_factor(solution_norm(n), n, d, log_q) > target_rhf:
        n += 1
    return n


# --------------------------------------------------------------------------- MLWE

@lru_cache(maxsize=1)
def lattice_estimator():
    """The pinned lattice estimator, imported from ``deps/`` on first use.

    Set the environment variable LWE_ESTIMATOR_DIR to use another copy.
    """
    directory = Path(os.environ.get("LWE_ESTIMATOR_DIR", DEFAULT_ESTIMATOR_DIR))
    path = directory / "estimator.py"
    if not path.is_file():
        raise RuntimeError(f"the bundled lattice estimator is missing ({path})")
    if hashlib.sha256(path.read_bytes()).hexdigest() != ESTIMATOR_SHA256:
        LOG.warning("%s is not the version used for the paper, results may differ", path)

    spec = importlib.util.spec_from_file_location("lwe_estimator_2022", path)
    estimator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(estimator)
    logging.getLogger("estimator").setLevel(logging.WARNING)  # it reports every estimate at INFO
    return estimator


def uniform_stddev(bound):
    """Standard deviation of the uniform distribution on {-bound, ..., bound}."""
    return RR(sqrt(((2 * bound + 1) ** 2 - 1) / 12))


@lru_cache(maxsize=None)
def mlwe_root_hermite_factor(noise, rank, d, log_q, uniform):
    """RHF of the best attack on MLWE, maximized over attacks and BKZ cost models.

    The instance has secret rank ``rank`` over a degree-``d`` ring, modulus
    ceil(2^log_q), and errors that are uniform on {-noise, ..., noise} if
    ``uniform`` is set, and Gaussian with standard deviation ``noise`` otherwise.
    The attacks are uSVP, decoding and dual, each with the enumeration and the
    sieving cost model of the estimator.
    """
    est = lattice_estimator()
    q = ceil(Integer(2) ** log_q)
    stddev = uniform_stddev(noise) if uniform else noise
    alpha = est.alphaf(est.sigmaf(stddev), q)

    rhfs = []
    for cost_model in (est.BKZ.enum, est.BKZ.sieve):
        costs = est.estimate_lwe(rank * d, alpha, q, reduction_cost_model=cost_model)
        rhfs += [costs[attack]["delta_0"] for attack in ("usvp", "dec", "dual")]
    return max(rhfs)


class MLWEDimensions:
    """Smallest secure MLWE ranks, from the precomputed table or the estimator.

    Table entries map (RHF, log q, noise) to the LWE dimension n d. The rules
    below reproduce the original script exactly:

    * a lookup uses round(noise), but a computed rank uses ceil(noise), so
      Gaussian lookups with a fractional part below 1/2 never hit the table,
    * a computed rank is added to the table for the rest of the run, keyed by
      ceil(log q) and ceil(noise), and the table on disk is never changed.
    """

    def __init__(self, path=DIMENSION_TABLE):
        self.table = {"uniform": {}, "gaussian": {}}
        with open(path, newline="") as f:
            for row in csv.DictReader(f):
                key = (RR(row["root_hermite_factor"]), Integer(row["log_q"]), Integer(row["noise"]))
                self.table[row["distribution"]][key] = Integer(row["dimension"])

    def rank(self, target_rhf, d, log_q, noise, uniform):
        """Smallest rank m such that MLWE with these parameters needs an RHF of at most target_rhf."""
        table = self.table["uniform" if uniform else "gaussian"]
        key = (round(target_rhf, 4), log_q, round(noise))
        if key in table:
            return int(ceil(table[key] / d))

        LOG.debug("  not in the table, running the lattice estimator: RHF %s, d = %s, log q = %s, noise %s (%s)",
                  target_rhf, d, log_q, noise, "uniform" if uniform else "Gaussian")
        log_q, noise = ceil(log_q), ceil(noise)

        def too_weak(rank):
            return mlwe_root_hermite_factor(noise, rank, d, log_q, uniform) > target_rhf

        # Search from total dimension 512: double, then go up in steps of 5, then of 1.
        rank = int(Integer(512) / d)
        while too_weak(rank):
            rank *= 2
        rank //= 2
        while too_weak(rank):
            rank += 5
        rank -= 5
        while too_weak(rank):
            rank += 1

        table[(round(target_rhf, 4), log_q, round(noise))] = rank * d
        return rank
