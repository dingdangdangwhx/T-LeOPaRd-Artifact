"""The parameter search of the performance section (``sec:performance``).

The widths depend on N = m + ell, the modulus depends on m and ell, and the
ranks m and ell depend on the modulus. Starting from m = ell = 1, we repeat

1. the widths s0, s, s1 of Theorem ``thm:MLWEP-to-HMLWEP`` with N = m + ell,
2. the smallest modulus q with p_u <= 2^-kappa,
3. the smallest m such that MLWE with Gaussian errors of width s0 is hard
   (server security), and the smallest ell such that MLWE with errors uniform
   in {-beta_r, ..., beta_r} is hard (client security),

until m and ell no longer change. The binding ranks n_s and n_c of the
commitments then follow from MSIS.
"""

from dataclasses import dataclass

from sage.all import RR, Integer, ceil, log, sqrt

from .hardness import smallest_secure_msis_rank
from .modulus import log_modulus
from .widths import gaussian_widths, max_queries_per_input, width_to_stddev

#: Bound beta_r on the entries of the client's randomness R.
BETA_R = Integer(1)

#: The PRF output has d h log2(p) bits, which must be at least this many.
MIN_OUTPUT_BITS = 128

MAX_ITERATIONS = 100


@dataclass(frozen=True)
class Parameters:
    prf: str
    mode: str
    kappa: int
    p: int
    d: int
    m: int
    ell: int
    h: int
    n_s: int
    n_c: int
    #: Base of the gadget matrix of the BP14 PRF, None for BLMR13.
    gadget_base: object
    s0: object
    s: object
    s1: object
    log_q: int


def raise_rounding_modulus(p, d, h):
    """Smallest admissible p, so that 1/p^(dh) <= 2^-128, keeping p if it is large enough."""
    if RR(d * h * log(p, 2)) < MIN_OUTPUT_BITS:
        p = ceil(Integer(2) ** (Integer(MIN_OUTPUT_BITS) / (d * h)))
    return p


def find_parameters(prf, mode, kappa, d, target_rhf, dimensions, beta_r=BETA_R, rule="both"):
    """Parameters for one setting, as in the tables of the performance section.

    ``prf`` is the PRF behind the message mapping H, "BP14" or "BLMR13".
    ``rule`` selects the modulus rule of ``t_leopard.modulus.log_modulus``.
    For BP14 the output has h = 1 ring elements, for BLMR13 h = m.
    """
    if prf not in ("BP14", "BLMR13"):
        raise ValueError(f"unknown PRF {prf!r}, expected 'BP14' or 'BLMR13'")
    q_max = max_queries_per_input(mode, kappa)
    p = Integer(4) if d < 128 else Integer(2)
    m_guess = ell_guess = Integer(1)

    for _ in range(MAX_ITERATIONS):
        h = Integer(1) if prf == "BP14" else m_guess
        p = raise_rounding_modulus(p, d, h)
        widths = gaussian_widths(q_max, d, m_guess + ell_guess, h)
        log_q = log_modulus(kappa, p, d, m_guess, ell_guess, h, widths.s, widths.s1, beta_r, rule=rule)
        m = dimensions.rank(target_rhf, d, log_q, width_to_stddev(widths.s0), uniform=False)
        ell = dimensions.rank(target_rhf, d, log_q, beta_r, uniform=True)
        if (m, ell) == (m_guess, ell_guess):
            break
        m_guess, ell_guess = m, ell
    else:
        raise RuntimeError(f"the search for m and ell did not converge for {prf}, {mode}, kappa = {kappa}")

    gadget_base = ceil(Integer(2) ** (log_q / m)) if prf == "BP14" else None

    # The scripts bound the norm of an MSIS solution that breaks binding by 2 sqrt(k),
    # where k counts the binding rows, the message rows and the randomness of the commitment.
    n_c = smallest_secure_msis_rank(target_rhf, d, log_q,
                                    lambda n: 2 * RR(sqrt(n + h * (ell + m) + 1 + ell)))
    key_rows = m if prf == "BP14" else 2 * m
    n_s = smallest_secure_msis_rank(target_rhf, d, log_q,
                                    lambda n: 2 * RR(sqrt(n + key_rows + m)))

    return Parameters(prf=prf, mode=mode, kappa=kappa, p=p, d=d, m=m, ell=ell, h=h,
                      n_s=max(m, n_s), n_c=max(ell, n_c), gadget_base=gadget_base,
                      s0=widths.s0, s=widths.s, s1=widths.s1, log_q=log_q)


def check_parameters(params, batch):
    """Sanity checks of the original script. Raises ValueError if one fails."""
    pr = params
    # Extraction from c_r and c_k (sec:com_instantiation) needs the messages to be small enough for q.
    if pr.prf == "BP14":
        limit = RR(Integer(2) ** pr.log_q / 4)
    else:
        limit = RR(Integer(2) ** (pr.log_q / 2) / 4)
    if (pr.n_c + (pr.h * (pr.ell + pr.m) + 1) * batch + pr.ell) * pr.d >= limit:
        raise ValueError(f"{pr}: q is too small to extract from c_r")
    if (pr.n_s + 2 * pr.m + pr.m) * pr.d >= limit:
        raise ValueError(f"{pr}: q is too small to extract from c_k")
    if RR(pr.d * pr.h * log(pr.p, 2)) < MIN_OUTPUT_BITS:
        raise ValueError(f"{pr}: 1/p^(dh) is not negligible")
