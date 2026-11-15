"""Gaussian widths of the server's errors, from Theorem ``thm:MLWEP-to-HMLWEP``.

The theorem reduces MLWE-PRF to Hint-MLWE-PRF. We take the smallest widths it
allows, with every condition met with equality:

* ``s0``, the width of the MLWE error in the reduction, from the first condition
  of ``eq:sig0-lobnd``,
* ``s``, the width of the pre-processing error e_s, from the ratio condition
  ``eq:sig0-to-s-upbnd``, unless the first condition of ``eq:slobnd`` asks for more,
* ``s1_bar`` and ``s1_B``, the widths of the two hint masks, from ``eq:sig1-to-s-bnd``.

The online error e'_s then has width ``s1 = sqrt(s1_bar^2 + s1_B^2)``.

All values are Gaussian widths. A width s has standard deviation s / sqrt(2 pi).
As in the performance section (``sec:performance``), the theorem's two error
widths are both set to s, and the reduction errors are spherical.
"""

from dataclasses import dataclass

from sage.all import RR, Integer, RealField, log, pi, sqrt

R256 = RealField(256)

#: Statistical distance per simulated sample, epsilon in the theorem.
EPSILON = Integer(2) ** -128

#: Slack epsilon_p of the theorem used for the paper's tables.
#: The value 0.02 gave the smallest widths. It is a 53-bit real, as in the
#: original script, which matters for the last printed digit.
PAPER_EPS_P = RR("0.02")


@dataclass(frozen=True)
class Widths:
    s0: object
    s: object
    s1_bar: object
    s1_B: object
    s1: object
    #: True if the eps_p condition of eq:slobnd, not the ratio, determines s.
    s_from_eps_p_bound: bool


def max_queries_per_input(mode, kappa):
    """Q_max of the theorem, the number of evaluations of one input.

    In the fully oblivious setting the total number of queries is limited to
    2^kappa, matching the correctness error. In the partially oblivious setting
    the server can limit the queries under one tag to 2^16 (``sec:performance``).
    """
    if mode == "OPRF":
        return Integer(2) ** kappa
    if mode == "POPRF":
        return Integer(2) ** 16
    raise ValueError(f"unknown mode {mode!r}, expected 'OPRF' or 'POPRF'")


def gaussian_widths(q_max, d, N, h, *, delta=Integer(1), eta=Integer(1),
                    eps_p=PAPER_EPS_P, beta_r=Integer(1)):
    """Smallest widths allowed by the theorem.

    ``N`` is the number of rows of the MLWE matrix, i.e. m + ell in the protocol,
    ``h`` the number of rows of the client's matrix R, and ``beta_r`` the bound on
    its entries. ``delta``, ``eta`` and ``eps_p`` are the free parameters of the
    theorem. Pass exact Sage numbers, since the arithmetic follows their type.
    """
    if eps_p > RR(1 - 1 / (1 + eta**2)):
        raise ValueError("eps_p is too large: the ratio bound eq:sig0-to-s-upbnd is not positive")

    s0 = RR(sqrt(RR(1 / pi) * (1 + 1 / delta) * R256(log(2 * N * d * (1 + 1 / EPSILON)))))

    ratio = RR(sqrt((1 - 1 / (1 + eta**2) - eps_p) / (1 + delta)))
    s = s0 / ratio
    s_eps_p_bound = RR(sqrt(log(2 * N * d + 4) / (pi * eps_p)))
    s_from_eps_p_bound = s_eps_p_bound > s
    if s_from_eps_p_bound:
        s = s_eps_p_bound

    s1_bar = eta * beta_r * d * RR(sqrt(h * N)) * s
    s1_B = eta * RR(sqrt(q_max)) * s
    s1 = RR(sqrt(s1_bar**2 + s1_B**2))
    return Widths(s0, s, s1_bar, s1_B, s1, s_from_eps_p_bound)


def width_to_stddev(width):
    """Standard deviation of a discrete Gaussian with the given width."""
    return RR(width / sqrt(2 * pi))
