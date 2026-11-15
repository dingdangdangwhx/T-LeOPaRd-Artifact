"""The modulus q.

The performance section (``sec:performance``) sets log q as small as possible
such that two conditions hold (``rule="both"``, the default).

1. The uniqueness bound of Theorem ``thm:uniqueness``,

       p_u = hd (2 B_f + 1) / (q / p),

   is at most 2^-kappa, where B_f is the worst-case bound ``Bf-worstcase`` on a
   coordinate of the final error e_f = e'_s - R e_s:

       B_f = (m + ell) d beta beta_r + beta_1.

2. The sufficient condition q/p >= 2^(kappa+2) hd (2 B_f + 1) of the correctness
   Lemma ``lemma:correctness`` holds with the Gaussian bound ``Bf-Gaussian``,

       B_f = sqrt(ln(2) (kappa + 3 + log2(hd)) / pi) (s1 + beta_r s sqrt((ell + m) d)),

   whose two-sided tail probability 2 exp(-pi B_f^2 / w^2) is at most the
   failure probability 2^-(kappa + 2 + log2(hd)) that the lemma requires.

Condition 1 alone (``rule="uniqueness"``) is the rule of the first version of
the paper. It is two bits below condition 2 with the worst-case bound, and the
Gaussian bound makes up for the two bits in all parameter sets of the paper
except the fully oblivious ones with kappa in {32, 64}, where condition 2 raises
log q by one and two bits.
"""

from sage.all import Integer, RealField, ceil, log, pi, round, sqrt

R256 = RealField(256)

MODULUS_RULES = ("both", "uniqueness")

#: The NIZKs prove ||e_s||_inf <= beta = 5 s and ||e'_s||_inf <= beta_1 = 5 s1.
#: A Gaussian sample of width w exceeds 5 w with probability about 2^-112.
TAIL_FACTOR = 5


def final_error_bound(m, ell, d, s, s1, beta_r):
    """Worst-case bound B_f on a coordinate of e_f = e'_s - R e_s."""
    beta = TAIL_FACTOR * s
    beta_1 = TAIL_FACTOR * s1
    return (m + ell) * d * beta * beta_r + beta_1


def gaussian_error_bound(kappa, d, m, ell, h, s, s1, beta_r):
    """Gaussian bound B_f of ``Bf-Gaussian`` on a coordinate of e_f, failing with probability 2^-(kappa+2+log2(hd))."""
    tail = R256(sqrt(log(2) * (kappa + 3 + log(h * d, 2)) / pi))
    return tail * (R256(s1) + beta_r * R256(s) * R256(sqrt((ell + m) * d)))


def uniqueness_modulus_bound(kappa, p, d, m, ell, h, s, s1, beta_r):
    """Smallest q with p_u <= 2^-kappa, as a real number."""
    bf = R256(final_error_bound(m, ell, d, s, s1, beta_r))
    return R256(p * R256(Integer(2) ** kappa) * h * d * (2 * bf + 1))


def correctness_modulus_bound(kappa, p, d, m, ell, h, s, s1, beta_r):
    """Smallest q allowed by ``eq:qpratio-stat`` with the Gaussian bound, as a real number."""
    bf = gaussian_error_bound(kappa, d, m, ell, h, s, s1, beta_r)
    return R256(p * R256(Integer(2) ** (kappa + 2)) * h * d * (2 * bf + 1))


def log_modulus(kappa, p, d, m, ell, h, s, s1, beta_r, rule="both"):
    """log2 of the smallest power-of-two modulus allowed by the rule.

    As in the original script, log2 of the bound on q is rounded to two decimals
    before it is rounded up to an integer.
    """
    bounds = [uniqueness_modulus_bound(kappa, p, d, m, ell, h, s, s1, beta_r)]
    if rule == "both":
        bounds.append(correctness_modulus_bound(kappa, p, d, m, ell, h, s, s1, beta_r))
    elif rule != "uniqueness":
        raise ValueError(f"unknown modulus rule {rule!r}, expected one of {MODULUS_RULES}")
    return ceil(round(log(max(bounds), 2), 2))
