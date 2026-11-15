"""Accepted coefficient-l2 budgets, separate from honest Gaussian widths.

For X sampled from D_{Z^n,s}, Poisson summation bounds the mgf of
2*pi*||X||^2/s^2 by that of chi-square(n). The usual Chernoff bound gives
Pr[||X||^2 > s^2/(2*pi)*(n+2*sqrt(n*t)+2*t)] <= exp(-t).
No Gaussian approximation and no infinity-to-l2-to-infinity round trip
is used here. This module does not certify the composed client proof.
"""
from dataclasses import dataclass
from sage.all import Integer, RealField, ceil, log, pi, sqrt

RF = RealField(256)


@dataclass(frozen=True)
class NoiseBudgets:
    preprocessing: object
    online: object
    digit_base: int
    digit_low: object
    digit_high: object
    tail_bits: int


def gaussian_l2(width, dimension, tail_bits=192):
    n, t = RF(dimension), RF(tail_bits)*RF(log(2))
    return RF(width)*sqrt((n+2*sqrt(n*t)+2*t)/(2*RF(pi)))


def noise_budgets(s, s1, rows, h, d, digit_base=None, tail_bits=192):
    """Two balanced digits; bounds include reconstruction relaxation.

    B1 >= ||d0|| + base*||d1|| covers every accepted digit witness.
    It is not the smaller norm of the honestly sampled online error.
    """
    bs = ceil(gaussian_l2(s, rows*d, tail_bits))
    honest = gaussian_l2(s1, h*d, tail_bits)
    base = Integer(digit_base) if digit_base is not None else max(
        Integer(2), ceil(sqrt(2*honest/sqrt(RF(h*d)))))
    if base < 2:
        raise ValueError("digit base must be at least two")
    low = ceil((base//2)*sqrt(RF(h*d)))
    high = ceil((honest+low)/base)
    return NoiseBudgets(bs, low+base*high, int(base), low, high, tail_bits)


def accepted_final_error(bs, b1, rows, d, beta_r, width, eta):
    return RF(width)*eta*(RF(beta_r)*sqrt(RF(rows*d))*bs+b1)
