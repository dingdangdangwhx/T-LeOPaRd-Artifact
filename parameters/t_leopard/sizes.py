"""Communication and public-key sizes in KB (1 KB = 8192 bits).

Per query, the client sends the commitment c_r to its randomness R in the
pre-processing phase, and C_x, the input commitment d_x and the proof pi_c
online. The server answers with v_k in the pre-processing phase, and with u_x
and the proof pi_s online. Its public key is the commitment c_k to the PRF key.

The commitments are the Regev-style commitments of ``sec:com_instantiation``
(Appendix F). Their binding part has n rows, from which the D low-order bits can
be dropped, and every message ring element takes one row, or two rows if its
coefficients are not small.
"""

from dataclasses import dataclass

from sage.all import Integer, round

BITS_PER_KB = 8192

#: D, the low-order bits dropped from the binding rows of c_k and c_r (``sec:performance``).
DROPPED_BITS = 12

#: Size assumed for the client's input commitment d_x, per query.
INPUT_COMMITMENT_KB = 1

#: Estimated LaBRADOR proof sizes per batch size, for both parties (``sec:performance``).
#: The original scripts also had an LNP22 estimate, but it was disabled, so
#: LaBRADOR is always used. The LNP22 estimate is in extras/lnp22_proof_size.py.
LABRADOR_PROOF_KB = {1: 45, 64: 79}
PROOF_SYSTEM = "labrador"


@dataclass(frozen=True)
class Costs:
    """Communication of one party for one batch of queries, in KB."""

    offline_per_query: object
    proof: object
    online_without_proof: object
    online_total: object


def kilobytes(bits):
    return round(Integer(bits) / BITS_PER_KB, 2)


def proof_bits(batch):
    if batch not in LABRADOR_PROOF_KB:
        raise ValueError(f"no LaBRADOR proof size for batch size {batch}, "
                         f"known sizes: {sorted(LABRADOR_PROOF_KB)}")
    return LABRADOR_PROOF_KB[batch] * BITS_PER_KB


def public_key_kb(prf, log_q, d, m, n_s):
    """Size of c_k. A BP14 key is small, a BLMR13 key is uniform and takes two rows per entry."""
    key_rows = {"BP14": m, "BLMR13": 2 * m}[prf]
    return kilobytes(n_s * d * (log_q - DROPPED_BITS) + key_rows * d * log_q)


def client_costs(log_q, d, m, ell, h, n_c, batch):
    """Client: c_r per query offline, then C_x, d_x and pi_c online."""
    c_r = n_c * d * (log_q - DROPPED_BITS) + (h * (ell + m) + 1) * d * log_q
    C_x = h * m * d * log_q * batch
    proof = proof_bits(batch)
    d_x = INPUT_COMMITMENT_KB * batch
    return Costs(
        offline_per_query=kilobytes(c_r),
        proof=kilobytes(proof),
        online_without_proof=d_x + kilobytes(C_x),
        online_total=d_x + kilobytes(C_x + proof),
    )


def server_costs(log_q, d, m, ell, h, batch):
    """Server: v_k per query offline, then u_x and pi_s online."""
    v_k = (ell + m) * d * log_q
    u_x = h * d * log_q * batch
    proof = proof_bits(batch)
    return Costs(
        offline_per_query=kilobytes(v_k),
        proof=kilobytes(proof),
        online_without_proof=kilobytes(u_x),
        online_total=kilobytes(proof + u_x),
    )
