#!/usr/bin/env python3
"""Recompute the q72 binary-opening geometry and conditional inner MSIS screen.

This is not a simulation proof, an end-to-end 128-bit certificate, or a modern
lattice attack survey. The commitment audit is a separate Sage calculation.
"""
import json
import math


def pair_budget(coeffs, base, logsigma):
    denominator = base * base
    return (coeffs * denominator // 4
            + (16 * coeffs * 2**(2*logsigma) + denominator - 1) // denominator
            + coeffs)


def calculate(base=40960):
    d, q, m, ell, depth, radix = 64, 2**72-107, 42, 45, 80, 319
    rank, qc = 14, 65437
    rhat = m + ell
    rho_dim = (2*rank+rhat)*d
    binary = rhat*d + rho_dim
    scalar = 0
    rbit = binary + scalar
    b_len = depth*m*9
    tail = rhat + rank + rhat + rank + rank + rhat
    m1 = b_len + depth + tail + rbit
    logsigma1, logsigma2 = 27, 17
    kmlwe, const_rows, eta, slack = 45, 4, 140, 3
    iterations = []
    k = 21
    for _ in range(100):
        m2 = kmlwe+k+const_rows+1
        budget = pair_budget(b_len*d, base, logsigma1)
        budget += pair_budget((depth+tail+rbit)*d, base, logsigma1)
        budget += pair_budget(m2*d, 2048, logsigma2)
        ptop = (q//2)//32768**4 + 2
        budget += depth*d*(4*32768**2//4+ptop*ptop)
        budget += 16*d
        responses = slack*math.sqrt(budget*(1+base*base))
        beta = 8*eta*responses
        next_k = math.ceil(math.log2(beta)**2/(4*d*math.log2(q)*math.log2(1.0044)))
        iterations.append(dict(kmsis=k, m2=m2, cap_squared=budget,
                               reconstructed_response_l2=responses,
                               msis_beta=beta, required_rank=next_k))
        if next_k == k:
            break
        k = next_k
    else:
        raise RuntimeError("parameter iteration did not converge")
    assert (m1, rbit) == (43551, 12928)
    assert beta < q
    result = dict(main_q=q, commitment_q=qc, commitment_rank=rank,
                  m=m, ell=ell, depth=depth, gadget_base=radix,
                  opening_l2_squared=rho_dim, opening_domain="{0,1}", coordinate_bits=1,
                  quotient_bits=0, slack_bits=0, m1=m1,
                  outer_response_base=base, outer_z2_base=2048,
                  outer_pbit_base=0, outer_pbit_digits=0,
                  outer_witness_rank=2*m1+400+2*m2,
                  outer_frontend_blocks=12,
                  outer_layout="linear-tail concatenation; norm-preserving permutation",
                  rbit_and_norm_scalar_rank=rbit, binary_rank=binary,
                  norm_scalar_rank=scalar, sigma1=1.55*2**logsigma1,
                  sigma2=1.55*2**logsigma2, eta=eta, actual_challenge_eta=127,
                  honest_s1_envelope_l2=1.55*2**logsigma1/(10*eta),
                  honest_s2_worst_l2=math.sqrt(m2*d),
                  sigma2_required=20*eta*math.sqrt(m2*d),
                  conditional_iterations=iterations,
                  status="CONDITIONAL_PARAMETER_SCREEN_NOT_SECURITY_CERTIFICATION",
                  outstanding=[
                      "Complete custom LNP/outer simulation and FS/query-loss analysis",
                      "Rigorous C arithmetic for the filter; ideal cardinality audited separately",
                      "JL and aggregation error accounting for the complete recursive proof",
                      "Justification of the honest-chain norm envelope/distribution",
                      "Independent review of the legacy lattice estimates",
                  ])
    arp_limit=2**25
    arp_bound=2*arp_limit
    assert rbit*arp_bound*(arp_bound+1) < q
    assert 1.55*2**20 >= 10*math.sqrt(128)*rbit
    result['arp']=dict(rows=128,sigma=1.55*2**20,limit=arp_limit,
        extracted_linf=arp_bound,nonnegative_sum_upper=rbit*arp_bound*(arp_bound+1),
        nowrap_margin=q/(rbit*arp_bound*(arp_bound+1)),
        message_commitment_rows=2,conditional_soundness_bits=128)
    assert result["sigma2"] >= result["sigma2_required"]
    return result


if __name__ == "__main__":
    print(json.dumps(calculate(), indent=2))
