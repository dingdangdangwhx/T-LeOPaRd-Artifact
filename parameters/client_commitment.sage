"""Binary client UMC screen. This does not certify the composed proof."""
import contextlib
import io
import json
import os
import sys
sys.path.insert(0,os.path.dirname(os.path.realpath(__file__)))
from client_screen import calculate
load(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                  "../third_party/estimator.py"))
q = Integer(65437)
k = 14
d = 64
message_rank = 87
dim = (2*k+message_rank)*d
opening_sq = Integer(dim)
mu = q//2
error_sq = (2*k*d+1)*opening_sq
assert is_prime(q)
assert 4*error_sq < mu**2
main_q = Integer(2)**72-107
assert is_prime(main_q)
P = (q-1)//2
residual = k*d*P+1+mu+P
# A nonzero challenge difference is a unit modulo q_c: all irreducible
# factors have degree 32, while its integer resultant is at most 254^64.
factor_degree = Mod(q,128).multiplicative_order()
ring = PolynomialRing(GF(q), 'X')
factor_degrees = [factor.degree() for factor, _ in (ring.gen()**d+1).factor()]
assert factor_degrees == [32,32]
challenge_eta = Integer(127)
assert factor_degree == 32
assert (2*challenge_eta)**d < q**factor_degree
screen = calculate()
assert screen['commitment_q']==q and screen['commitment_rank']==k
cap = Integer(screen['conditional_iterations'][-1]['cap_squared'])
response_bound = ceil(sqrt(9*cap*(1+Integer(screen['outer_response_base'])**2)))
weak_quotient_bound = 2*response_bound
# If cbar*v = W (mod main_q), LNP22 extraction bounds W, not v.
# Integer lifting of cbar*(A*r+e+mu*R-C)-q_c*W now has no wrap.
lifted_bound = 2*challenge_eta*ceil(sqrt(d))*residual+q*weak_quotient_bound
assert 2*lifted_bound < main_q
fake_quotient = inverse_mod(q,main_q)
for coeff in range(1,17):
    value = (coeff*fake_quotient)%main_q
    value = min(value,main_q-value)
    assert value > weak_quotient_bound
result = {"q":int(q), "rank":k, "message_rank":message_rank,
          "opening_domain":"coefficientwise {0,1}; L2^2 <= dimension follows exactly",
          "opening_l2_squared":int(opening_sq),
          "opening_l2":float(sqrt(opening_sq)),
          "extraction_error_linf":float(sqrt(error_sq)), "scale":int(mu),
          "extraction_margin":float(mu/2-sqrt(error_sq)),
          "quotient_bits":0, "unmultiplied_residual_bound":int(residual),
          "challenge_eta":int(challenge_eta), "factor_degree":int(factor_degree),
          "factor_degrees":list(map(int,factor_degrees)),
          "conditional_response_bound":int(response_bound),
          "weak_quotient_bound":int(weak_quotient_bound),
          "lift_residual_bound":int(lifted_bound),
          "lifting_margin_ratio":float(main_q/(2*lifted_bound)),
          "security_scope":"algebra checks and legacy lattice estimator only"}
result["lwe"] = {}
for label, model in [("enum", BKZ.enum), ("sieve", BKZ.sieve)]:
    with contextlib.redirect_stdout(io.StringIO()):
        est = estimate_lwe(k*d, alphaf(sigmaf(RR(0.5)),q),q,
                           secret_distribution=(0,1),reduction_cost_model=model)
    result["lwe"][label] = {
        name:{"log2_rop":float(log(cost["rop"],2)),"delta":float(cost["delta_0"])}
        for name,cost in est.items() if "delta_0" in cost}
    assert result["lwe"][label]
    assert min(c["log2_rop"] for c in result["lwe"][label].values()) >= 128
beta = 2*sqrt(opening_sq)
delta = 2**(log(beta,2)**2/(4*k*d*log(q,2)))
result["sis_collision_l2"] = float(beta)
result["sis_legacy_delta"] = float(delta)
assert beta < q and delta <= 1.0044
print(json.dumps(result,indent=2,default=int))
