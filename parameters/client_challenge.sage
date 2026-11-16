"""Interval Chernoff lower bound for the ideal eta-filtered challenge set.

The C implementation uses long-double cosine with a 1e-6 margin. This script
certifies the mathematical filter at 126.999, not a C libm error theorem.
"""
import json
F = RealIntervalField(256)
t = F(1)/10
threshold = F(126999)/1000
bad = F(0)
for root in range(32):
    theta = F.pi()*(2*root+1)/64
    mgf = sum((t*a).exp() for a in range(-8,9))/17
    for i in range(1,32):
        weight = 2*(theta*i).cos()
        mgf *= sum((t*weight*a).exp() for a in range(-8,9))/17
    bad += 2*(-t*threshold).exp()*mgf
accepted = 1-bad
assert accepted.lower() > 0
bits = (accepted*Integer(17)**32).log2()
assert bits.lower() > 129
print(json.dumps({
    "ideal_operator_norm_threshold":float(threshold.center()),
    "union_bad_probability_upper":float(bad.upper()),
    "accepted_fraction_lower":float(accepted.lower()),
    "log2_challenge_count_lower":float(bits.lower()),
    "scope":"Ideal real-arithmetic filter; not a simulation or FS theorem"
},indent=2))
