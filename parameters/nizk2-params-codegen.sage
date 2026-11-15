import sys
import hashlib

assert len(sys.argv) == 2
params_file = sys.argv[1]

load(params_file)

assert vname != ""
assert deg == 64
assert len(dim) == 2
assert len(wpart) == len(wl2)
assert len(wpart) == len(wbin)

name = f"_{vname}"
d = 64
k = Integer(deg / d)
message_indices = list(globals().get("message_indices", []))
short_indices = [j for part in wpart for j in part]
wdim = dim[1]
assert sorted(short_indices + message_indices) == list(range(wdim))

fixed_q = mod
m1 = k * len(short_indices)
l = k * len(message_indices)
# Main-ring equations use the toolbox's exact algebraic interface. Only an
# explicit small-modulus quotient needs an approximate range proof.
nprime = k * globals().get("arp_rows", 0)
Bprime = globals().get("arp_bound", 0)
statistical_security = globals().get("statistical_security", 130)

wl2_ = wl2.copy()
for i in range(len(wpart)):
    if wl2_[i] == 0:
        wl2_[i] = sqrt(len(wpart[i]) * deg)
alpha = ceil(sqrt(sum(x ** 2 for x in wl2_)))

nbin = 0
for i in range(len(wpart)):
    if wbin[i] == 1:
        nbin += k * len(wpart[i])

n = []
B = []
for i in range(len(wpart)):
    if wl2[i] > 0:
        n += [k * len(wpart[i])]
        B += [wl2[i]]

verbose = 1
code = 1
loaded = 1
codegen_err = 0
load("lnp-tbox-codegen.sage")
assert codegen_err == 0
assert eknow < mpf(2)**(-128)
print(f"// exact main-ring equations: {dim[0] - nprime}; ARP quotient rows: {nprime}")
if globals().get("double_commitment", False):
    # Trapdoor S,E have coefficients at most one. Bound E*r-S*e_A+e_j
    # for every opening in the accepted Euclidean ball, not just sampled ones.
    scale = Integer(fixed_q).isqrt()
    extraction_error = sqrt(2 * M * deg + 1) * wl2[0]
    assert 2 * extraction_error < scale
    assert 2 * (scale + 1) * extraction_error < fixed_q
    print(f"// full-message double commitment: a={scale}, error<={extraction_error}")
    print(f"// extraction: 2E<a and 2(a+1)E<q; key block is unbounded")
if globals().get("small_commitment", False):
    # Include the ARP's extracted (not honest) quotient bound.
    assert small_residual_bound + COM_Q * psi * arp_bound < fixed_q / 2
    print(f"// small commitment: q={COM_Q}, rank={COM_RANK}, scale={COM_SCALE}")
    print(f"// quotient bound={arp_bound}, extracted quotient bound={psi * arp_bound}")
    print(f"// extraction error <= {extraction_error}; exact main-ring equations required")
    print(f"// lifted no-wrap lhs={small_residual_bound + COM_Q * psi * arp_bound}, rhs={fixed_q / 2}")

s1_indices = []
Ps_indices = []
Es_indices = []

for i in range(len(wpart)):
    if wbin[i] == 1:
        Ps_indices += list(range(len(s1_indices), len(s1_indices) + len(wpart[i])))
    if wl2[i] > 0:
        Es_indices += [list(range(len(s1_indices), len(s1_indices) + len(wpart[i])))]
    s1_indices += wpart[i]

tmp = []
for i in Ps_indices:
    tmp += list(range(i * k, i * k + k))
Ps_indices = tmp
Ps_indices.sort()

tmp2 = []
for j in Es_indices:
    tmp3 = []
    for i in j:
        tmp3 += list(range(i * k, i * k + k))
    tmp2 += [tmp3]
Es_indices = tmp2
Es_indices.sort()

out = ""

# Match LaZer's dense, row-major encoding, including coder_enc_end per group.
# These selectors depend only on the compiled parameters, never on a request.
qbits = Integer(mod).nbits()
if qbits % 8 == 0 and Es_indices:
    selector_hash = hashlib.shake_128()
    polybytes = int(d * (qbits // 8))
    for indices in Es_indices:
        encoded = bytearray(int(len(indices) * len(s1_indices) * polybytes) + 1)
        for row, col in enumerate(indices):
            encoded[int((row * len(s1_indices) + col) * polybytes)] = 1
        encoded[-1] = 1
        selector_hash.update(encoded)
    digest = selector_hash.digest(32)
    out += "#define NIZK2_FIXED_L2_HASH_AVAILABLE 1\n"
    out += "static const unsigned char nizk2_fixed_l2_hash[32] = { "
    out += ", ".join(str(x) for x in digest) + " };\n"

if Ps_indices == []:
    matPs_nrows = 0
    vname_Ps = "NULL"
else:
    out += f"static const unsigned int {vname}_Ps[{len(Ps_indices)}] = {intlist2intarray(Ps_indices)};\n"
    matPs_nrows = len(Ps_indices)
    vname_Ps = f"{vname}_Ps"

vname_Es = []
for i in range(len(Es_indices)):
    _Es_indices = Es_indices[i]
    if _Es_indices == []:
        vname_Es += ["NULL"]
    else:
        out += f"static const unsigned int {vname}_Es{i}[{len(_Es_indices)}] = {intlist2intarray(_Es_indices)};\n"
        vname_Es += [f"{vname}_Es{i}"]

if Es_indices == []:
    strEs = "NULL"
    strEs_nrows = "NULL"
else:
    strEs = f"{vname}_Es"
    strEs_nrows = f"{vname}_Es_nrows"
    out += f"static const unsigned int *{vname}_Es[{len(Es_indices)}] = {{ "
    for i in range(len(Es_indices)):
        out += f"{vname_Es[i]}, "
    out += "};\n"
    out += f"static const unsigned int {vname}_Es_nrows[{len(Es_indices)}] = {intlist2intarray([len(i) for i in Es_indices])};\n"

if message_indices:
    out += f"static const unsigned int {vname}_m_indices[{len(message_indices)}] = {intlist2intarray(message_indices)};\n"
    message_ptr = f"{vname}_m_indices"
else:
    message_ptr = "NULL"
out += f"""
{int_t(f"{vname}_p", mod, q_nlimbs)}
{int_t(f"{vname}_pinv", 1, q_nlimbs)}
static const unsigned int {vname}_s1_indices[{len(s1_indices)}] = {intlist2intarray(s1_indices)};
static const lin_params_t {vname} = {{{{ {name}, {deg}, {vname}_p, {vname}_pinv, {k}, {vname}_s1_indices, {len(s1_indices)}, {message_ptr}, {len(message_indices)}, {vname_Ps}, {matPs_nrows}, {strEs}, {strEs_nrows}, NULL, NULL, 1 }}}};
"""

printc(out)
