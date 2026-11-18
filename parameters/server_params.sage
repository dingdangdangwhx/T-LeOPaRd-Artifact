import sys

if len(sys.argv) != 2 or sys.argv[1] not in ('partial72', 'plain72'):
    raise ValueError('usage: sage server_params.sage partial72|plain72')
profile = sys.argv[1]
name = '_params'
d = 64
fixed_q = 2**72 - 107
opening, rlen, keyrank, elen = 174, 45, 42, 87
m1, l = rlen + elen + 2, keyrank
nbin, nprime, Bprime = 0, 0, 0
Bopen = 2*sqrt(opening*d)
B0, B1 = (6072, 6073) if profile == 'partial72' else (12144, 12146)
n = [opening, elen, 1, 1]
B = [Bopen, 2156, B0, B1]
# The joint opening constraint also bounds its retained r subvector.
alpha = ceil(sqrt(sum(b*b for b in B)))
statistical_security = 130
verbose, code, loaded, codegen_err = 1, 1, 1, 0
load('lnp-tbox-codegen.sage')
if codegen_err or eknow >= mpf(2)**(-128):
    raise RuntimeError('LNP22 parameter screen failed')
scale = Integer(fixed_q).isqrt()
extraction_error = sqrt(2*rlen*d+1)*Bopen
assert 2*extraction_error < scale
assert 2*(scale+1)*extraction_error < fixed_q
print('#define NIZK2_IMPLICIT_OPENING 1')
print('#define NIZK2_KEY_BITS 0')
print('#define NIZK2_KEY_DIGITS 0')
print('static const lin_params_t params = {{ .tbox_params = _params, '
      '.ns1_indices = %d, .nm_indices = %d }};' % (m1, l))
