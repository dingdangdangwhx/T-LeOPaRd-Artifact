from sage.all import sqrt

vname = "params"
deg = 64
M, N, H = 45, 42, 1
mod = 2**72 - 107
double_commitment = True
opening = 2*M + 2*N
key = opening
error = key + N
digits = error + M + N
dim = (M + 2*N + M + N + H, digits + 2*H)
message_indices = list(range(key, error))
wpart = [list(range(opening)), list(range(error, digits)),
         list(range(digits, digits+H)), list(range(digits+H, digits+2*H))]
wl2 = [2*sqrt(opening*deg), 2156, 6072, 6073]
wbin = [0, 0, 0, 0]
