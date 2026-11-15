"""Parameters of the threshold scheme T-LeOPaRd.

The search is the one of ``t_leopard.parameters``, with these changes for the
threshold setting. Labels refer to the LaTeX sources of the threshold paper.

1. The modulus q is the smallest power of two that satisfies two conditions
   (``rule="both"``, the default). The first is the threshold correctness lemma
   (``lemma:correctness``, ``eq:qpratio-stat``):

       q/p >= 2^(kappa+2) hd (2 B_f + 1),

   with the Gaussian bound ``Bf-Gaussian`` on a coordinate of the final error,
   which sums the errors of W share components weighted by reconstruction
   coefficients of l1-norm at most eta:

       B_f = sqrt(ln(2) (kappa + 3 + log2(hd)) / pi) sqrt(W) eta (s1 + beta_r s sqrt((ell + m) d)).

   The second is p_u <= 2^-kappa for the uniqueness bound p_u of
   ``thm:uniqueness``, as in LeOPaRd. That bound uses the worst-case B_f of
   ``Bf-worstcase`` and a union bound over all quorums, and the request privacy
   bound of ``thm:request-priv`` contains the term 3 Q_RO_pk p_u. With
   ``rule="correctness"`` or ``rule="uniqueness"`` only one condition is imposed.

2. A share component is held by nu servers, which compute their pre-processing
   values and responses with the same matrix A_ep and independent errors. The
   reduction ``thm:MLWEPs-to-HMLWEPsI`` of the paper is the LeOPaRd reduction with
   s replaced by s / sqrt(nu) in its conditions, so for a fixed s0 the widths s
   and sbar1 grow by sqrt(nu). The PRF-side width sB does not grow: it is the
   nu = 1 value s / sqrt(nu), the smallest one the theorem allows with s_{0,B} = s0.
   The bound Q_max on the responses for one input is per server, so the nu holders
   of the key-carrying component together answer up to nu Q_max requests for one
   input, and the reduction of ``thm:pseudorandomness`` makes one more PRF query
   with the same matrix for the Eval oracle. Hence the theorem is used with
   Q_B = nu Q_max + 1 queries per matrix (``prf_queries_per_matrix``), and
   s1B = eta sqrt(nu Q_max + 1) sB grows by sqrt(nu) as well.

3. With more than one component, the shares of the short key are uniform modulo q,
   so their commitments use the encoding for large messages, two rows per entry.

4. The servers pre-process once per epoch and answer up to K requests per epoch with
   the same MLWE sample (``sessions``). The protocol relies on the K-use adaptive
   assumption nu-aMLWE-PRF-RU* (``sec:assumptions``), whose reduction multiplies the
   width sbar1 of the MLWE hint masks by sqrt(K) (``thm:MLWEPs-to-HMLWEPsI``, K identity
   hints per error) and by sqrt(1 + zeta) (``thm:HMLWEPsI-to-aHMLWEPs``), where zeta is
   the smallest value allowed by ``eq:rerand-cond`` for identity hint masks of width
   sbar' = sqrt(K) eta s. ``adaptive=False`` leaves out the factor sqrt(1 + zeta), which
   is the chain for masks fixed before the samples are returned (a variant that the
   paper does not state since the restructuring of its Section 3 on 2026-11-11).

5. The free parameters ``eta``, ``delta`` and ``eps_p`` of the theorem are passed
   through to the widths. The paper uses eta = delta = 1 and eps_p = 0.02 in all rows.

The NIZKs are left out of the sizes. Only the BP14 message mapping (h = 1) and
replicated secret sharing are supported.
"""

from dataclasses import dataclass

from sage.all import RR, Integer, RealField, binomial, ceil, log, pi, sqrt

from .hardness import smallest_secure_msis_rank
from .modulus import final_error_bound
from .proof_bounds import noise_budgets, accepted_final_error
from .parameters import BETA_R, MAX_ITERATIONS, MIN_OUTPUT_BITS, raise_rounding_modulus
from .sizes import DROPPED_BITS, kilobytes
from .widths import EPSILON, PAPER_EPS_P, gaussian_widths, max_queries_per_input, width_to_stddev

R256 = RealField(256)

MODULUS_RULES = ("both", "correctness", "uniqueness")


@dataclass(frozen=True)
class ReplicatedSharing:
    """Replicated (T, N) secret sharing (``fig:replicated-ss``).

    There is one component per set of N - T + 1 servers, held by each server in
    the set, and every component enters the reconstruction once with coefficient 1.
    T = N is additive sharing, and N = T = 1 is the single-server scheme LeOPaRd.
    """

    N: int
    T: int

    def __post_init__(self):
        if not 1 <= self.T <= self.N:
            raise ValueError(f"need 1 <= T <= N, got T = {self.T}, N = {self.N}")

    @property
    def components(self):
        return binomial(self.N, self.T - 1)

    @property
    def components_per_server(self):
        return binomial(self.N - 1, self.T - 1)

    @property
    def width(self):
        """Reconstruction width W, the same for every quorum."""
        return self.components

    @property
    def eta(self):
        """Bound on the l1-norm of the reconstruction coefficients."""
        return Integer(1)

    @property
    def multiplicity(self):
        """nu of ``lemma:ss-sim``, the number of servers that hold a component."""
        return self.N - self.T + 1

    @property
    def quorums(self):
        """Number of admissible quorums, |Q_{N,T}|."""
        return sum(binomial(self.N, k) for k in range(self.T, self.N + 1))

    @property
    def short_shares(self):
        """With a single component, the share is the short key itself."""
        return self.components == 1


@dataclass(frozen=True)
class ThresholdParameters:
    sharing: ReplicatedSharing
    mode: str
    kappa: int
    rule: str
    p: int
    d: int
    m: int
    ell: int
    h: int
    n_s: int
    n_c: int
    gadget_base: int
    s0: object
    s: object
    s1: object
    #: Responses per epoch and server, K.
    sessions: object
    #: Slack zeta of the rerandomization step, None without the adaptive factor.
    zeta: object
    log_q: int
    log2_p_u: object
    norm_model: str = "linf"
    B_s: object = None
    B_1: object = None
    error_digit_base: object = None


def prf_queries_per_matrix(sharing, q_max):
    """The bound Q_x^infty of ``thm:pseudorandomness`` on the PRF queries with one matrix B.

    ``q_max`` is the limit of one server on its responses for one input, Q_S in the paper
    (``rem:qmax-count``). The reduction makes one PRF query per response and secret-carrying
    share component. There are nu such components among the honest servers
    (``lemma:ss-sim``), so at most nu q_max queries have the same matrix B, and the Eval
    oracle adds one.
    """
    return sharing.multiplicity * q_max + 1


def gaussian_error_bound(kappa, d, m, ell, h, s, s1, beta_r, sharing):
    """B_f of ``Bf-Gaussian`` for errors of widths s and s1 in each component.

    A coordinate of the final error is Gaussian with parameter at most
    w = sqrt(W) eta (s1 + beta_r s sqrt((ell + m) d)). The two-sided tail bound
    Pr[|x| >= t] <= 2 exp(-pi t^2 / w^2) is at most p_e = 2^-(kappa + 2 + log2(hd)),
    as ``lemma:correctness`` requires, for t = sqrt(ln(2) (kappa + 3 + log2(hd)) / pi) w.
    """
    tail = RR(sqrt(log(2) * (kappa + 3 + log(h * d, 2)) / pi))
    return tail * RR(sqrt(sharing.width)) * sharing.eta * (s1 + beta_r * s * RR(sqrt((ell + m) * d)))


def adaptive_mask_slack(s, N, d, sessions, eta=Integer(1)):
    """Smallest zeta of ``eq:rerand-cond`` for identity hint masks of width sbar' = sqrt(K) eta s.

    The second condition reads sbar' >= sqrt(1 + 1/zeta) eta_eps(Z^(Nd)), with the bound
    eta_eps(Z^(Nd)) <= sqrt(ln(2Nd(1 + 1/eps))/pi) of ``lemma:smupbnd``.
    """
    s_prime = RR(sqrt(sessions)) * eta * s
    eta_eps = RR(sqrt(R256(log(2 * N * d * (1 + 1 / EPSILON))) / pi))
    if s_prime <= eta_eps:
        raise ValueError(f"identity hint masks of width {s_prime} are too narrow for the rerandomization")
    return RR(1 / ((s_prime / eta_eps) ** 2 - 1))


def worst_case_error_bound(d, m, ell, s, s1, beta_r, sharing, norm_model="linf", h=1):
    """B_f of ``Bf-worstcase`` with the NIZK bounds beta = 5 s and beta_1 = 5 s1."""
    if norm_model == "l2-digits":
        bounds=noise_budgets(s,s1,m+ell,h,d)
        return accepted_final_error(bounds.preprocessing,bounds.online,m+ell,d,
                                    beta_r,sharing.width,sharing.eta)
    if norm_model != "linf":
        raise ValueError("unknown norm model")
    return sharing.width * sharing.eta * final_error_bound(m, ell, d, s, s1, beta_r)


def correctness_modulus_bound(kappa, p, d, m, ell, h, s, s1, beta_r, sharing, norm_model="linf"):
    """Lower bound p 2^(kappa+2) hd (2 B_f + 1) on q from ``eq:qpratio-stat``, with B_f of ``Bf-Gaussian``."""
    bf = R256(gaussian_error_bound(kappa, d, m, ell, h, s, s1, beta_r, sharing))
    return p * R256(Integer(2) ** (kappa + 2)) * h * d * (2 * bf + 1)


def uniqueness_modulus_bound(kappa, p, d, m, ell, h, s, s1, beta_r, sharing, norm_model="linf"):
    """Lower bound on q for p_u <= 2^-kappa in ``thm:uniqueness``, with B_f of ``Bf-worstcase``."""
    bf = R256(worst_case_error_bound(d, m, ell, s, s1, beta_r, sharing, norm_model, h))
    return sharing.quorums * p * R256(Integer(2) ** kappa) * h * d * (2 * bf + 1)


def log_modulus(rule, kappa, p, d, m, ell, h, s, s1, beta_r, sharing, norm_model="linf"):
    """log2 of the smallest power-of-two modulus allowed by the chosen rule."""
    bounds = {
        "correctness": correctness_modulus_bound,
        "uniqueness": uniqueness_modulus_bound,
    }
    if rule == "both":
        q_bound = max(bound(kappa, p, d, m, ell, h, s, s1, beta_r, sharing, norm_model) for bound in bounds.values())
    elif rule in bounds:
        q_bound = bounds[rule](kappa, p, d, m, ell, h, s, s1, beta_r, sharing, norm_model)
    else:
        raise ValueError(f"unknown modulus rule {rule!r}, expected one of {MODULUS_RULES}")
    return ceil(log(q_bound, 2))


def log2_uniqueness_bound(log_q, p, d, m, ell, h, s, s1, beta_r, sharing, norm_model="linf"):
    """log2 of p_u = |Q_{N,T}| hd (2 B_f + 1) / (q/p), per query to RO_pk."""
    bf = R256(worst_case_error_bound(d, m, ell, s, s1, beta_r, sharing, norm_model, h))
    return log(sharing.quorums * h * d * (2 * bf + 1) * p, 2) - log_q


def threshold_widths(sharing, q_max, d, rows, h, *, sessions, eta=Integer(1), delta=Integer(1),
                     eps_p=PAPER_EPS_P, adaptive=True, scale_widths=True):
    """The widths (s0, s, s1) and the slack zeta for an MLWE matrix with ``rows`` = m + ell rows.

    ``q_max`` bounds the responses of one server for one input, and ``sessions`` is the
    bound K on the responses of a server per epoch. See items 2 and 4 of the module docstring.
    """
    width_factor = RR(sqrt(sharing.multiplicity)) if scale_widths else RR(1)
    q_B = prf_queries_per_matrix(sharing, q_max) if scale_widths else q_max
    widths = gaussian_widths(q_max, d, rows, h, delta=delta, eta=eta, eps_p=eps_p)
    s = width_factor * widths.s
    s1_bar = width_factor * widths.s1_bar
    s1_B = eta * RR(sqrt(q_B)) * widths.s  # s_B = widths.s, the nu = 1 value
    s1_bar = RR(sqrt(sessions)) * s1_bar  # K identity hints per error, thm:MLWEPs-to-HMLWEPsI
    zeta = adaptive_mask_slack(s, rows, d, sessions, eta=eta) if adaptive else None
    if adaptive:
        s1_bar = RR(sqrt(1 + zeta)) * s1_bar
    return widths.s0, s, RR(sqrt(s1_bar**2 + s1_B**2)), zeta


def find_threshold_parameters(sharing, mode, kappa, d, target_rhf, dimensions,
                              rule="both", scale_widths=True, beta_r=BETA_R, adaptive=True,
                              sessions=Integer(2) ** 10, eta=Integer(1), delta=Integer(1),
                              eps_p=PAPER_EPS_P, norm_model="linf"):
    """Parameters of T-LeOPaRd for one setting, by the fixed-point search of LeOPaRd.

    ``sessions`` is the bound K on the responses of a server per epoch, and ``eta``,
    ``delta`` and ``eps_p`` are the free parameters of the width theorem. The bound
    Q_max on the responses for one input is per server. ``scale_widths=False`` gives
    the LeOPaRd widths, with Q_max PRF queries per matrix and no factor sqrt(nu).
    """
    q_max = max_queries_per_input(mode, kappa)
    p = Integer(4) if d < 128 else Integer(2)
    h = Integer(1)
    m_guess = ell_guess = Integer(1)

    for _ in range(MAX_ITERATIONS):
        p = raise_rounding_modulus(p, d, h)
        s0, s, s1, zeta = threshold_widths(sharing, q_max, d, m_guess + ell_guess, h, sessions=sessions, eta=eta,
                                           delta=delta, eps_p=eps_p, adaptive=adaptive, scale_widths=scale_widths)
        log_q = log_modulus(rule, kappa, p, d, m_guess, ell_guess, h, s, s1, beta_r, sharing, norm_model)
        m = dimensions.rank(target_rhf, d, log_q, width_to_stddev(s0), uniform=False)
        ell = dimensions.rank(target_rhf, d, log_q, beta_r, uniform=True)
        if (m, ell) == (m_guess, ell_guess):
            break
        m_guess, ell_guess = m, ell
    else:
        raise RuntimeError(f"the search for m and ell did not converge for {sharing}, {mode}, kappa = {kappa}")

    # Binding of the commitments, with the MSIS norm bound of t_leopard.parameters.
    key_rows = m if sharing.short_shares else 2 * m
    n_c = smallest_secure_msis_rank(target_rhf, d, log_q, lambda n: 2 * RR(sqrt(n + h * (ell + m) + 1 + ell)))
    n_s = smallest_secure_msis_rank(target_rhf, d, log_q, lambda n: 2 * RR(sqrt(n + key_rows + m)))

    budgets=noise_budgets(s,s1,m+ell,h,d) if norm_model=="l2-digits" else None
    return ThresholdParameters(
        sharing=sharing, mode=mode, kappa=kappa, rule=rule, p=p, d=d, m=m, ell=ell, h=h,
        n_s=max(m, n_s), n_c=max(ell, n_c), gadget_base=ceil(Integer(2) ** (log_q / m)),
        s0=s0, s=s, s1=s1, sessions=sessions, zeta=zeta, log_q=log_q,
        log2_p_u=log2_uniqueness_bound(log_q, p, d, m, ell, h, s, s1, beta_r, sharing, norm_model),
        norm_model=norm_model, B_s=budgets.preprocessing if budgets else None,
        B_1=budgets.online if budgets else None, error_digit_base=budgets.digit_base if budgets else None,
    )


def check_threshold_parameters(params):
    """The sanity checks of ``t_leopard.parameters.check_parameters``. Raises ValueError if one fails."""
    pr = params
    if (pr.n_c + pr.h * (pr.ell + pr.m) + 1 + pr.ell) * pr.d >= RR(Integer(2) ** pr.log_q / 4):
        raise ValueError(f"{pr}: q is too small to extract from c_r")
    share_limit = RR(Integer(2) ** pr.log_q / 4) if pr.sharing.short_shares else RR(Integer(2) ** (pr.log_q / 2) / 4)
    if (pr.n_s + 2 * pr.m + pr.m) * pr.d >= share_limit:
        raise ValueError(f"{pr}: q is too small to extract from c_k,i")
    if RR(pr.d * pr.h * log(pr.p, 2)) < MIN_OUTPUT_BITS:
        raise ValueError(f"{pr}: 1/p^(dh) is not negligible")


@dataclass(frozen=True)
class ThresholdCosts:
    """Sizes in KB, without NIZKs. Server values are per server, client values are sent once."""

    public_key: object
    #: v_k,i for every component of the server, published once per epoch.
    server_preprocessing_per_epoch: object
    server_online: object
    #: c_r and C_x, sent with every request.
    client_online: object


def threshold_costs(params, security):
    """Sizes of the messages of ``fig:NLOPRF-setup-preproc`` and ``fig:NLOPRF``.

    The public key holds a commitment c_k,i per share component and the hash c_s
    of the seed (2 lambda bits). A server publishes v_k,i for each component it
    holds once per epoch and sends u_x,i for each of them per request. The client
    sends c_r and C_x with every request to all servers it contacts.
    """
    pr, sh = params, params.sharing
    key_rows = pr.m if sh.short_shares else 2 * pr.m
    c_k = pr.n_s * pr.d * (pr.log_q - DROPPED_BITS) + key_rows * pr.d * pr.log_q
    c_r = pr.n_c * pr.d * (pr.log_q - DROPPED_BITS) + (pr.h * (pr.ell + pr.m) + 1) * pr.d * pr.log_q
    return ThresholdCosts(
        public_key=kilobytes(sh.components * c_k + 2 * security),
        server_preprocessing_per_epoch=kilobytes(sh.components_per_server * (pr.ell + pr.m) * pr.d * pr.log_q),
        server_online=kilobytes(sh.components_per_server * pr.h * pr.d * pr.log_q),
        client_online=kilobytes(c_r + pr.h * pr.m * pr.d * pr.log_q),
    )
