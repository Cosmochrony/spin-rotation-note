# Spinorial-Higgs transversality: exact verification script.
# Reproduction code for "Spin as a Transformation Class of Rotating Wave Modes" (v1.4),
# Section "Downstream boundary: chiral transversality and Higgs-typed composite channels".
#
# Scope, as corrected by Jerome on first review of the underlying reconnaissance note
# (fermionic-matter/front-spinorial-higgs-transversality-recon.md, working history, superseded here
# by the published paper as the citable source):
#   - Block 1 is a standard control (chiral 2x2 Hamiltonian, Bloch-vector precession, decoupling,
#     rephasing invariance). It reproduces known identities and is NOT the discriminating test.
#   - Block 2 is the actual open gate: the exact representation content of the composite condensate
#     psibar_R psi_L in each Yukawa sector (e, d, u) under Spin(3,1) x SU(3)_c x SU(2)_L x U(1)_Y,
#     including a genuine symbolic re-confirmation of the Lorentz-scalar claim already asserted (by
#     citation only) in the note, not a restatement by assertion.
#   - No photon/W/Z code: that mechanism is already-known physics illustrating section 3's principle,
#     not a test of the new intuition. Only the composite gate is the point of this script.
#
# Every check below is exact (sympy rationals / symbolic matrices), no floating-point sampling.

import sympy as sp

# ---------------------------------------------------------------------------
# Block 1 -- standard chiral 2x2 control
# ---------------------------------------------------------------------------

print("=" * 78)
print("BLOCK 1 -- standard chiral control (sanity check, not the discriminating test)")
print("=" * 78)

tau1 = sp.Matrix([[0, 1], [1, 0]])
tau2 = sp.Matrix([[0, -sp.I], [sp.I, 0]])
tau3 = sp.Matrix([[1, 0], [0, -1]])
I2 = sp.eye(2)
taus = [tau1, tau2, tau3]

p_abs, eta = sp.symbols("p_abs eta", real=True, positive=False)
mr, mi = sp.symbols("m_r m_i", real=True)  # Re(m), Im(m)
m = mr + sp.I * mi

# H_eta = -eta*p*tau3 + Re(m)*tau1 - Im(m)*tau2  (section 1 of the note)
H_eta = -eta * p_abs * tau3 + mr * tau1 - mi * tau2
B_eta = sp.Matrix([mr, -mi, -eta * p_abs])  # Bloch vector, exactly as specified

# (1) H_eta^2 = (|p|^2 + |m|^2) I, exactly, using eta^2 = 1
H2 = sp.expand(H_eta * H_eta)
H2 = H2.subs(eta ** 2, 1)
target = (p_abs ** 2 + mr ** 2 + mi ** 2) * I2
check_H2 = sp.simplify(H2 - target)
print("\n(1) H_eta^2 - (|p|^2+|m|^2) I  [should be the zero matrix]:")
print(check_H2)
assert check_H2 == sp.zeros(2, 2), "H_eta^2 identity FAILED"

# (2) eigenvalues +- sqrt(|p|^2+|m|^2)
E = sp.sqrt(p_abs ** 2 + mr ** 2 + mi ** 2)
eigs = H_eta.subs(eta ** 2, 1).eigenvals()
print("\n(2) Eigenvalues of H_eta (symbolic, eta treated as a sign not yet fixed):")
for val, mult in eigs.items():
    print("   ", sp.simplify(val), " (mult", mult, ")")
# Direct check against +-E after fixing eta = +1 (eta=-1 is the same statement, B_eta flips sign
# only in the tau3 component and the spectrum is eta-independent since it enters only as eta^2=1).
eigs_eta1 = sp.simplify(H_eta.subs(eta, 1)).eigenvals()
assert any(sp.simplify(v - E) == 0 for v in eigs_eta1.keys()), "missing +E eigenvalue"
assert any(sp.simplify(v + E) == 0 for v in eigs_eta1.keys()), "missing -E eigenvalue"
print("    Confirmed: {+E, -E} with E = sqrt(|p|^2+|m|^2).")

# (3) Chiral precession as an EXACT operator identity: i[H_eta, tau_i] = 2 (B_eta x tau)_i
print("\n(3) Precession identity  i*[H_eta, tau_i] == 2*(B_eta x tau)_i  for i=1,2,3:")


def cross_tau(B, taus_):
    Bx, By, Bz = B
    tx, ty, tz = taus_
    return [
        By * tz - Bz * ty,
        Bz * tx - Bx * tz,
        Bx * ty - By * tx,
    ]


rhs_list = cross_tau(B_eta, taus)
all_ok = True
for i in range(3):
    lhs = sp.I * (H_eta * taus[i] - taus[i] * H_eta)
    lhs = sp.expand(lhs).subs(eta ** 2, 1)
    rhs = sp.expand(2 * rhs_list[i])
    diff = sp.simplify(lhs - rhs)
    ok = diff == sp.zeros(2, 2)
    all_ok &= ok
    print(f"    i={i + 1}: difference = 0 ? {ok}")
assert all_ok, "precession identity FAILED"
print("    => d<s>/dt = 2 B_eta x s holds exactly (Heisenberg picture, d/dt tau_i = i[H,tau_i]).")

# (4) Decoupling at m=0: H_eta is diagonal, eigenstates are pure L, R.
H0 = H_eta.subs({mr: 0, mi: 0})
print("\n(4) H_eta at m=0 (should be diagonal):")
print(H0)
assert H0[0, 1] == 0 and H0[1, 0] == 0, "decoupling at m=0 FAILED"

# (5) Rephasing invariance: U^dagger H(m) U = H(m * e^{i(beta-alpha)}), U = diag(e^{i a}, e^{i b}).
alpha, beta = sp.symbols("alpha beta", real=True)
U = sp.diag(sp.exp(sp.I * alpha), sp.exp(sp.I * beta))
Udag = sp.diag(sp.exp(-sp.I * alpha), sp.exp(-sp.I * beta))
# H_eta = Re(m) tau1 - Im(m) tau2 - eta p tau3 has off-diagonal entries (upper-right, lower-left)
# = (m, conj(m)); verify this reconciliation exactly before using it.
H_m_generic = -eta * p_abs * tau3 + sp.Matrix([[0, m], [sp.conjugate(m), 0]])
check_form = sp.simplify(H_m_generic - H_eta.subs(eta ** 2, 1))
assert check_form == sp.zeros(2, 2), "tau-basis / (m, conj m) reconciliation FAILED"

lhs_rephase = sp.simplify(Udag * H_m_generic * U)
m_rotated = m * sp.exp(sp.I * (beta - alpha))
H_rotated = -eta * p_abs * tau3 + sp.Matrix([[0, m_rotated], [sp.conjugate(m_rotated), 0]])
diff_rephase = sp.simplify(sp.expand(lhs_rephase - H_rotated))
print("\n(5) U^dagger H(m) U - H(m * e^{i(beta-alpha)})  [should be zero]:")
print(diff_rephase)
assert diff_rephase == sp.zeros(2, 2), "rephasing identity FAILED"
print("    => arg(m) is the basis-dependent orientation of the transverse chiral component.")
print("       For a single Dirac sector it can be removed by relative chiral rephasing")
print("       (beta - alpha is free); only rephasing-invariant phases involving several")
print("       sectors (e.g. a CKM-type phase across Yukawa matrices) can be physical.")

print("\nBlock 1: all standard identities verified exactly.\n")

# ---------------------------------------------------------------------------
# Block 2 -- the composite-condensate gate (the actual discriminating test)
# ---------------------------------------------------------------------------

print("=" * 78)
print("BLOCK 2 -- exact representation audit of psibar_R psi_L per sector")
print("=" * 78)

# --- (a) Lorentz-scalar check, independent symbolic re-confirmation --------
#
# Two-component (Weyl) infinitesimal SL(2,C) transformations, generic rotation angle theta and
# boost rapidity beta (both real, treated to first order / as exact generator algebra):
#   delta psi_L = ( i*theta.sigma/2 + beta.sigma/2 ) psi_L      (left-handed: (1/2,0))
#   delta psi_R = ( i*theta.sigma/2 - beta.sigma/2 ) psi_R      (right-handed: (0,1/2))
# i.e. both chiralities share the SAME rotation generator (spin-1/2 under the compact SU(2)) and
# carry OPPOSITE-sign boost generators. This is the standard, textbook fact; verified here exactly
# rather than merely cited.
print("\n(a) Lorentz-scalar content of psibar_R psi_L (= psi_R^dagger psi_L, two-component form):")

sx = sp.Matrix([[0, 1], [1, 0]])
sy = sp.Matrix([[0, -sp.I], [sp.I, 0]])
sz = sp.Matrix([[1, 0], [0, -1]])
sigmas = [sx, sy, sz]

th = sp.symbols("theta_1 theta_2 theta_3", real=True)
bt = sp.symbols("beta_1 beta_2 beta_3", real=True)

gen_L = sp.zeros(2)
gen_R = sp.zeros(2)
for k in range(3):
    gen_L += sp.I * th[k] * sigmas[k] / 2 + bt[k] * sigmas[k] / 2
    gen_R += sp.I * th[k] * sigmas[k] / 2 - bt[k] * sigmas[k] / 2

# psi_R^dagger psi_L as an abstract bilinear: represent psi_L, psi_R by symbolic 2-vectors and
# compute delta(psi_R^dagger psi_L) = (delta psi_R)^dagger psi_L + psi_R^dagger (delta psi_L),
# and check the RESULT IS ZERO AS AN OPERATOR (i.e. for arbitrary psi_L, psi_R), by checking the
# operator identity  gen_R^dagger + gen_L == 0  (since delta(psi_R)^dagger = psi_R^dagger gen_R^dagger).
lorentz_scalar_defect = sp.simplify(sp.expand(gen_R.H + gen_L))
print("    gen_R^dagger + gen_L  [should be the zero matrix for ALL theta, beta]:")
print(lorentz_scalar_defect)
assert lorentz_scalar_defect == sp.zeros(2, 2), "Lorentz-scalar identity FAILED"
print("    => psi_R^dagger psi_L is Lorentz-invariant: rotation terms match (same-sign, both spin-1/2")
print("       under the compact SU(2)); boost terms cancel (opposite-sign chiral boost generators).")
print("    This independently re-confirms (does not merely restate) the corpus fact cited in section 2")
print("    of the companion note (front-graviweak-order-parameter-recon.md, section O.3).")

# --- (b) SU(3)_c color-singlet check for the quark channels -----------------
print("\n(b) SU(3)_c singlet content of the quark-antiquark color contraction:")

# Standard Gell-Mann matrices (generators T^a = lambda^a / 2).
lam = [
    sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]),
    sp.Matrix([[0, -sp.I, 0], [sp.I, 0, 0], [0, 0, 0]]),
    sp.Matrix([[1, 0, 0], [0, -1, 0], [0, 0, 0]]),
    sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]),
    sp.Matrix([[0, 0, -sp.I], [0, 0, 0], [sp.I, 0, 0]]),
    sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]),
    sp.Matrix([[0, 0, 0], [0, 0, -sp.I], [0, sp.I, 0]]),
    sp.Matrix([[1, 0, 0], [0, 1, 0], [0, 0, -2]]) / sp.sqrt(3),
]
T = [m_ / 2 for m_ in lam]

# The contraction dbar^i Q^i (sum over color index i) transforms as
#   delta(dbar^i Q^i) = -i*eps^a * dbar^j (T^a)_{ji} Q^i + i*eps^a * dbar^i (T^a)_{ij} Q^j
# (antifundamental minus, fundamental plus). This is exactly zero for every generator a because
# (T^a)_{ji} = (T^a)_{ij} (Gell-Mann matrices are Hermitian, real part symmetric where it matters
# after relabelling i<->j) -- verified here as an EXACT matrix identity: trace consistency via
# sum_{i,j} dbar^j (T^a)_{ji} Q^i - sum_{i,j} dbar^i (T^a)_{ij} Q^j == 0 for arbitrary dbar, Q,
# i.e. the operator identity  T^a - (T^a)^T == 0  is NOT what we need (T^a is Hermitian, not
# symmetric in general for the imaginary ones) -- the correct operator statement is:
#   delta(dbar Q) = dbar (i eps^a T^a)^T-part cancellation via T^a Hermitian: (T^a)^dagger = T^a.
# Concretely: writing the singlet contraction as a matrix trace, the color-singlet combinator is
# invariant iff  T^a^T + T^a == 0 is NOT required; instead the exact invariance to check is on the
# bilinear dbar^i Q^i = Tr(Qbar^T ... ), most directly verified via: for the antifundamental index
# transforming with -(T^a)^T (the standard antifundamental representation rule), invariance of the
# contracted scalar dbar^i Q^i requires  (T^a)^T == T^a^*, i.e. T^a Hermitian -- which holds by
# construction. We verify this defining Hermiticity property exactly for all eight generators, which
# is precisely the condition making the trace/contraction an exact singlet.
color_singlet_ok = True
for a, Ta in enumerate(T):
    defect = sp.simplify(Ta.H - Ta)  # Hermiticity: T^a^dagger - T^a == 0
    ok = defect == sp.zeros(3, 3)
    color_singlet_ok &= ok
print(f"    All 8 Gell-Mann generators exactly Hermitian (T^a^dagger = T^a): {color_singlet_ok}")
assert color_singlet_ok, "color generator Hermiticity FAILED"
# Direct explicit check on a symbolic color vector: delta(dbar^i Q^i) with antifundamental rule
# dbar -> dbar (I - i eps.T^a)   [row vector, i.e. dbar_j -> dbar_i (delta_ij - i eps^a (T^a)_ij)]
# and fundamental rule Q -> (I + i eps.T^a) Q, to first order in eps:
Qcol = sp.Matrix(sp.symbols("Q1 Q2 Q3"))
Dbar = sp.Matrix([sp.symbols("Db1 Db2 Db3")]).T  # column, will use as row via transpose below
eps_a = sp.symbols("eps_a")
all_generators_ok = True
for a_index, Ta in enumerate(T):
    delta_Q = sp.I * eps_a * Ta * Qcol
    delta_Dbar_row = -sp.I * eps_a * (Dbar.T * Ta)  # antifundamental: dbar -> dbar - i eps dbar.T^a
    delta_singlet = sp.simplify(sp.expand((delta_Dbar_row * Qcol)[0] + (Dbar.T * delta_Q)[0]))
    ok = delta_singlet == 0
    all_generators_ok &= ok
    print(f"    delta(Dbar^i Q^i) under generator T^{a_index + 1}  [should be 0]: {delta_singlet}")
assert all_generators_ok, "explicit color-singlet invariance check FAILED for some generator"
print("    => dbar^i Q^i (and ubar^i Q^i) is an exact SU(3)_c singlet: the color trace absorbs the")
print("       full 3bar (x) 3 = 1 (+) 8 decomposition down to its singlet piece.")

# --- (c) SU(2)_L doublet content + exact hypercharge arithmetic ------------
print("\n(c) SU(2)_L representation content and hypercharge (exact fractions):")

# Standard-Model hypercharges (Y with Q = T3 + Y convention), as sp.Rational, exact.
Y = {
    "L_L": sp.Rational(-1, 2),  # lepton doublet (nu_L, e_L)
    "e_R": sp.Rational(-1, 1),
    "Q_L": sp.Rational(1, 6),  # quark doublet (u_L, d_L)
    "u_R": sp.Rational(2, 3),
    "d_R": sp.Rational(-1, 3),
}

channels = {
    "ebar_R L_L": {"conj_field": "e_R", "other": "L_L", "color": False},
    "dbar_R Q_L": {"conj_field": "d_R", "other": "Q_L", "color": True},
    "ubar_R Q_L": {"conj_field": "u_R", "other": "Q_L", "color": True},
}

expected_type = {"ebar_R L_L": "H", "dbar_R Q_L": "H", "ubar_R Q_L": "Htilde"}

results = {}
for name, info in channels.items():
    Y_conj = -Y[info["conj_field"]]  # hypercharge flips sign under Dirac/charge conjugation
    Y_other = Y[info["other"]]
    Y_total = sp.nsimplify(Y_conj + Y_other)
    higgs_type = "H" if Y_total == sp.Rational(1, 2) else (
        "Htilde" if Y_total == sp.Rational(-1, 2) else "NEITHER"
    )
    # SU(2)_L content: singlet (x) doublet = doublet, always (exact, dimension count 1*2=2).
    su2_content = "doublet (1 (x) 2 = 2)"
    # electric charges of the two doublet slots, Q = T3 + Y
    Q_upper = sp.Rational(1, 2) + Y_total
    Q_lower = sp.Rational(-1, 2) + Y_total
    neutral_slot = "upper" if Q_upper == 0 else ("lower" if Q_lower == 0 else None)
    results[name] = dict(Y=Y_total, type=higgs_type, Q_upper=Q_upper, Q_lower=Q_lower,
                          neutral_slot=neutral_slot, color=info["color"])
    print(f"    {name:14s}  Y = {Y_total!s:6s}  type = {higgs_type:8s}"
          f"  Q(upper,lower) = ({Q_upper!s},{Q_lower!s})  neutral slot: {neutral_slot}"
          f"  color-singlet available: {info['color']}")
    assert higgs_type == expected_type[name], f"hypercharge/type mismatch for {name}"
    assert neutral_slot is not None, f"no neutral component found for {name}"

print("\n    Target table reproduced exactly:")
print("    channel        Y      type")
for name in channels:
    r = results[name]
    print(f"    {name:14s} {r['Y']!s:6s} {r['type']}")

# --- (d) Unbroken U(1)_em generator on the neutral VEV direction ------------
print("\n(d) Exact check that the neutral-VEV direction is annihilated by Q = T3 + Y:")

T3 = sp.diag(sp.Rational(1, 2), sp.Rational(-1, 2))
for name, r in results.items():
    Qgen = T3 + r["Y"] * sp.eye(2)
    v = sp.Matrix([1, 0]) if r["neutral_slot"] == "upper" else sp.Matrix([0, 1])
    Qv = sp.simplify(Qgen * v)
    annihilated = Qv == sp.zeros(2, 1)
    print(f"    {name:14s}  Q . v_neutral = {list(Qv)}  -> unbroken U(1)_em: {annihilated}")
    assert annihilated, f"U(1)_em NOT unbroken on the neutral slot of {name}"

print("\nBlock 2, checks (a)-(d): all exact identities verified.")

# ---------------------------------------------------------------------------
# The discriminating question
# ---------------------------------------------------------------------------

print("\n" + "=" * 78)
print("DISCRIMINATING QUESTION: single composite order, or three uncorrelated bilinears?")
print("=" * 78)
print(
"""
What Block 2 established, exactly:
  - each of the three channels (ebar_R L_L, dbar_R Q_L, ubar_R Q_L) IS separately a legal
    Lorentz scalar x (color singlet where applicable) x SU(2)_L doublet x correct hypercharge
    object, matching the target table component-by-component;
  - each channel has an exact electrically-neutral slot on which SU(2)_L x U(1)_Y breaks to the
    correct unbroken U(1)_em.

What Block 2 does NOT supply, and what no computation in this script can supply without further
input: a canonical map identifying the three channels as conjugate realizations of ONE composite
order. The Standard-Model gauge group alone provides no generator, symmetry, or relation connecting
the lepton sector (L_L, e_R) to either quark sector (Q_L, u_R, d_R): they are representation-
compatible, independently typed objects, not images of one another under any transformation in
Spin(3,1) x SU(3)_c x SU(2)_L x U(1)_Y.

VERDICT (typing only, no dynamics assumed):
  Branch B is REPRESENTATIONALLY POSSIBLE for all three sectors (typing does not close it).
  No canonical identifying mechanism between the three channels is found within the audited gauge
  content alone (a genuine gap, not a proof of impossibility -- a flavor symmetry, a GUT embedding,
  or some other external structure could still supply one, but none is derived here).

  This does NOT leave condensation dynamics as the sole remaining obstruction. Three open items
  remain, and must not be collapsed into one:
    1. the ambient doubling M_L ~= C^2 (Branch B's own precondition, conditional, not tested here);
    2. a canonical mechanism relating/selecting the three bilinears as a single common order (the
       gap this script exposes -- not supplied by the audited gauge group);
    3. the condensation dynamics itself and the selection of the electroweak scale.

  Conditional on item 1, this audit leaves item 3 (condensation dynamics) as the next obstruction to
  producing AT LEAST ONE composite Higgs channel. Producing ONE UNIVERSAL Higgs order shared by all
  three Yukawa sectors additionally requires item 2, which the audited gauge group does not supply.

  Clean result of this front:  spinorial chirality => Higgs-typed composite channels are allowed.
  NOT (yet, and not by anything shown here):  spinorial rotation => a unique Higgs condensate.
"""
)
