# Spin as a Transformation Class of Rotating Wave Modes

J. Beau, Independent Researcher, France

## Status

Preprint, v1.2. Companion conceptual/structural audit of the fermionic-matter sub-programme.
DOI: [10.5281/zenodo.21380026](https://doi.org/10.5281/zenodo.21380026)
v1.1 adds the soldering-audit section: the Veronese obstruction lemma and the spinorial-soldering
indeterminacy proposition.
v1.2 is a rigour consolidation after review: surjectivity argument completed in the real-form
classification proof; the indeterminacy proposition's falsifiability clause restated exactly; the
three missing data qualified as independently missing, with the joint-sufficiency question left open.

## Abstract

This note examines the hypothesis that the wave function is a combination of rotating and
non-rotating wave modes, spin expressing the rotation of these components rather than the
rotation of a particle-like object.
For the photon the hypothesis is nearly literal: the transverse components of a circularly
polarized wave rotate, and the associated angular momentum is mechanically measurable.
A scalar counterexample shows that the rotation of components is not sufficient: spin is fixed
by the way the components transform under rotations of physical space.
The resulting hierarchy is: the relative phase is a measurable coordinate within a sector, the
weight m is an axial winding degree, the spin j labels the irreducible representation of SU(2),
and univalence (-1)^{2j} separates superselection sectors as a rule on the algebra of
observables.
The formulation avoids the historical superluminal-velocity objection, which constrains rigid
bodies rather than internal rotations of wave components.
Any rotational ontology of spin must involve a relational, extended, or topologically attached
wave configuration, never a point particle endowed with a mere internal axis.

## Epistemic status of the Cosmochrony-facing claims

1. **Exact (proved)**: the central grading on tensor powers of the admissible module,
   with the central element acting on V_rho^(tensor n) as (-1)^n.
2. **Diagnostic**: Lorentzian spin and the weak doublet must act on distinct tensor factors.
3. **Correction**: in Q7, e_0 is the spin-weight m = 0 within the j = 1 module, not a j = 0
   state (matching surgical fix applied to Q7.tex).
4. **Conditional**: the physical univalence interpretation of the tensor grading, conditional
   on the Lorentz-spin identification of V_rho.
5. **Open**: the representational problem (constructing the full Spin(3) action rho_rot with
   rho_rot(-I) = (-1)^n) and the stronger dynamical problem (whether the cascade canonically
   traces an axial one-parameter orbit, with axis and angle relation both derived).
6. **Proved (v1.1, unconditional)**: the Veronese-null versus compact-real-form lemma — nonzero
   pure tensors w w^T lie on the null cone of the invariant form, so they meet no compact real
   form W_phi; the first internal real-form selector candidate (O28 outer products) is excluded.
7. **Conditional on the listed corpus inventory (v1.1)**: the spinorial-soldering indeterminacy
   proposition — the derived data fix the complex carrier and its univalence but determine
   neither the spatial realisation, nor a physical rotation action, nor its scale; the missing
   soldering is typed as three data (phase-sensitive real-form selector, independent rotation
   action rho_sp, Casimir-independent normalisation), with rho_sp structurally blocking.

## Position in the programme

Upstream input: O29 / Q7 (admissibility construction selects V_rho isomorphic to C^2).
Downstream: Q14 and the fermionic-matter sub-programme (univalence, spinorial carrier,
Lorentz / weak SU(2) separation).

## Build

```
./compile.sh
```

Output: `out/SpinRotationNote.pdf`.
