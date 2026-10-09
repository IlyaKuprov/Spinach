# kernel/relaxation.m

- Signature: `R=relaxation(spin_system,euler_angles)`

## Purpose

Builds a relaxation superoperator by accumulating the terms selected in `spin_system.rlx.theories`. It can also add dissipative bosonic-mode terms from `spin_system.inter.modes` in Liouville-space formalisms.

## Physical / mathematical content

- The source dispatches among several named models, including extended `t1_t2`, Redfield and Nakajima–Zwanzig rotational-diffusion treatments, a user-rate Lindblad model, and other named relaxation models. A returned `R` is a dissipative superoperator; when assembling a Liouvillian manually, retain the documented convention `L=H+1i*R+1i*K`, not `H+R+K`.
- `euler_angles` specifies orientation in the active ZYZ convention, in radians, relative to the input orientation. It affects theories that use relaxation-rate anisotropy; the header explicitly notes that Redfield is orientation-independent here. Mode dissipators are independent of these spin Euler angles.
- The Redfield and Nakajima–Zwanzig choices are alternative evaluations of the same kernel and are rejected when selected together. The source also enforces model/formalism combinations: for example, extended `t1_t2` is restricted to `sphten-liouv`, while diagonal retention is rejected for `zeeman-liouv`.
- For SRSK, the source recursively builds an extended `t1_t2` contribution with the local equilibrium set to `'zero'`; it leaves thermalisation and mode dissipation to the outer call. The outer call applies the selected thermalisation to the accumulated spin relaxation and appends bosonic-mode terms once. The source's SRSK rate-report text labels the displayed `R1` and `R2` increments in Hz.

## Numerical / algorithmic content

The function starts with a zero superoperator and adds the enabled spin-relaxation contributions. For Redfield, it checks that rotational correlation times are nonzero and checks the stated `T1,2 >> tau_c` condition against the constructed relaxation rates. Positive bosonic-mode damping or dephasing adds thermalised GKSL terms in Liouville space; outside the supported Liouville formalisms the source reports that those mode terms are not added. When no spin theory is selected and no dissipative mode term is present, the result is set to zero.

Retention by longitudinal order or base frequency uses each substance’s local descriptor. Diagonal retention and uniform damping exempt the unit coordinate of every block separately.

Nottingham's four-level electron manifold is implemented for a single substance with exactly two electrons. `relaxation` rejects every segmented Nottingham descriptor with `Spinach:relaxation:nottinghamSubstance`, including separate two-electron substances; nucleus-only, spin-free, and split-electron partners are not assigned partial models. The existing `create` restriction of two electrons overall is unchanged.

## Parameters / inputs

- `spin_system` — Spinach system containing the selected relaxation theories and their parameters.
- `euler_angles` — optional orientation specification for theories that support anisotropic rates; active ZYZ angles in radians.

The function header does not state units for the relaxation-rate or correlation-time fields, so this entry does not assign them. Its explicit unit note is that the Euler angles are in radians.

## Outputs

- `R` — relaxation superoperator. Spinach contexts include relaxation and kinetics superoperators in the total Liouvillian automatically; use the `1i*R` convention when assembling one manually.

## Source

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/relaxation.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=relaxation.m)

SRSK normalised source vectors use unweighted `coil_state`, so zero concentration never makes their norm vanish. Diagonal retention and uniform damping use geometric units independent of concentration. IME requests unit-concentration equilibrium shapes before adding unit-column sources; propagated populations provide the concentration weighting.

Normalised SRSK vectors explicitly request the `exact` method of the four-argument unweighted `coil_state` primitive.

## Explicit NZ evaluation point

`inter.nz_shift` must be an explicit finite scalar with non-negative real part. The former `'chem'` request is rejected by `create`: a general reaction network does not specify a unique scalar lifetime. For a legacy exponential radical-pair model use the summed channel rates; for the legacy Haberkorn/Jones–Hore scalar approximation use half their sum, stating that approximation explicitly in the example. `relaxation` uses the supplied scalar without reading retired radical-pair fields.
