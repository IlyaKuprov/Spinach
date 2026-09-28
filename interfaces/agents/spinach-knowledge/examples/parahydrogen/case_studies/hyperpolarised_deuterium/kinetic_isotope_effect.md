# examples/parahydrogen/case_studies/hyperpolarised_deuterium/kinetic_isotope_effect.m

- Signature: `kinetic_isotope_effect()`

## Purpose

Simulates the evolution of spin-state populations during ortho-deuterium bubbling in the presence of a parahydrogenation catalyst, followed by a 45-degree deuterium pulse and acquisition.

## Physical / mathematical content

- The spin system contains four `2H` nuclei. The experimental chemical shifts are 4.55, 4.55, −13.5, and −16.5; the specified scalar couplings are 12.0 between spins 1 and 2 and 0.24 between spins 3 and 4.
- DFT-derived nuclear quadrupolar interaction tensors are specified for spins 3 and 4, along with their Cartesian coordinates; the two D2 spins have no coordinates specified.
- The chemical-kinetics model uses parts `[1 2]` and `[3 4]`, rates `[-1 5000; 1 -5000]`, and initial concentrations `[1 0]`.
- At a magnetic field of 7.05, the calculation uses the `sphten-liouv` formalism without basis approximation. Relaxation combines `redfield` and `t1_t2` models, with zero equilibrium, secular retention, correlation times of 1 and 400 ps, longitudinal rates `[0.04 0.04 0 0]`, and transverse rates `[8.00 8.00 0 0]`.

## Numerical / algorithmic content

- `deut_pair(spin_system,1,2)` supplies the singlet, triplet, quintet, and associated coherence operators. The initial state is `unit_state(spin_system)`; detection operators are made traceless.
- After the `nmr` assumptions are applied, the Hamiltonian, relaxation superoperator, and free-kinetics operator are assembled. Bubbling is modelled by `magpump` with a traceless pumped state comprising `S` and all five `Q` components and an assumed pumping rate of `1e-1` (the source notes that this rate needs refinement).
- The state evolves under bubbling for 7 seconds in 1,000 steps of 0.007 seconds, then under free kinetics for 30 seconds in 1,000 steps of 0.03 seconds. A 45-degree pulse about the deuterium y axis is applied to the resulting trajectory before the transition-coherence observables are projected out.

## Implementation structure

- The first plot shows selected singlet, triplet, and quintet state coefficients against time, with labels `$S$`, `$T_{\pm1}$`, `$T_{0}$`, `$Q_{\pm2}$`, `$Q_{\pm1}$`, and `$Q_{0}$`.
- The second plot shows the post-pulse coherences labelled `$T_{1} \rightarrow T_{0}$`, `$Q_{1} \rightarrow Q_{0}$`, and `$Q_{2} \rightarrow Q_{1}$`.