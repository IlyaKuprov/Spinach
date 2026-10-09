# experiments/nmr_liquids/deptq.m

- Signature: `fid=deptq(spin_system,parameters,H,R,K)`
- Source: [`experiments/nmr_liquids/deptq.m`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_liquids/deptq.m)

## Purpose and sequence

This is the fixed-first-proton-pulse, DEPTQ135-style variant; unlike `dept.m`, the source notes that quaternary-carbon signals appear. Starting from isotropic thermal equilibrium, it applies a carbon pulse, evolves under the working J coupling, and creates phase-alternated carbon/proton pulse branches. The branches are combined, followed by the adjustable proton editing pulse `parameters.beta`, another J interval, a carbon detection pulse, and proton decoupling. This is the implemented sequence design, not a claim of a measured spectrum or a run-verified result.

The Liouvillian is `L=H+1i*R+1i*K`; the J-evolution interval is `abs(1/(2*parameters.J))`. Detection uses the `L+` state for the first listed spin, and acquisition returns its one-dimensional FID. The dwell time is `1/parameters.sweep` seconds and acquisition requests `parameters.npoints` points.

## Parameters and inputs

- `parameters.sweep`: positive real scalar sweep width in Hz.
- `parameters.npoints`: positive integer point count.
- `parameters.spins`: two distinct isotope labels in a cell array. The first is the carbon channel used for carbon pulses and detection; the second is the proton channel used for proton pulses and decoupling. The source example is `{'13C','1H'}`.
- `parameters.J`: non-zero real working coupling in Hz.
- `parameters.beta`: finite real proton editing-pulse angle in radians.
- `H`, `R`, and `K`: same-sized numeric Hamiltonian, relaxation, and kinetics matrices supplied by the context function. The function requires the `sphten-liouv` formalism.

## Output and reference

- `fid`: free induction decay.
- The source recommends using `dilute.m` to generate carbon isotopomers.
- [DEPTQ paper, DOI 10.1006/jmre.1998.1595](https://doi.org/10.1006/jmre.1998.1595)
- [Spinach Wiki: `deptq.m`](https://spindynamics.org/wiki/index.php?title=deptq.m)
