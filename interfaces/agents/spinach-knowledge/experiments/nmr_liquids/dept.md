# experiments/nmr_liquids/dept.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_liquids/dept.m) · [Spinach Wiki: dept.m](https://spindynamics.org/wiki/index.php?title=dept.m)

- Signature: `fid=dept(spin_system,parameters,H,R,K)`

## Purpose

DEPT (distortionless enhancement by polarisation transfer) pulse sequence, citing [DOI 10.1016/0022-2364(82)90286-4](https://doi.org/10.1016/0022-2364(82)90286-4). It is a one-dimensional liquid-NMR experiment; the source documents multiplicity-dependent phase behaviour, not measured output from a run.

## Inputs

- `parameters.sweep`: positive scalar F1 sweep width in Hz.
- `parameters.npoints`: positive integer number of F1 points.
- `parameters.spins`: two different isotope names in a cell array ordered {F1, F2}; the header example is `{'13C','1H'}`.
- `parameters.J`: nonzero working scalar coupling in Hz; the source sets each J interval to `abs(1/(2*parameters.J))`.
- `parameters.beta`: finite real selection-pulse angle in radians.
- `H`, `R`, and `K`: Hamiltonian, relaxation, and kinetics matrices supplied by the context function; the source combines them as `H + 1i*R + 1i*K`. The implementation requires the `sphten-liouv` formalism.

## Sequence and detection

The source begins at isotropic thermal equilibrium, applies a 90-degree x pulse to spin 2, then evolves for a J interval. It applies a 90-degree x pulse to spin 1 and constructs two phase alternatives from x and y pi rotations on spin 2. After another J interval, it combines the alternatives using opposite spin-1 pi rotations, applies the `parameters.beta` pulse about y on spin 2, and performs a third J interval. It decouples spin 2 and detects on spin 1 using that spin's `L+` state. F1 acquisition uses timestep `1/parameters.sweep` and `parameters.npoints` points.

The source's explanatory note states that DEPT135 gives CH and CH3 signals opposite in phase to CH2; DEPT90 gives only CH signals; DEPT45 gives positive CH, CH2, and CH3 signals; and quaternary carbons do not appear. It recommends isotope dilution to generate carbon isotopomers; see `dilute.m`.

## Output and scope

- `fid`: one-dimensional free induction decay detected on spin 1.

The phase and multiplicity statements above are sequence notes in the source, not reported measurements or run-verified results for a particular system.
