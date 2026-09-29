# experiments/nmr_liquids/dqf_cosy.m

- Signature: `fid=dqf_cosy(spin_system,parameters,H,R,K)`
- Source: [`experiments/nmr_liquids/dqf_cosy.m`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_liquids/dqf_cosy.m)

## Purpose and sequence

This phase-sensitive double-quantum-filtered COSY implementation starts from `Lz` magnetisation on the selected isotope, applies an `Lx` 90-degree pulse, and records an indirect-dimension (F1) trajectory. Two second-pulse branches, `Lx` and `Ly`, form the States quadrature components. An exact analytical coherence-order projection retains orders +2 and -2 in each branch; the code does not implement this filter as a phase cycle or a gradient-selection block. A third `Lx` 90-degree pulse precedes direct-dimension (F2) evolution and detection with the selected isotope's `L+` coil state. This describes the parameterised sequence, not a measured spectrum or run-verified output.

The Liouvillian is `L=H+1i*R+1i*K`; both dimensions use dwell time `1/parameters.sweep` seconds.

## Parameters and inputs

- `parameters.sweep`: positive real scalar sweep width in Hz, applied to both dimensions.
- `parameters.npoints`: two positive integer point counts, ordered F1 then F2.
- `parameters.spins`: one-element cell array naming an isotope present in the system (for example, `{'1H'}` or `{'13C'}`); the selected isotope must have at least two spins.
- `H`, `R`, and `K`: same-sized numeric Hamiltonian, relaxation, and kinetics matrices supplied by the context function. The function requires the `sphten-liouv` formalism.

## Outputs and references

- `fid.cos` and `fid.sin`: the two FID components for hypercomplex processing.
- [Double-quantum-filtered COSY reference, DOI 10.1016/0006-291X(83)91225-1](https://doi.org/10.1016/0006-291X(83)91225-1)
- [COSY reference, DOI 10.1021/ja00388a062](https://doi.org/10.1021/ja00388a062)
- [Spinach Wiki: `dqf_cosy.m`](https://spindynamics.org/wiki/index.php?title=dqf_cosy.m)
