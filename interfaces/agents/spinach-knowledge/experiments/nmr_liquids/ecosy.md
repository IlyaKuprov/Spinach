# experiments/nmr_liquids/ecosy.m

- Signature: `fid=ecosy(spin_system,parameters,H,R,K)`
- Source: [`experiments/nmr_liquids/ecosy.m`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_liquids/ecosy.m)

## Purpose and sequence

This phase-sensitive E.COSY implementation begins with an `Lx` initial (post-pulse) state on the selected isotope and records its indirect-dimension (F1) trajectory. The second pulse is represented by States quadrature branches using `Lx` and `Ly`. Each branch is projected onto coherence orders +/-2 through +/-6; the code weights these orders 1, 2, 4, 6, and 9, respectively. A third `Lx` pulse is applied to the cosine branch and `Ly` to the sine branch, followed by direct-dimension (F2) propagation and observation using the selected isotope's `L+` coil state. These are code-defined sequence operations, not a measured spectrum or a run-verified result.

The Liouvillian is `L=H+1i*R+1i*K`; both dimensions use dwell time `1/parameters.sweep` seconds.

## Parameters and inputs

- `parameters.sweep`: positive real scalar sweep width in Hz, for both dimensions.
- `parameters.npoints`: two positive integer point counts, ordered F1 then F2.
- `parameters.spins`: one-element cell array naming an isotope present in the system (for example, `{'1H'}` or `{'13C'}`).
- `H`, `R`, and `K`: same-sized numeric Hamiltonian, relaxation, and kinetics matrices supplied by the context function. The function requires the `sphten-liouv` formalism.

## Outputs and references

- `fid.cos` and `fid.sin`: real and imaginary States-quadrature FID components, as identified by the source header.
- [E.COSY reference, DOI 10.1021/ja00308a042](https://doi.org/10.1021/ja00308a042)
- [E.COSY reference, DOI 10.1063/1.451421](https://doi.org/10.1063/1.451421)
- [E.COSY reference, DOI 10.1016/0022-2364(87)90102-8](https://doi.org/10.1016/0022-2364(87)90102-8)
- [Spinach Wiki: `ecosy.m`](https://spindynamics.org/wiki/index.php?title=ecosy.m)
