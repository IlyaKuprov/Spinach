# experiments/nmr_solids/cn2d_sq.m

[Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_solids/cn2d_sq.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=cn2d_sq.m)

## Experiment and coherence pathway

This function implements the single-quantum version of the 13C-detected 14N-13C MAS two-dimensional correlation experiment described by Jarvis, Haies, Williamson, and Carravetta ([DOI: 10.1039/c3cp50787d](https://doi.org/10.1039/c3cp50787d)). It takes `parameters.rho0` as the initial state and `parameters.coil` as the detection state.

`parameters.spins` is a two-entry cell array with 14N first and 13C second. The code obtains `L+` operators for these nuclei, forms their x/y components, and embeds the operators in the Fokker-Planck spatial dimension with `kron(speye(parameters.spc_dim),...)`. It evolves the initial state under 14N RF at `rf_pwr` for `rf_dur`: the cosine component uses the 14N x operator and the sine component uses the y operator. It then selects 13C single-quantum coherence (`-1` or `+1`) and 14N single-quantum coherence (`-1` or `+1`). The indirect dimension is propagated as a trajectory, followed by a 13C x-pi pulse; the second dimension is observed against the supplied coil. The output separates the sine and cosine States-quadrature components.

The source labels this a MAS experiment, while its explicit spatial field is `spc_dim`; it does not accept a rotor rate or field orientation as a sequence parameter. The supplied `H`, `R`, `K`, and `spin_system` carry the modeled dynamics. The function body defines finite-duration 14N RF evolution, coherence selection, 13C refocusing, and detection.

## Inputs and units

Call signature: `fid=cn2d_sq(spin_system,parameters,H,R,K)`.

- `parameters.spins`: `{'14N','13C'}` order, with exactly two entries.
- `parameters.spc_dim`: Fokker-Planck spatial dimension used to lift spin operators.
- `parameters.sweep`: two sweep widths in Hz.
- `parameters.npoints`: numbers of points in the two dimensions.
- `parameters.rho0`: initial state.
- `parameters.coil`: detection state.
- `parameters.rf_pwr`: 14N RF power in Hz.
- `parameters.rf_dur`: 14N RF pulse duration in seconds.
- `H`, `R`, `K`: Hamiltonian, relaxation, and kinetics matrices supplied by the context function.

The consistency checks require `sphten-liouv` formalism and matrix-valued `H`, `R`, and `K` of the same dimension. Sweep and point-count vectors have two entries, with the spin list in the documented 14N-then-13C order; RF power and pulse duration are positive real scalars.

## Output

`fid.sin` and `fid.cos` are the two States-quadrature components. Each is generated over the two configured dimensions: the first sweep/point-count entry is the indirect coherence evolution and the second is 13C detection. The source propagates `npoints(i)-1` steps after the starting state in each dimension.
