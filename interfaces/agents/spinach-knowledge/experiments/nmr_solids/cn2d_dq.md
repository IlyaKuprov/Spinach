# experiments/nmr_solids/cn2d_dq.m

[Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_solids/cn2d_dq.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=cn2d_dq.m)

## Experiment and coherence pathway

This function implements the double-quantum version of a 13C-detected 14N-13C MAS two-dimensional correlation experiment described by Jarvis, Haies, Williamson, and Carravetta ([DOI: 10.1039/c3cp50787d](https://doi.org/10.1039/c3cp50787d)). It accepts the initial state and detection state as `parameters.rho0` and `parameters.coil`; it does not generate those states itself.

`parameters.spins` is a two-entry cell array with 14N first and 13C second. The routine obtains `L+` operators for those nuclei, converts them to x/y operators, and embeds them in the Fokker-Planck spatial dimension using `kron(speye(parameters.spc_dim),...)`. Preparation applies a finite 14N RF evolution of duration `rf_dur` at `rf_pwr`: the cosine component uses `Nx`, and the sine component uses `(Nx+Ny)/sqrt(2)`. The routine then selects 13C single-quantum coherence (`-1` or `+1`) together with 14N double-quantum coherence (`-2` or `+2`). The first sweep dimension is evolved as a trajectory, a 13C x-pi pulse is applied, and the second dimension is acquired against the supplied coil state. The returned sine and cosine components are the two States-quadrature components.

The source describes this as a MAS experiment, but this function's explicit spatial control is `spc_dim`; it has no sequence parameter for a rotor rate or field orientation. `H`, `R`, `K`, and `spin_system` are supplied by the context, so any modeled orientation or MAS dependence must be represented there. The function itself does not define a DANTE pulse train, REDOR rotor-synchronised recoupling block, overtone-CP step, or pseudocontact-tensor calculation; those interpretations are not supported by this mapped source.

## Inputs and units

Call signature: `fid=cn2d_dq(spin_system,parameters,H,R,K)`.

- `parameters.spins`: `{'14N','13C'}` order, with exactly two entries.
- `parameters.spc_dim`: Fokker-Planck spatial dimension used to lift spin operators.
- `parameters.sweep`: two sweep widths in Hz.
- `parameters.npoints`: numbers of points in the two dimensions.
- `parameters.rho0`: initial state.
- `parameters.coil`: detection state.
- `parameters.rf_pwr`: 14N RF power in Hz.
- `parameters.rf_dur`: 14N RF pulse duration in seconds.
- `H`, `R`, `K`: Hamiltonian, relaxation, and kinetics matrices supplied by the context function.

The consistency checks require `sphten-liouv` formalism and matrix-valued `H`, `R`, and `K` of the same dimension. The two-element sweep and point-count vectors and the 14N/13C spin ordering are part of this sequence's contract; RF power and pulse duration are positive real scalars.

## Output

`fid.sin` and `fid.cos` contain the sine and cosine States-quadrature data, respectively, over the two configured dimensions: the first sweep/point-count entry for the indirect coherence evolution and the second for 13C detection. The code advances each dimension with `npoints(i)-1` propagation steps after its initial state.
