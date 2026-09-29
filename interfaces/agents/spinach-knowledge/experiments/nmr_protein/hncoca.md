# experiments/nmr_protein/hncoca.m

[Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_protein/hncoca.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=hncoca.m)

## What the sequence represents

This is a phase-sensitive three-dimensional HN(CO)CA sequence for 1H, 13C, and 15N-labelled proteins. The documented frequency dimensions are F1 = N, F2 = CA, and F3 = H. The source cites the reported experiment ([DOI: 10.1007/BF01874573](https://doi.org/10.1007/BF01874573)) and the bidirectional-propagation method ([DOI: 10.1016/j.jmr.2014.04.002](https://doi.org/10.1016/j.jmr.2014.04.002)).

This is the HN-to-N-to-carbonyl-to-CA correlation pathway, with the final proton channel detected. In the implementation, PDB atom labels select `H`, `N`, `C` (carbonyl), and `CA`; the source comment asks that labels such as `CA` and `HA` be provided through `sys.labels`. The code uses `spin_system.comp.labels` for the atom selections, and obtains isotope-wide proton and nitrogen pulse operators with `operator(...,'L+','1H')` and `operator(...,'L+','15N')`. Carbonyl and alpha-carbon operators are selected from the `C` and `CA` label masks. These are ideal Cartesian pulses applied with `step`.

The routine creates an NH `Lz` initial state and NH `L+` detection state when `parameters.rho0` or `parameters.coil` is not supplied. It selects positive and negative 15N coherence in F1, then constructs the four sign combinations by stitching forward state trajectories with backward-propagated detection trajectories. This is why the return value contains four States-quadrature components rather than one ordinary spectrum.

## Inputs and timing

Call signature: `fid=hncoca(spin_system,parameters,H,R,K)`.

- `parameters.npoints`: three integers ordered `[t1 t2 t3]`.
- `parameters.sweep`: three sweep widths in Hz, ordered `[f1 f2 f3]`.
- `parameters.tau`: four delays in seconds. The source gives `[2.25e-3, 2.75e-3, 8.00e-3, 7.00e-3]` as reasonable values.
- `parameters.rho0`: optional initial state; defaults to NH-proton `Lz` state.
- `parameters.coil`: optional detection state; defaults to NH-proton `L+` state.
- `H`, `R`, `K`: Hamiltonian, relaxation, and kinetics matrices supplied by the context function.

The routine requires the `sphten-liouv` formalism and matching matrix dimensions for `H`, `R`, and `K`. It checks for three entries in `npoints` and `sweep`, four in `tau`, and positive real delay values. It is hard-wired to the stated protein nuclei and atom-label conventions.

## Output

`fid` contains `pos_pos`, `pos_neg`, `neg_pos`, and `neg_neg` for States quadrature. The four 3D arrays are permuted with `[3 2 1]`, giving dimension order `[t3,t2,t1]` and nominal extents `[npoints(3),npoints(2),npoints(1)]`; these correspond to H, CA, and N, respectively.
