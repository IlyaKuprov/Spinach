# kernel/contexts/powder.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/contexts/powder.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=powder.m)

## Contract

`powder(spin_system,pulse_sequence,parameters,assumptions)` is the static-powder interface; the source directs MAS calculations to `singlerot`. The pulse sequence must be a function handle, and the assumptions string is passed to `assume`. The context constructs the spin Hamiltonian in the chosen basis, the relaxation and kinetics operators, and any requested initial/detection information, then runs the sequence independently at each orientation.

There is no spatial-motion axis in this context: `parameters.spc_dim=1`, and `parameters.spn_dim=size(I,1)` for the Hamiltonian object `I` constructed in the selected spin basis. An orientation grid loaded from `parameters.grid` supplies `alphas`, `betas`, `gammas`, and `weights`. The three Euler angles rotate anisotropic spin-system contributions for each powder point; a function-handle initial state receives the angles in the ZYZ active convention. When `parameters.sum_up` is true, the outputs are summed with the grid weights; otherwise `answer` is a cell array with one pulse-sequence result per orientation. `sph_grid` returns the grid angles and weights.

## Parameters and assumptions

- `parameters.grid` is required and names a spherical grid file in `kernel/grids`.
- `parameters.spins` lists spin channels, for example `{'1H','13C'}`; matching `parameters.offset` entries are transmitter offsets in Hz. If offsets are omitted, they default to zero.
- `parameters.rho0` may be a spin state or a function handle of the three Euler angles. The source also supports the `iso_eq` and `aniso_eq` equilibrium requests; they cannot be combined with each other or with a supplied `rho0`.
- `parameters.needs` may request `iso_eq`, `aniso_eq`, or `zeeman_op`; the latter supplies the orientation-specific Hermitian Zeeman Hamiltonian as `localpar.hzeeman` to the sequence. The two equilibrium requests are mutually exclusive. Omitted `parameters.decouple` defaults to no decoupling.
- `parameters.rframes` can request rotating-frame transformations, for example `{{'13C',2},{'14N',3}}` specifies second order for carbon-13 and third order for nitrogen-14; the source header requires the respective spins to use laboratory-frame assumptions.
- The context evaluates orientations in parallel when available unless `parameters.serial` disables that parallelisation.

## Source-supported example

`examples/nqr/pure_nqr_iodine.m` runs a static powder NQR calculation for one 127I nucleus (spin 5/2), using a 560 MHz quadrupole interaction, asymmetry 0.01, the `rep_2ang_200pts_sph` grid, and 512 points. It requests orientation-dependent equilibrium with `parameters.needs={'aniso_eq'}`.

## State-dependent chemistry boundary

This context rejects a function handle returned by `kinetics` with `Spinach:powder:stateDependentKinetics`. Multi-reactant or callback-rate reaction records require a custom pulse sequence using `step`/`iserstep`, rather than static context assembly; see `examples/kinetics/nonlinear/bimolecular_closures.m` and `examples/microfluidics/reacting_flow_nmr.m`. Constant matrix kinetics remain supported.
