# examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/mqmas_drifts.m

- Signature: `drifts=mqmas_drifts(spin_system,parameters)`
- Status: historical documentation. The corresponding `.m` file is absent from the current checkout; this description is based on `3975f139^:examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/mqmas_drifts.m` and does not establish current behavior.

## Purpose and output

The historical function builds drift Hamiltonians for a quadrupolar nucleus under magic-angle spinning, resolved for each powder-grid orientation and each initial rotor phase. It uses the second-order rotating-frame Hamiltonian and holds each Hamiltonian constant over its rotor-phase tick.

The output `drifts` is a cell array over the ensemble: grid orientations are in the outer index and initial rotor phases in the inner index. Each element is a cell array of `parameters.n_slices` Hamiltonian matrices, one per pulse tick.

## Inputs and constraints

- `parameters.spins`: a one-element cell array containing an isotope string present in the spin system (the source gives `{'27Al'}` as an example).
- `parameters.axis`: normalised three-element row vector specifying the spinning axis.
- `parameters.grid`: two-angle powder-grid name. Its weights must be uniform; the historical implementation checks deviations from the first weight against (10^{-12}), because the ensemble average is unweighted.
- `parameters.n_ticks`: positive integer number of rotor-phase ticks per rotor period.
- `parameters.n_phases`: number of initial rotor phases; it must divide `n_ticks`.
- `parameters.n_slices`: positive integer number of rotor ticks in the pulse.

The caller must supply a spin system with laboratory-frame assumptions already set (the source says to call `assume()` first). The grid is loaded from Spinach's `kernel/grids` directory.

## Historical construction

Tick phases are sampled at midpoints, (2pi(k-1/2)/n_{ticks}). The source composes the crystallite, rotor-axis, and rotor-phase Wigner rotations, adds the rotated anisotropic tensor components to the isotropic Hamiltonian, symmetrises the result, and applies the second-order rotating-frame transformation. For a given initial rotor phase, the slice sequence is selected by a cyclic shift of `n_ticks/n_phases` ticks. Orientations are evaluated with `parfor`.

For two-angle grids, the source uses all three Euler angles: the azimuth is in the third angle, the first two are zero, and the rotor phase supplies the remaining rotation. The function's path places it in the `Smelko_ChemRxiv_2026` case-study directory, but the historical source itself supplies no article title, DOI, or bibliographic citation. It ends with the attribution “Everything should be made as simple as possible, but not simpler.” — Albert Einstein.
