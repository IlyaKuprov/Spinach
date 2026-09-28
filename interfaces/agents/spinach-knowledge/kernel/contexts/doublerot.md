# kernel/contexts/doublerot.m

- Signature: `[answer,sph_grid]=doublerot(spin_system,pulse_sequence,...`

## Purpose

Implements the double-angle-spinning context. In Liouville space it builds the Fokker–Planck evolution generator; in Hilbert space it builds a stack of spin Hamiltonians, one for each pair of rotor phases on the two-rotor phase grid. It passes the generator or Hamiltonian stack to the supplied pulse-sequence function handle.

## Parameters / inputs

- `pulse_sequence` — function handle for a pulse sequence in the experiments directory.
- `assumptions` — string passed to `assume.m` when the Hamiltonian is built.
- `parameters.rate_outer`, `parameters.rate_inner` — outer and inner rotor spinning rates in Hz.
- `parameters.axis_outer`, `parameters.axis_inner` — normalized three-element vectors specifying the outer and inner rotor axes.
- `parameters.rank_outer`, `parameters.rank_inner` — maximum harmonic ranks retained for the outer and inner rotors. Increase them until convergence; the source notes that the rank is approximately the number of spinning sidebands in the spectrum.
- `parameters.rframes` — rotating-frame specification; see the header of `rotframe.m`. When used, the assumptions for the affected spins should be in the laboratory frame.
- `parameters.grid` — spherical-grid file name from the kernel `grids` directory. Use a two-angle grid in Liouville space and a three-angle grid in Hilbert space.
- `parameters.needs` — cell array of sequence requirements. `'iso_eq'` requests the thermal-equilibrium state of the isotropic Hamiltonian in `parameters.rho0`.
- `parameters.serial` — when true, disables automatic parallelisation.
- `parameters.sum_up` — defaults to 1 and returns the powder average; set to 0 to return each orientation's result in a cell array.
- Other `parameters` subfields may be required by the pulse sequence; consult its documentation.

The wrapper sets `parameters.spc_dim` to the rotor-trajectory space dimension and `parameters.spn_dim` to the spin-dynamics matrix dimension before calling the sequence.

## Outputs

- `answer` — the weighted powder average, or, when `parameters.sum_up` is 0, the pulse-sequence outputs for individual orientations.
- `sph_grid` — the spherical grid used in the calculation.

## Notes

- Arbitrary-order rotating-frame transformations, including infinite order, are supported; see the header of `rotframe.m`.
- The state projector assumes a powder; single-crystal DOR is not currently supported.
- Parallel processing through MATLAB's Distributed Computing Toolbox is supported, with system orientations evaluated in parallel.