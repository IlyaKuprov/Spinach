# kernel/contexts/powder.m

- Signature: `[answer,sph_grid]=powder(spin_system,pulse_sequence,...`

## Purpose

Static powder interface to pulse sequences. It generates a Liouvillian superoperator, initial state, and coil state, then passes them to the pulse-sequence function. It evaluates the sequence over orientations in a spherical averaging grid and returns either the weighted powder average or the result at each orientation. Syntax: `[answer,sph_grid]=powder(spin_system,pulse_sequence,... parameters,assumptions)`.

Use `singlerot` for MAS simulations.

## Parameters / inputs

- `pulse_sequence` — pulse-sequence function handle. See the experiments directory for pulse sequences shipped with Spinach.
- `parameters.spins` — cell array of active spins, in channel order, e.g. `{'1H','13C'}`. Required; spins must be present in the system.
- `parameters.offset` — numeric transmitter offsets in Hz, one per active spin. Defaults to zero offsets if omitted.
- `parameters.grid` — required name of the spherical averaging grid file in the kernel grids directory.
- `parameters.rframes` — rotating-frame specifications, each containing a spin and transformation order. For example, `{{'13C',2},{'14N',3}}` requests second-order carbon-13 and third-order nitrogen-14 transformations. Defaults to no additional rotating frames. When used, the assumptions on those spins should be laboratory frame.
- `parameters.decouple` — defaults to `{}` when omitted.
- `parameters.needs` — optional cell array requesting information for the sequence:
  - `'zeeman_op'` — lab-frame Zeeman operator, passed as `parameters.hzeeman`.
  - `'iso_eq'` — isotropic lab-frame thermal equilibrium, passed as `parameters.rho0`.
  - `'aniso_eq'` — orientation-specific thermal equilibrium from the full anisotropic lab-frame Hamiltonian, passed as `parameters.rho0`.
  - `'iso_eq'` and `'aniso_eq'` cannot be requested together. Neither can be requested when `parameters.rho0` is specified.
- `parameters.rho0` — initial state; may be a function handle of the three Euler angles in ZYZ active convention. A function handle is evaluated separately at each orientation.
- `parameters.serial` — if true, disables powder-grid parallelisation.
- `parameters.sum_up` — defaults to true. If false, returns the pulse-sequence output for each orientation rather than the powder average.
- `parameters.verbose` — defaults to `0`; when absent or zero, pulse-sequence output is silenced during the orientation loop.
- `parameters.*` — additional fields may be required by the pulse sequence; see its documentation page.
- `assumptions` — context-specific assumptions such as `'nmr'`, `'epr'`, or `'labframe'`; see the pulse-sequence header. Must be a character string.

## Processing

Before applying the user-supplied assumptions, the function conditionally builds the lab-frame Zeeman operator and lab-frame Hamiltonian components needed for requested equilibrium calculations. It then builds Hamiltonian components and the kinetics superoperator, and applies transmitter offsets to the isotropic Hamiltonian component. For each grid orientation it assembles the Hamiltonian, applies requested rotating-frame transformations, obtains the relaxation superoperator, and calls `pulse_sequence(spin_system,parameters,H,R,K)` with the orientation-specific parameters. Requested Zeeman operators and anisotropic equilibrium states are evaluated at each orientation.

By default, orientations are evaluated using MATLAB parallel processing; setting `parameters.serial` to true disables powder-grid parallelisation. The implementation uses MATLAB's Distributed Computing Toolbox to evaluate different system orientations in parallel on different labs.

Arbitrary-order rotating-frame transformations, including infinite order, are supported; see the header of `rotframe.m` for further information.

Source: https://spindynamics.org/wiki/index.php?title=powder.m

## Outputs

- `answer` — weighted powder average of the pulse-sequence outputs. If `parameters.sum_up` is false, this is a cell array containing the output at each orientation.
- `sph_grid` — grid structure containing the Euler angles (`alphas`, `betas`, `gammas`) and weights for each point.