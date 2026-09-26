# examples/optimal_control/case_studies/Smelko_ChemRxiv_2026/mqmas_drifts.m

- Signature: `drifts=mqmas_drifts(spin_system,parameters)`

## Purpose

Drift Hamiltonians of a quadrupolar nucleus under magic angle spinning, resolved in rotor phase, for every combination of a two-angle powder grid orientation and an initial rotor phase. The Hamiltonian is taken to second order in the rotating frame of the nucleus, and is held constant within each rotor phase tick. Syntax: drifts=mqmas_drifts(spin_system,parameters) Parameters: parameters.spins - the nucleus, a cell array with one isotope string, e.g. {'27Al'} parameters.axis - spinning axis, a normalised row vector with three elements parameters.grid - two-angle powder grid name; the grid must have uniform weights because the ensemble average in optimcon is unweighted parameters.n_ticks - rotor phase ticks per rotor period parameters.n_phases - number of initial rotor phases, must be a divisor of parameters.n_ticks parameters.n_slices - number of ticks in the pulse Outputs: drifts - cell array over the ensemble, grid orientations in the outer index and initial rotor phases in the inner index; each element is a cell array of n_slices Hamiltonian matrices, one per tick of the pulse Note: the crystallite orientation uses all three Euler angles of the grid, as in singlerot.m; two-angle grids keep the azimuth in the third angle and have zero first angles, which the rotor phase then supplies. The spin system must carry laboratory frame assumptions, call assume() first.

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Parameters / inputs

- parameters.spins -the nucleus, a cell array with one
- isotope string, e.g. {'27Al'}
- parameters.axis -spinning axis, a normalised row
- vector with three elements
- parameters.grid -two-angle powder grid name; the grid
- must have uniform weights because the
- ensemble average in optimcon is unweighted
- parameters.n_ticks -rotor phase ticks per rotor period
- parameters.n_phases -number of initial rotor phases, must
- be a divisor of parameters.n_ticks
- parameters.n_slices -number of ticks in the pulse

## Outputs

- drifts -cell array over the ensemble, grid orientations in
- the outer index and initial rotor phases in the in-
- ner index; each element is a cell array of n_slices
- Hamiltonian matrices, one per tick of the pulse
- Note: the crystallite orientation uses all three Euler angles of the
- grid, as in singlerot.m; two-angle grids keep the azimuth in
- the third angle and have zero first angles, which the rotor
- phase then supplies. The spin system must carry laboratory
- frame assumptions, call assume() first.

## Implementation structure

- Drift Hamiltonians of a quadrupolar nucleus under magic angle spin-
- ning, resolved in rotor phase, for every combination of a two-angle
- powder grid orientation and an initial rotor phase. The Hamiltonian
- is taken to second order in the rotating frame of the nucleus, and
- is held constant within each rotor phase tick. Syntax:
- drifts=mqmas_drifts(spin_system,parameters)
- parameters.spins -the nucleus, a cell array with one
- isotope string, e.g. {'27Al'}
- parameters.axis -spinning axis, a normalised row
- vector with three elements
- parameters.grid -two-angle powder grid name; the grid
- must have uniform weights because the
