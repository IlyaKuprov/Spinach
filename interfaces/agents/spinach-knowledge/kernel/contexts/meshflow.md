# kernel/contexts/meshflow.m

- Signature: `answer=meshflow(spin_system,pulse_sequence,parameters)`

## Purpose

First draft of the magnetohydrodynamics context for microfluidic simulations. It assembles evolution generators and passes them to a pulse-sequence function handle.

## Inputs

- `spin_system` must include mesh indexing and Voronoi tessellation information.
- `pulse_sequence` is a function handle. Pulse sequences that ship with Spinach are in the experiments directory.
- `parameters` must provide Hamiltonian, relaxation, and kinetics phantoms, plus initial-state and detection-state phantoms:
  - `H_ph` and `H_op`; `R_ph` and `R_op`; `K_ph` and `K_op`. Each phantom is a collection of spatial maps with the sample voxel-grid dimensions. `R_op` entries are relaxation superoperators.
  - `rho0_ph` and `rho0_st` specify spatial maps and corresponding spin states for the initial condition. `coil_ph` and `coil_st` specify detection maps and spin states; detection phantoms allow different voxels to be detected at different angles and sensitivities. Spin states are obtained from `state()`.
  - Additional `parameters` subfields may be required by the pulse sequence.

## Processing

Consistency checks are performed by `grumble` before generator construction. The function builds the Hamiltonian (`H`), relaxation superoperator (`R`), and kinetics superoperator (`K`) from their phantoms, and constructs the spatial flow/diffusion generator (`F`) and dummy gradient generator (`G`). It forms `parameters.rho0` and `parameters.coil` by summing Kronecker products of each spatial phantom map with its corresponding spin state.

The spatial dimension is `spc_dim`, the number of Voronoi cells (`spin_system.mesh.vor.ncells`); the spin dimension is `spn_dim`, the number of rows in `spin_system.bas.basis`. The total Fokker–Planck dimension is `spc_dim*spn_dim`. The function stores `spc_dim` and `spn_dim` in `parameters`. If polyadic support is disabled, `H`, `R`, `K`, `G`, and `F` are inflated.

The pulse sequence is called with `spin_system`, `parameters`, `H`, `R`, `K`, `G`, and `F`.

## Output

Returns whatever the pulse sequence returns.

Source: https://spindynamics.org/wiki/index.php?title=meshflow.m