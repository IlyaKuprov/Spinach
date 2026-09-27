# experiments/hyperpol/dnp_field_scan.m

- Signature: `dnp=dnp_field_scan(spin_system,parameters,H,R,K)`

## Purpose

Magnetic-field scan of a steady-state DNP experiment. Returns the steady-state expectation values of the states specified in `parameters.coil` at each supplied magnetic field offset.

## Physical / mathematical content

- The calculation uses `H` (Hamiltonian), `R` (relaxation superoperator), and `K` (kinetics superoperator). The relaxation superoperator must **not** be thermalized for this calculation.
- The thermal equilibrium state and relaxation superoperator are assumed unchanged across the sweep; **do not use this function for broad magnetic-field sweeps**.

## Numerical / algorithmic content

- Supports the `sphten-liouv` and `zeeman-liouv` formalisms. It damps the trace direction in `R`, checks that `R` is sufficiently nonsingular, and forms the Liouvillian from `H`, `R`, and `K` with microwave and frequency-offset terms.
- Computes `b=R*parameters.rho0`, then solves for the steady state at each field offset in a `parfor` loop. The `'backslash'` method uses MATLAB's linear solver; `'gmres'` uses ILU-preconditioned GMRES.

## Parameters / inputs

- `parameters.mw_pwr` — microwave power, Hz.
- `parameters.mw_frq` — microwave frequency offset from the free-electron frequency at the reference B0 field, Hz.
- `parameters.fields` — vector of magnetic-field offsets from the reference B0 field, Tesla.
- `parameters.rho0` — equilibrium state at the reference B0 field.
- `parameters.coil` — coil state vector or a horizontal stack of coil state vectors.
- `parameters.mw_oper` — microwave irradiation operator.
- `parameters.ez_oper` — electron Lz operator.
- `parameters.method` — `'backslash'` for MATLAB's linear equation solver or `'gmres'` for ILU-preconditioned GMRES.
- `H` — Hamiltonian matrix, received from the context function.
- `R` — relaxation superoperator, received from the context function.
- `K` — kinetics superoperator, received from the context function.

## Output

- `dnp` — array of steady-state expectation values for the states specified in `parameters.coil` at each supplied field.

## Citation and link

- ilya.kuprov@weizmann.ac.il
- <https://spindynamics.org/wiki/index.php?title=dnp_field_scan.m>