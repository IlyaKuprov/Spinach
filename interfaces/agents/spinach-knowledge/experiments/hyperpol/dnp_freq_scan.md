# experiments/hyperpol/dnp_freq_scan.m

- Signature: `dnp=dnp_freq_scan(spin_system,parameters,H,R,K)`

## Purpose

Microwave frequency scan steady-state DNP experiment. Returns the steady-state expectation values for the states specified in `parameters.coil` at each supplied microwave irradiation frequency.

## Physical / mathematical content

- The calculation combines the Hamiltonian `H`, relaxation superoperator `R`, kinetics superoperator `K`, and microwave irradiation operator to obtain a steady state at each frequency.
- The Fokker–Planck path represents microwave phase on a grid and couples its Fourier derivative to the spin dynamics. The Liouville-space path applies an electron-frequency offset through `parameters.ez_oper`.
- The relaxation superoperator must **not** be thermalized for this calculation (`inter.equilibrium='zero'`).

## Numerical / algorithmic content

- Supported methods are `fp-backs`, `fp-gmres`, `lvn-backs`, and `lvn-gmres`. The `backs` methods use MATLAB backslash; the `gmres` methods use GMRES. The Liouville-space GMRES path uses an incomplete-LU preconditioner.
- Frequencies are processed in a parallel loop. The Fokker–Planck path averages the solved state over the microwave phase grid before evaluating the observables.
- The function checks matrix dimensions and input consistency, damps the trace direction of `R`, and checks its conditioning before solving. Fokker–Planck methods require `labframe` assumptions; Liouville-space methods require `esr` assumptions.

## Parameters / inputs

- `parameters.mw_pwr` — microwave power, rad/s.
- `parameters.mw_frq` — row vector of microwave frequency offsets (rad/s) relative to the reference g-factor.
- `parameters.g_ref` — reference g-factor around which frequency offsets are specified.
- `parameters.rho0` — thermal equilibrium state.
- `parameters.coil` — coil state vector or a horizontal stack thereof.
- `parameters.mw_oper` — microwave irradiation operator.
- `parameters.ez_oper` — electron Lz operator; required for the Liouville-space methods.
- `parameters.method` — calculation method: `fp-backs`, `fp-gmres`, `lvn-backs`, or `lvn-gmres`.
- `parameters.nphases` — number of microwave phase grid points for the Fokker–Planck path.
- `H` — Hamiltonian matrix, received from the context function.
- `R` — relaxation superoperator, received from the context function.
- `K` — kinetics superoperator, received from the context function.

## Outputs

- `dnp` — an array of steady-state expectation values for the states specified in `parameters.coil` at each supplied microwave frequency; rows correspond to frequencies and columns to coil states.

## Reference

- <https://spindynamics.org/wiki/index.php?title=dnp_freq_scan.m>

## Authors

- ilya.kuprov@weizmann.ac.il
- alexander.karabanov@nottingham.ac.uk
- walter.kockenberger@nottingham.ac.uk
- mariagrazia.concilio@sjtu.edu.cn