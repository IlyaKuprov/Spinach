# experiments/hyperpol/dnp_time_dep.m

- Signature: `answer=dnp_time_dep(spin_system,parameters,H,R,K)`

## Purpose

Propagate a time-dependent spin state under microwave irradiation and return its projection onto the requested detection coil state(s) at every time step.

## Method

After converting to the Liouville representation when needed, the function adds the microwave power term `mw_pwr*mw_oper` and the frequency-offset term `-mw_off*ez_oper` to `H`. It calls `evolution` with the combined generator `H + 1i*R + 1i*K`, initial state `rho0`, coil state(s), time step, and number of steps, using multichannel output. The accepted formalisms are `sphten-liouv` and `zeeman-liouv`.

## Parameters / inputs

- `parameters.mw_pwr` — microwave power, in rad/s.
- `parameters.mw_off` — microwave frequency offset from the free electron, in rad/s.
- `parameters.rho0` — thermal-equilibrium initial state.
- `parameters.coil` — coil state vector or horizontal stack of coil states.
- `parameters.mw_oper` — microwave-irradiation operator.
- `parameters.ez_oper` — electron `Lz` operator.
- `parameters.dt` — time step, in seconds.
- `parameters.nsteps` — number of time steps.
- `H` — Hamiltonian matrix supplied by the context function.
- `R` — relaxation superoperator supplied by the context function.
- `K` — kinetics superoperator supplied by the context function.

## Output

- `answer` — matrix of projections of the trajectory onto each supplied coil state at each time step.

## Note

The relaxation superoperator must be thermalised for the selected type of calculation.

- Source documentation: <https://spindynamics.org/wiki/index.php?title=dnp_time_dep.m>
