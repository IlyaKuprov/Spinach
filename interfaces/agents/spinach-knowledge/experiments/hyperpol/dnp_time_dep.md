# experiments/hyperpol/dnp_time_dep.m

- MATLAB implementation: [experiments/hyperpol/dnp_time_dep.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/hyperpol/dnp_time_dep.m)

- Signature: `answer=dnp_time_dep(spin_system,parameters,H,R,K)`

## Purpose and physical scope

This routine produces time-domain `coil`-projection trajectories under a fixed microwave drive and offset. It converts the supplied generator and parameters to Liouville space, adds the microwave and electron-offset terms to the Hamiltonian, and propagates the initial state under the Hamiltonian, relaxation and kinetics. Electron-nuclear hyperfine effects are represented only when included in the supplied spin-system model. This routine returns signal trajectories; it does not reconstruct an image or execute an ESEEM/ENDOR sequence.

## Inputs

- `parameters.mw_pwr`: microwave power in radians per second; `parameters.mw_off`: microwave-frequency offset from the free-electron frequency in radians per second.
- `parameters.rho0`: initial state; `parameters.coil`: detection-state vector or horizontal stack.
- `parameters.mw_oper`: microwave irradiation operator; `parameters.ez_oper`: electron `Lz` operator.
- `parameters.dt`: time step in seconds; `parameters.nsteps`: number of `evolution` steps.
- H, R and K: Hamiltonian, relaxation and kinetics matrices supplied by the context function.

## Propagation and output axes

The routine calls Liouville-space conversion, updates the Hamiltonian with the microwave and offset operators, and calls `evolution` in 'multichannel' mode. The returned matrix has detection channels as rows and trajectory samples as columns; time spacing is dt. Its values are projections of the evolving state onto the supplied `coil` states. After `sim2liouv` conversion, the propagator supports `sphten-liouv` and `zeeman-liouv`. A `zeeman-hilb` density-matrix input is also accepted: `sim2liouv` first converts its generators, states and operators to `zeeman-liouv` before the formalism guard. Unlike the two steady-state scan routines, its source requires a thermalised relaxation superoperator.

## Source-coded numerical example

`examples/dnp_liq/odnp_liquid_2.m` models liquid-phase Overhauser DNP after a perfect electron ESR inversion pulse. It sets `parameters.mw_pwr`=0, `parameters.mw_off`=0, `parameters.dt`=1e-6 seconds and `parameters.nsteps`=1e3. Its plot uses a 0-to-1000-microsecond axis with 1001 samples; the pulse prepares the input state before this propagator call. This is a simulation setup, not experimental data.

## Source and attribution

- Source: `experiments/hyperpol/dnp_time_dep.m`
- <https://spindynamics.org/wiki/index.php?title=dnp_time_dep.m>
- Source attribution: ilya.kuprov@weizmann.ac.il
