# experiments/esr_dipolar/eseem.m

- MATLAB implementation: [experiments/esr_dipolar/eseem.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_dipolar/eseem.m)

Source: https://spindynamics.org/wiki/index.php?title=eseem.m

`fid=eseem(spin_system,parameters,H,R,K)`

## What it calculates

This routine simulates an ideal hard-pulse ESEEM echo and returns its time-domain signal. In the physical experiment, ESEEM modulation can report nuclear-frequency structure coupled to the electron through hyperfine interactions. ENDOR is a distinct pulse experiment and is not implemented by this function. The modulation in a calculation depends on the caller-supplied spin system and Hamiltonian; the routine does not prescribe a hyperfine coupling value.

This is a single spin-dynamics sequence, not a field sweep, DNP or other hyperpolarisation workflow, or spatial imaging routine. No measured signal is read or claimed.

## Inputs

- `spin_system` — Spinach spin system. The routine accepts `sphten-liouv` and `zeeman-liouv` formalisms.
- `parameters.npoints` — number of sampled points; the source checks that it has one element.
- `parameters.timestep` — sampling time step in seconds; the source checks that it has one element.
- `parameters.rho0` — initial state.
- `parameters.coil` — detection state.
- `parameters.screen` — optional screen state, documented as the Hermitian conjugate of the detection state; omitted values default to `[]`.
- `parameters.pulse_op` — caller-supplied pulse operator, dimension-matched to `H`. Its flip-angle arguments below are in radians.
- `H`, `R`, `K` — dimension-matched numeric matrices supplied by the context: Hamiltonian, relaxation, and kinetics contributions. After `sim2liouv` conversion as needed, the routine forms `L=H+1i*R+1i*K`.

## Sequence and propagation

The supplied pulse operator acts on `rho0` with a `pi/2` rotation. The first evolution uses `timestep/2`, `npoints-1` steps, trajectory mode, and `screen`; a `pi` pulse is then applied, followed by a second `timestep/2` evolution over `npoints-1` steps in refocus mode with `coil` as observable. The code detects `full(coil'*rho_stack)` and transposes the result to form `fid`. The pulse operator is provided by the caller rather than constructed in this routine.

## Output and scope

`fid` is the time-domain echo signal, not a spectrum. The source does not return a separate time vector or name its array axes; `timestep` supplies the sampling interval. The source-backed numerical pulse examples are `pi/2` and `pi` radians, with two half-step evolution periods. It gives no concrete experimental parameter set, hyperfine value, or DOI.
