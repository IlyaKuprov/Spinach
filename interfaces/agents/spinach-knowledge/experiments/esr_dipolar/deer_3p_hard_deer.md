# experiments/esr_dipolar/deer_3p_hard_deer.m

- MATLAB implementation: [experiments/esr_dipolar/deer_3p_hard_deer.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_dipolar/deer_3p_hard_deer.m)

- Signature: `deer=deer_3p_hard_deer(spin_system,parameters,H,R,K)`

## Purpose and sequence

This implementation returns a three-pulse DEER trace: the pump pulse perturbs the probe-spin evolution, and the probe-coil signal is the measured trace. The probe excitation starts with a `pi/2` pulse; after the first evolution interval the pump operator applies a `pi` pulse, evolution is refocused across the corresponding interval, and a final probe `pi` pulse precedes the detection evolution. Each interval is represented by `parameters.stepsize` and `parameters.nsteps`. The function implements this three-pulse timing only; it does not describe a four-pulse DEER sequence.

## Inputs and pulse operators

`H`, `R`, and `K` must be same-sized matrices and form `H + 1i*R + 1i*K`. The parameter structure requires `rho0` (initial state), `coil_prob` (probe detection state), `stepsize` (time increment), `nsteps` (positive integer interval length), `ex_prob` and `ex_pump` (caller-supplied probe and pump excitation operators), and `output` set to `'brief'` or `'detailed'`. The pulse operators determine which spin or transition is affected; the routine does not construct the dipolar coupling, which must be represented in the supplied Hamiltonian if required by the model. The source does not state a unit for `stepsize`.

Hard pulses are appropriate only for spin-1/2 systems in the source's stated contract; higher-spin cases require transition-selective pulse operators. The function does not build those operators from spin labels.

## Output modes and detection

Both modes return a structure containing `deer_trace`, the probe-coil-detected trace. The trace is normalised by the coil norm; the Hilbert-space branch evaluates a trace against `coil_prob`, while the other branch uses the coil projection directly.

For `output='detailed'`, also provide `ex_hard`, `spectrum_sweep` (EPR sweep width in Hz), `spectrum_nsteps` (FID sample count), and `coil_pump`. The structure then additionally contains `hard_pulse_fid`, `prob_pulse_fid`, and `pump_pulse_fid`, the FIDs after the corresponding nonselective, probe-selective, and pump-selective `pi/2` pulses. Their sampling interval is `1/spectrum_sweep`.

[Source page](https://spindynamics.org/wiki/index.php?title=deer_3p_hard_deer.m)
