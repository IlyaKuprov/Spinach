# experiments/imaging/grad_echo.m

- Signature: `fid=grad_echo(spin_system,parameters,H,R,K,G,F)`

## Purpose

Gradient echo pulse sequence. Call this function from the `imaging()` context, which supplies `H`, `R`, `K`, `G`, and `F`.

## Physical / mathematical content

The sequence assembles the Liouvillian `L=H+F+1i*R+1i*K`. It applies a hard 90-degree pulse to `parameters.rho0`, evolves the resulting state under an X gradient, then detects the echo under an X gradient of opposite sign.

## Numerical / algorithmic content

- Construct the pulse operator from `L+` for `parameters.spins{1}` and `parameters.npts`, then apply the pulse with `step`.
- Evolve under `L+parameters.g_amp*G{1}` for `parameters.g_n_steps` steps of duration `parameters.g_step_dur`, retaining the final state.
- Evolve under `L-parameters.g_amp*G{1}` for `2*parameters.g_n_steps` steps of the same duration, using `parameters.coil` to acquire the observable signal.

## Parameters / inputs

- `parameters.g_amp` — gradient amplitude in T/m; a real scalar.
- `parameters.g_step_dur` — gradient step duration; a positive real scalar.
- `parameters.g_n_steps` — number of gradient steps in the initial evolution; a positive integer.
- `parameters.spins` — nonempty cell array of character strings; its first entry selects the spin for the pulse operator.
- `parameters.npts` — positive integer used to construct the pulse operator.
- `parameters.rho0` — numeric initial state.
- `parameters.coil` — numeric detection operator.
- `H`, `R`, `K`, and `F` — numeric matrices of the same dimensions.
- `G` — cell array containing at least one gradient operator; the sequence uses `G{1}`.

The spin-system formalism must be `sphten-liouv` or `zeeman-liouv`.

## Outputs

- `fid` — time-domain echo signal.

## Reference

- <https://spindynamics.org/wiki/index.php?title=grad_echo.m>