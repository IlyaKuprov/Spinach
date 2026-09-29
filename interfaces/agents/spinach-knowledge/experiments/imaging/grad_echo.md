# experiments/imaging/grad_echo.m

- Signature: `fid=grad_echo(spin_system,parameters,H,R,K,G,F)`.
- Canonical implementation: `experiments/imaging/grad_echo.m` — https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/grad_echo.m.

## Contract and sequence

Call from `imaging()`, which supplies `H`, `R`, `K`, `G`, and `F`. The source forms `L=H+F+1i*R+1i*K`, obtains the spin operator for the first entry of `parameters.spins`, and applies a hard 90-degree pulse about `Ly` to the caller's initial state `parameters.rho0`. It evolves under the positive `G{1}` gradient for `g_n_steps` time steps, then acquires the echo under the opposite-sign gradient using the coil state as detector. This is a gradient-echo FID sequence; the returned data are not a reconstructed image or a k-space matrix.

The spin system must use `sphten-liouv` or `zeeman-liouv`; `H`, `R`, `K`, and `F` must be same-size matrices, and `G{1}` must be compatible with `H`. No DNP polarisation is generated or quantified here: `rho0` is an input.

## Parameters and units

- `parameters.rho0`: state vector matching the dimension of `H`; `parameters.coil`: numeric detection vector matching that dimension.
- `parameters.spins`: nonempty cell array of character strings; the first entry selects the spin operator.
- `parameters.npts`: positive spatial sample-count vector used to replicate the operator.
- `parameters.g_amp`: real scalar gradient amplitude in T/m, applied along `G{1}`.
- `parameters.g_step_dur`: time step in seconds; `parameters.g_n_steps`: positive integer number of dephasing steps.
- `G` must contain at least one gradient operator.

## FID dimensions and source-derived numerical facts

For `N=g_n_steps`, the first gradient evolution spans `N` steps. Observable evolution runs for `2*N` steps at spacing `g_step_dur`; the evolution routine includes the initial point, so a single coil yields a one-dimensional FID of `2*N+1` samples. The acquisition interval is `2*N*g_step_dur` seconds. For the smallest accepted count, `g_n_steps=1`, the returned trace has three samples at intervals of `g_step_dur` and spans `2*g_step_dur` seconds. This is a shape/timing example only, not a simulated or measured signal.

## References

- [Spinach documentation: `grad_echo.m`](https://spindynamics.org/wiki/index.php?title=grad_echo.m).
- [Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/grad_echo.m).
