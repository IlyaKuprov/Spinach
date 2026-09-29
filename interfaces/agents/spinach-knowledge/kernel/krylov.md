# kernel/krylov.m

- Signature: `answer=krylov(spin_system,L,coil,rho,timestep,nsteps,output)`
- Direct MATLAB source: [`kernel/krylov.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/krylov.m)
- Existing Wiki: [`krylov.m`](https://spindynamics.org/wiki/index.php?title=krylov.m)

## Purpose and propagation

Propagates state vectors without constructing the full propagator; the source recommends it when the propagator does not fit in memory but `L` does, while warning that it may be slow. A propagation step applies the matrix-exponential action to the state (Spinach convention: `rho_next=exp(-1i*L*timestep)*rho`); this implementation advances states through calls to `step(spin_system,L,rho,timestep)` ([helper source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/step.m)). The source comment identifies the implemented approach as a reordered Taylor process rather than a Krylov-subspace/Arnoldi iteration. In `zeeman-wavef`, the source header identifies `L` as the Hamiltonian and `rho` as a wavefunction; `zeeman-hilb` is unsupported.

The routine has no frequency-unit conversion. `timestep` must be consistent with the units/convention of `L`; the source does not declare a standalone Hz-versus-angular-frequency conversion or normalisation of the input/output state.

## Inputs and guards

The formalism guard accepts `sphten-liouv`, `zeeman-liouv`, and `zeeman-wavef`. The source checks that `L`, `coil`, `timestep`, and `nsteps` are numeric, that `rho` is numeric or a cell array, and that `output` is a character value in a whitelist. These are type/formalism checks; this function does not add explicit dimensional, finiteness, sign, or integer checks in the visible guards. GPU-enabled systems move `L`, `rho`, and `coil` to the GPU; otherwise `rho` is made full on the CPU.

## Output modes and shapes

- `final`: advances once by `timestep*nsteps`; returns the gathered final state, with the state dimensions of `rho`.
- `trajectory`: returns state vectors at the initial point and after each step, allocated as `size(rho,1)×(nsteps+1)`.
- `refocus`: sets the step count to `size(rho,2)`; successively advances columns 2:end, 3:end, and so on, then returns the resulting `rho` matrix. This schedules the initial vector at zero steps, the next at one step, etc., as used for the second indirect-dimension evolution after a refocusing pulse.
- `observable`: records `coil'*rho` initially and after every step, returning `(nsteps+1)×size(rho,2)` for a single detection vector.
- `multichannel`: records all coil overlaps, returning `size(coil,2)×(nsteps+1)×size(rho,2)`; the source header notes that destination-state screening can be less efficient with multiple destinations.

The input whitelist also contains `total`, but the dispatch switch has no `total` case; any unmatched dispatch reaches the source's `unknown output option` error. This distinction is based on the source's guard and dispatch, not a runtime test.
