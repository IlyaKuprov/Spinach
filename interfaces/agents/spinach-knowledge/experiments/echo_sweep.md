# experiments/echo_sweep.m

- MATLAB implementation: [experiments/echo_sweep.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/echo_sweep.m)

- Signature: `echo=echo_sweep(spin_system,parameters,H,~,~)`

## Purpose and pulse sequence

This is a two-pulse, echo-detected frequency sweep for pulsed EPR in Hilbert space. The carrier sweep samples the echo response across resonance offsets; for an anisotropic spin packet, orientation-dependent Zeeman and hyperfine terms shift the resonance during the sequence. The source describes static samples and magic-angle spinning (MAS), with a spinning P1 centre in diamond as its motivating case. The equal-duration finite pulses are separated by `tau` (from the end of pulse one to the start of pulse two); the echo is integrated over `echo_win` after pulse two. The carrier offset is stepped across the requested sweep. Pulse, delay, and echo-window durations are rounded to whole propagation steps.

## Inputs and operator construction

The caller supplies `rho0`, the density-matrix initial state, and `coil`, the detection state. `spins` is a one-element cell array naming the pulsed spin. The scalar timing/frequency settings are `pulse_dur` (seconds), `pulse_frq` (Hz), `tau` (seconds), `echo_win` (seconds), and `timestep` (seconds). `rate` is the spinning rate in Hz, with zero denoting a static sample; `nphases` is the number of start rotor phases averaged; `sweep` is the carrier sweep width in Hz; and `npoints` is the number of carrier offsets (an integer greater than two). The grumbler also requires `pulse_dur` and `echo_win` each to be at least `timestep`, and positive `nphases` not to exceed `spc_dim`. The caller supplies `spc_dim`, the rotor-stack size, and `H`, a cell array of Hamiltonian matrices indexed by rotor phase. The implementation constructs the selected spin's `Lx` and `Lz` operators for the pulse and carrier-offset terms. The rotor-stack elements must commute with the pulsed spin's `Lz`, because the offset is applied as a separate propagator; the source notes the `esr` assumption as a case satisfying this condition. It recommends matching rotor-stack resolution to the time step at the fastest rate, with context `max_rank` about `1/(2*abs(rate)*timestep)`. The source header calls the two ignored positional inputs `R` and `K`; the implementation signature discards them (`~,~`), and neither input is used.

## Returned spectrum and scope

`echo` is a complex column vector with `npoints` elements, one integrated signal per carrier offset. For each offset, the function sums the detected echo-window samples, multiplies by `timestep`, and averages over the requested start rotor phases. The source's Hamiltonian stack and rotor parameters define the spinning model; this function does not set up an instrument pulse program or supply a phase cycle. Its carrier-sweep output is not a time-domain acquisition trace. The [MAS diamond P1 example](../examples/esr_sol_pulsed/mas_diamond_p1.md) invokes this routine through `singlerot` with `@echo_sweep`.

[Source page](https://spindynamics.org/wiki/index.php?title=echo_sweep.m)
