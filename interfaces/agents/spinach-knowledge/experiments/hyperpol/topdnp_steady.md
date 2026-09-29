# experiments/hyperpol/topdnp_steady.m

- Signature: `dnp=topdnp_steady(spin_system,parameters,H,R,K)`
- Context: call from a powder context; `H`, `R`, and `K` are supplied by that context.
- Canonical implementation: [`experiments/hyperpol/topdnp_steady.m`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/hyperpol/topdnp_steady.m)

## What it computes

This is the stroboscopic steady-state implementation of the time-optimised pulsed DNP experiment cited in the source. It scans the microwave resonance offsets in `parameters.el_offs`, builds one repeated TOP block and its post-contact shot-spacing delay, solves for the periodic fixed point, and reports its overlap with the electron detection state `parameters.coil`. It models the spin-system dynamics represented by `H`, relaxation `R`, and kinetics `K`; the source requires finite-temperature thermalisation of `R`. This is a computed observable, not a claim of measured polarisation or a numerical result.

## Sequence and units

The source forms `L=H+1i*R+1i*K`. For each offset it adds `2*pi*(el_offs+addshift)*Ez` and electron microwave drive `2*pi*irr_powers*Ex`; the offset, centre shift, and microwave amplitude are specified in Hz and converted by `2*pi` to angular-frequency terms. The electron operators are `Ez` and `Ex`.

The block propagator is `propagator(L_curr,delay_dur)*propagator(L1,pulse_dur)`: with Spinach's left-acting propagators this applies the +X microwave pulse first, then the free/contact delay. That block is repeated `nloops` times, followed by evolution under `L_curr` for `shot_spacing`. The code finds the resulting stroboscopic steady state with `steady(...,'newton')`. `pulse_dur`, `delay_dur`, and `shot_spacing` are seconds; `nloops` is a positive integer.

## Inputs and output

Required experiment fields are `irr_powers` (microwave nutation frequency, Hz), `coil` (detection state), `pulse_dur` and `delay_dur` (seconds), `nloops`, `shot_spacing` (seconds), `addshift` (Hz), and numeric offset vector `el_offs` (Hz). The matrices `H`, `R`, and `K` must be dimension-compatible context matrices. The solver initialises one complex output element per offset; `size(dnp)` follows `size(parameters.el_offs)`. Each element is `coil'*rho_ss` at that offset. GPU use is optional when enabled by the Spinach system settings.

## References

- Source paper: [DOI 10.1126/sciadv.aav6909](https://doi.org/10.1126/sciadv.aav6909).
- Spin Dynamics Wiki: [`topdnp_steady.m`](https://spindynamics.org/wiki/index.php?title=topdnp_steady.m).
- MATLAB implementation: [`experiments/hyperpol/topdnp_steady.m`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/hyperpol/topdnp_steady.m).
