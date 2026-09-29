# experiments/hyperpol/xixdnp_steady.m

- Signature: `dnp=xixdnp_steady(spin_system,parameters,H,R,K)`
- Context: call from a powder context; `H`, `R`, and `K` are supplied by that context.
- Canonical implementation: [`experiments/hyperpol/xixdnp_steady.m`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/hyperpol/xixdnp_steady.m)

## What it computes

This is the steady-state counterpart of the TPPM/XiX DNP protocol. For every microwave offset it builds the two-pulse contact block, repeats the block, appends the shot-spacing interval, finds the stroboscopic steady state, and returns the detection-state overlap. It is not a transient curve and does not take an initial density operator. As in the source comments, `R` must be thermalised to a finite temperature. Values are model outputs; this description asserts no measured polarisation or unrun result.

## Pulse block and units

The code forms `L=H+1i*R+1i*K`, then adds the microwave offset and centre shift `2*pi*(el_offs+addshift)*Ez`. `L1` is the +X electron drive; `L2` drives in the phase-set direction `Ex*cos(phase)+Ey*sin(phase)`. Both `el_offs`, `addshift`, and `irr_powers` are frequencies in Hz converted to angular-frequency terms with `2*pi`; `phase` is in radians. The block propagator is `propagator(L2,pulse_dur)*propagator(L1,pulse_dur)`, so it applies the +X pulse followed by the phase-set pulse, each with duration `pulse_dur` seconds. It repeats for `nloops`, then evolves under `L_curr` for `shot_spacing` seconds before the Newton steady-state solve. `nloops` is a positive integer.

## Inputs and output

Required fields are `irr_powers` (Hz), `coil` (detection state), `pulse_dur` (seconds), `phase` (radians), positive integer `nloops`, `shot_spacing` (seconds), scalar centre shift `addshift` (Hz), and offset vector `el_offs` (Hz). Compatible matrices `H`, `R`, and `K` come from the powder context. There is no `rho0` input: the state is found by the steady-state solver. The complex output `dnp` preserves the shape of `el_offs`; each entry is `coil'*rho_ss` for that offset.

## References

- Source paper: [DOI 10.1021/jacs.1c09900](https://doi.org/10.1021/jacs.1c09900).
- Spin Dynamics Wiki: [`xixdnp_steady.m`](https://spindynamics.org/wiki/index.php?title=xixdnp_steady.m).
- MATLAB implementation: [`experiments/hyperpol/xixdnp_steady.m`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/hyperpol/xixdnp_steady.m).
