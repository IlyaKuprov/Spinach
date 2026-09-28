# experiments/nmr_liquids/mqs.m

- Signature: `fid=mqs(spin_system,parameters,H,R,K)`

## Purpose

Two-dimensional multiple-quantum NMR using the non-refocused multiple-quantum / MaxQ variant. Use `mqs_refocus.m` when post-mixing refocusing is required. The source cites [this paper](https://doi.org/10.1002/cphc.201800667) and [this communication](https://doi.org/10.1039/d1cc03079e).

## Sequence and signal

The function should be called through `liquid.m`, which supplies `H`, `R`, and `K`. It applies a 90-degree pulse, two evolution periods separated by a 180-degree pulse, and selects the requested multiple-quantum coherence order. The selected state is evolved in F1, mixed with the configured flip angle, and directly acquired in F2. The output is a 2D magnitude-mode FID.

## Syntax

```matlab
fid=mqs(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- `parameters.sweep`: [F1 F2] sweep widths, Hz.
- `parameters.npoints`: [F1 F2] numbers of points.
- `parameters.spins`: working spins, e.g. `{'1H','1H'}`.
- `parameters.angle`: flip angle, radians.
- `parameters.delay`: J-coupling evolution delay, seconds.
- `parameters.mqorder`: coherence order to select.
- `parameters.rho0`: initial state.
- `parameters.coil`: detection state.

## Output

- `fid`: 2D magnitude-mode free induction decay.

[Spin Dynamics Wiki page](https://spindynamics.org/wiki/index.php?title=mqs.m).
