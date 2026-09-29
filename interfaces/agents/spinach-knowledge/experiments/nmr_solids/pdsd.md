# experiments/nmr_solids/pdsd.m

Source: [canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_solids/pdsd.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=pdsd.m)

## Purpose and inputs

A simplified model of the PDSD experiment with NOESY-type quadrature detection and a four-step phase cycle, called from the `singlerot` context. Signature: `fid=pdsd(spin_system,parameters,H,R,K)`. `H`, `R`, and `K` are square numeric matrices of equal size.

- `parameters.sweep` is a positive real sweep width in Hz; the indirect and direct dwell is `1/sweep`.
- `parameters.npoints` is a two-element vector of positive integers for indirect and direct samples.
- `parameters.tmix` is a non-negative mixing duration in seconds.
- `parameters.rate` is a real scalar in Hz, used as the amplitude of the proton irradiation term during mixing; `spc_dim` is a positive integer MAS spatial dimension.

## Sequence outline

The operators and states are hard-coded for `13C` and `1H`: the initial state is a spatially averaged `13C` `Ly` state, detection uses spatially averaged `13C` `L+`, and proton decoupling is requested for the indirect and direct evolution. For each of four phase-cycle pathways, the code evolves the indirect trajectory under the decoupled generator, applies one of the listed transverse `13C` rotations, evolves for `tmix` under `L+2*pi*rate*Hx` (`Hx` is the `1H` `Lx` control), then applies the pathway’s third rotation. It decouples `1H` again for detection evolution and subtracts paired pathway signals.

The outputs `fid.cos=fids{1}-fids{3}` and `fid.sin=fids{2}-fids{4}` are 2D quadrature components, sampled with `npoints(1)` indirect and `npoints(2)` direct points. This source initialises directly on `13C`; it contains no CP-transfer block. It is a simplified sequence model, not a runtime or experimental-validation claim.