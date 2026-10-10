# examples/nmr_liquids/cosy90_sucrose.m

- MATLAB implementation: [examples/nmr_liquids/cosy90_sucrose.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/cosy90_sucrose.m)

- Signature: `cosy90_sucrose()`

## Purpose

A liquid-state proton COSY-90 simulation for sucrose using magnetic parameters imported from a vacuum DFT calculation. The source estimates a calculation time of minutes.

## Spin system and basis

The example parses the sucrose DFT log at ../standard_systems/sucrose.log and calls `g2spinach` to import hydrogen nuclei as `1H`. It passes `options.min_j=2.0` (the helper defines this as a scalar-coupling threshold in Hz) and `options.no_xyz=1`; it supplies 31.8 ppm as the reference-shielding argument to `g2spinach`. The imported system is then assigned a field of 5.9 T. The Liouville-space basis uses IK-2, scalar-coupling connectivity, proximity level 1, and the greedy system-building option; the source also sets a proximity cutoff of 4.0. Zero track elimination is also explicitly enabled with `zte` in `sys.enable`.

## COSY acquisition and processing

The pulse angle is pi/2, the offset is 800 Hz, and the sweep width is 1700 Hz. The FID has 512 by 512 sampled points and is zero-filled to 2048 by 2048 before the 2D FFT. Both axes are in ppm. Cosine apodisation is applied on both time dimensions, and the plotted array is the real part of the shifted spectrum.

## Interpretation and scope

The plotted spectrum is generated from the magnetic parameters imported from the DFT log, rather than from shifts and couplings tabulated directly in this function. The source does not state the DFT method in the function body or present an experimental comparison, so the calculation should be read as a simulation using that supplied log and reference-shielding input, not as a measured sucrose spectrum.
