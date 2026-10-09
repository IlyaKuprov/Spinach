# examples/nmr_solids/redor_curve.m

- Signature: `redor_curve()`
- Source: [examples/nmr_solids/redor_curve.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/redor_curve.m)

## Purpose

Computes and plots a REDOR dephasing curve for a two-spin 13C–15N pair with Fokker–Planck MAS propagation. The source estimates a calculation time of seconds.

## Spin pair and rotor sampling

The two spins have zero scalar Zeeman entries, with coordinates `[0,0,0]` and `[0,0,1.47]`; the field parameter is 14.1. No units are stated for the field or separation. The full spherical-tensor Liouville basis is used without approximation. The simulation uses `singlerot` and the REDOR callback, rate 10000, rotor axis `[sqrt(2/3), 0, sqrt(1/3)]`, maximum rank 9, and grid `leb_2ang_rank_23`. No gradient is configured.

## REDOR evolution and observable

The initial state and receiver are both the 13C `Lx` state. The source evaluates cycle counts 0 through 48 and plots `real(curve(3,:)./curve(1,:))`, labelled as the normalised REDOR difference. The horizontal axis is labelled REDOR evolution time in rotor cycles. The example sets no explicit pulse-duration or pulse-amplitude parameters.
