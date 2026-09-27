# examples/optimal_control/bloch_siegert/ramsey_shifts.m

- Signature: `ramsey_shifts()`

## Purpose

The current example lives at `examples/fundamentals/ramsey_shifts.m`; this knowledge entry retains its earlier path for lookup. A constant proton-channel drive shifts the off-resonant 13C and 15N Zeeman frequencies in a three-spin system. For either nucleus, the analytic shift is `delta_n=w_n*w1n^2/(w_n^2-w_c^2)`, where `w_n` is its signed Zeeman frequency, `w_c` the signed carrier frequency, and `w1n` the proton-drive amplitude scaled by that nucleus’s magnetogyric-ratio ratio. The negative 15N magnetogyric ratio gives an opposite phase direction to 13C.

## Method

The script propagates the driven state for 20 ms with the Ramsey correction enabled and compares the two accumulated phases with the analytic expression. It also checks the opposite signs, the fourfold phase change on doubling drive amplitude, and the twofold change on halving the magnetic field. These are checks in the example code, not independently rerun results.
