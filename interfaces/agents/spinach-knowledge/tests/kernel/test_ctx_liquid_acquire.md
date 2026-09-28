# tests/kernel/test_ctx_liquid_acquire.m

- Signature: `result=test_ctx_liquid_acquire()`

## Purpose

Checks that `liquid()` forwards the offset Liouvillian correctly by comparing its FID with a direct call to `acquire()` using the same spin system and dynamics objects.

## Physical / mathematical content

The fixture is one `1H` spin at 14.1 T with zero isotropic Zeeman coupling, in the `sphten-liouv` formalism with no approximation and projection `+1`. The test focuses on the wrapper/context path rather than a new physical model.

## Numerical / algorithmic content

The four-point acquisition starts with and detects `L+`, uses zero sweep center offset parameter except for `parameters.offset=25`, a 1000 Hz sweep, and empty `decouple`, `needs`, and `rframes` fields. The context result is `liquid(spin_system,@acquire,parameters,`nmr`)`. The reference path applies `assume(...,`nmr`)`, builds `H`, `R`, and `K`, applies `frqoffset()` to `H`, then calls `acquire(spin_system,parameters,H,R,K)`.

## Check

The full context and direct FIDs must agree to absolute and relative tolerance `1e-12`.
