# kernel/includes/redfield_integral_serial.m

- Signature: `(script file)`

## Purpose

This serial include evaluates Bloch-Wangsness-Redfield and Nakajima-Zwanzig relaxation integrals from the `relaxation.m` theory blocks. It is used when `relaxation.m` is not at the top of the parallelisation call stack; otherwise the asynchronous include is used. The implementation follows the notation of [the cited paper](http://dx.doi.org/10.1016/j.jmr.2010.12.004) and uses the faster auxiliary-matrix integral method described in [the second cited paper](http://dx.doi.org/10.1063/1.4928978).

## Theory parameters

- `rlx_onshell`: true selects the back-rotated kernel, which reduces to Redfield theory at zero shift; false selects the Nakajima-Zwanzig resolvent kernel.
- `rlx_shift`: the Laplace evaluation point, in Hz. Redfield theory is the on-shell form at zero shift.

## Algorithm

The code sums over spherical-rank projection pairs and correlation-function exponentials, applies chemical-species state masks, evaluates the integrals with `expmint`, and accumulates the result in `R`.

## Source documentation

https://spindynamics.org/wiki/index.php?title=redfield_integral_serial.m
