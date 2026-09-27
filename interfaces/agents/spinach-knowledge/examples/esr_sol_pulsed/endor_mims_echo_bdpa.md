# examples/esr_sol_pulsed/endor_mims_echo_bdpa.m

- Signature: `endor_mims_echo_bdpa()`

## Purpose

Stimulated-echo stage of Mims ENDOR on BDPA, simulated without applying the nuclear pulse. The example plots the echo response over the interval from 0 to `2τ`, with `τ = 200 ns`.

## Physical / mathematical content

- The source defines the electron–proton coupling matrices and uses a `sphten-liouv` basis without approximation.

## Numerical / algorithmic content

- The powder calculation uses a 400-point spherical grid and the `endor_mims_echo` sequence. The plotted quantity is the real echo intensity as a function of time.
