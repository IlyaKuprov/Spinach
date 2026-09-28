# examples/dnp_sol/solid_effect_timedep_2.m

- Signature: `solid_effect_timedep_2()`

## Purpose

Simulates solid-effect DNP time dependence and steady state for one electron and a tilted chain of protons. The header estimates minutes with a Tesla A100 GPU. It describes distances as 4+n^2 Å for n=2:6, but the executable loop uses n=2:7, producing six protons at 8, 13, 20, 29, 40, and 53 Å; the source initializes seven total spins.

## Model and method

The field is 3.4 T. The chain is rotated using Euler angles [pi/6, pi/7, pi/8]. Weizmann relaxation uses IME equilibrium, secular retention, temperature 4.2, and specified electron/nuclear rates; distance-dependent R1d and R2d entries of 0.1 are assigned symmetrically between adjacent proton pairs. The basis is `sphten-liouv` with `IK-0`, five-spin-order inter-level restriction, and projections [+2, +1, 0, -1, -2].

The microwave power is 250 kHz and the nuclear frequency is 144.76 MHz. The calculation selects second-order Krylov–Bogolyubov theory, 0.01 s steps, and 1000 steps, then evaluates time dependence and steady state. The code explicitly disables Krylov and leaves the GPU-enable line commented out; the header's A100 runtime estimate should not be read as evidence that this invocation enables GPU execution.

## Outputs

The time-domain plots show the electron and six proton longitudinal expectation values on logarithmic time axes. The real steady-state `Tr(Sz*rho)` values for all spins are printed.
