# examples/dnp_sol/solid_effect_timedep_2.m

- Signature: `solid_effect_timedep_2()`
- Source: [`examples/dnp_sol/solid_effect_timedep_2.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/solid_effect_timedep_2.m)

## Purpose

Extends the solid-effect DNP time-evolution example to a larger, tilted electron–proton system, using second-order Krylov–Bogolyubov average-Hamiltonian theory and a basis restricted to five-spin orders. The header estimates minutes on a Tesla A100 GPU. There is an important source mismatch: the header describes `n=2:6`, but the executable loop uses `n=2:7`; the latter creates six protons at 8, 13, 20, 29, 40, and 53 Å, not five. The actual loop therefore builds seven spins including the electron.

## Spin system and relaxation

The field setting is `sys.magnet=3.4` (the source labels it a magnetic field without stating a unit). The electron is at the origin; all six proton coordinates are [0,0,4+n^2] multiplied by `euler2dcm(pi/6,pi/7,pi/8)` for `n=2:7`. Relaxation uses the Weizmann model, secular retention, IME equilibrium, and temperature 4.2. The source assigns `weiz_r1e=1e2`, `weiz_r1n=0.1`, `weiz_r2e=1e5`, and `weiz_r2n=1e3`; symmetric entries for adjacent nuclei in the 7-by-7 R1/R2 dipolar arrays (pairs 2–3 through 6–7) are set to 0.1. The rate and temperature units are not specified.

The basis uses `sphten-liouv` with `IK-0` approximation, interaction level 5, and projections [+2, +1, 0, -1, -2]. The experiment parameters are `mw_pwr=2*pi*250e3`, `nuclear_frq=2*pi*144.76e6`, theory `kb_second_order`, 0.01 s time step, and 1000 steps. The code explicitly sets `sys.disable={'krylov'}`; the line enabling the GPU is commented out. Thus the header's A100 timing note is not, by itself, evidence that this invocation enables GPU execution.

## Computation and output

The function constructs the system and basis, calls `solid_effect` for time dependence, and plots real longitudinal expectation values on a logarithmic time axis from 0 to 10 s: the electron in one panel and the six proton channels in the other. It then requests a steady state and prints the real `Tr(Sz*rho)` values for all seven spins. It takes no arguments and returns no explicit output; the source does not save the arrays or figure. It depends on Spinach's model/basis, solid-effect, and plotting routines being available on the MATLAB path.
