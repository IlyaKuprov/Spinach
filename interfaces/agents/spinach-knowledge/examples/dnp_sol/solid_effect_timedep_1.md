# examples/dnp_sol/solid_effect_timedep_1.m

- Signature: `solid_effect_timedep_1()`
- Source: [`examples/dnp_sol/solid_effect_timedep_1.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/solid_effect_timedep_1.m)

## Purpose

Demonstrates solid-effect DNP time evolution and steady state for one electron and three protons arranged in a tilted linear chain. The source places the protons 7, 10, and 14 Å from the electron at the origin and estimates the calculation time as seconds.

## Spin system and relaxation

The field setting is `sys.magnet=3.4` (the source labels it a magnetic field without giving a unit). The chain is rotated by `euler2dcm(pi/6,pi/7,pi/8)`; the spin order is `E, 1H, 1H, 1H`, with coordinates [0,0,0], [0,0,7], [0,0,10], and [0,0,14] multiplied by that rotation. Relaxation uses the Weizmann model, secular retention, IME equilibrium, and temperature 4.2. The assigned rates are `weiz_r1e=1e2`, `weiz_r1n=0.1`, `weiz_r2e=1e5`, and `weiz_r2n=1e3`. The 4-by-4 dipolar rate matrices start at zero; symmetric entries for the two adjacent proton pairs (spins 2–3 and 3–4) are set to 0.1 for both R1 and R2. Units for the rate values and temperature are not stated in the source.

The basis is `sphten-liouv` with no approximation and projections [-2, -1, 0, 1, 2]. The experiment sets `mw_pwr=2*pi*250e3`, `nuclear_frq=2*pi*144.76e6`, second-order Krylov–Bogolyubov theory (`kb_second_order`), a 0.01 s time step, and 1000 steps.

## Computation and output

After `create` and `basis`, the function calls `solid_effect` first with `calc_type='time_dependence'`, then with `calc_type='steady_state'`. The first result is plotted as real longitudinal expectation values against 1001 times from 0 to 10 s: electron separately and the three protons together, with logarithmic time axes. The second result replaces the first and its real values are printed under the label `Tr(Sz*rho) on all spins`. The function takes no arguments, returns no explicit output, and does not save the figure or numerical arrays. Running it requires the Spinach model/basis, solid-effect, and plotting routines on the MATLAB path.
