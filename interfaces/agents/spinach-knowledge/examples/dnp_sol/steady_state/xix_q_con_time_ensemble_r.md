# examples/dnp_sol/steady_state/xix_q_con_time_ensemble_r.m

- MATLAB implementation: [examples/dnp_sol/steady_state/xix_q_con_time_ensemble_r.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_q_con_time_ensemble_r.m)

- Signature: `xix_q_con_time_ensemble_r()`

## Purpose

Calculates the steady-state proton signal versus XiX contact time averaged over electron–proton distance. This is the distance-only variant: the source fixes the electron nutation frequency at 18 MHz and does not scan an electron B1 ensemble or electron T1 values.

## Model and settings

The function builds a Q-band electron–`1H` pair at `sys.magnet=1.2142` (1.2142 T), with electron Zeeman principal values `[2.00319 2.00319 2.00258]` and proton shift `[0 0 5]` (described as a ppm guess). Both Euler-angle triplets are `[0 10 0]` degrees, converted to radians. Spin temperature is 80 K; the distance quadrature is `gaussleg(3.5,20,3)` Å, and each node sets the pair coordinates to `[0 0 0]` and `[0 0 r]`.

The basis is `sphten-liouv` with no approximation; propagator chop tolerance is `1e-12`, and hygiene is disabled. Relaxation uses `t1_t2`, diagonal terms, and `dibari` equilibrium. Electron R1/R2 are set to 1000/200000, proton R2 to 50000, and proton R1 is an orientation- and distance-dependent handle calling `r1n_dnp(sys.magnet,inter.temperature,2.00230,1e-3,52,r(n),bet)`. The source does not annotate units for these rate values.

## Experiment and scan

The proton detector is `coil_state(spin_system,'Lz','1H','exact')`. XiX settings are `parameters.spins={'E','1H'}`, fixed electron nutation frequency 18 MHz, grid `rep_2ang_800pts_sph`, pulse duration 48 ns, and second-pulse phase π. The source sets `addshift=-13e6`, `el_offs=61e6`, and scans `nloops=1:64`, with shot spacing set to 153 μs minus the total pulse duration. Total contact time is `2*nloops*48e-9` s (96 ns to 6.144 μs).

## Calculation and output

For each distance node, each contact-time point is evaluated by `powder(spin_system,@xixdnp_steady,localpar,'esr')`. The code then performs a distance average with quadrature weights multiplied by the radial (r^2) Jacobian and normalises by the weighted (r^2) sum. It plots the real proton (L_z) expectation value against total contact time in μs and saves `xix_q_con_time_ensemble_r.fig`; the function has no explicit MATLAB output. The source comment estimates calculation time as minutes, not a measured runtime here.
