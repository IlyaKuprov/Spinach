# examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_b1_r.m

- Signature: `xix_q_field_profile_ensemble_b1_r()` (no arguments)
- Source: [MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_b1_r.m)
- Method: [steady-state XiX implementation](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/hyperpol/xixdnp_steady.m), which cites [10.1021/jacs.1c09900](https://doi.org/10.1021/jacs.1c09900).

## Purpose and protocol

Calculates a steady-state XiX DNP microwave-offset profile with two explicit ensemble integrations: electron–proton distance and electron nutation frequency (B1). It is the distance-plus-B1 variant of the fixed-distance B1 example. The source estimates minutes of calculation time.

## Setup and scan

The script uses an `E`/`1H` pair, `sys.magnet=1.2142`, spin temperature `80`, trityl electron Zeeman principal values `[2.00319 2.00319 2.00258]`, proton shift guess `[0 0 5]` ppm, and Euler inputs `(pi/180)*{[0 10 0],[0 0 10]}`. At each sampled distance, the proton coordinate is set to z=`r(n)`; here `gaussleg(3.5,20,3)` is explicitly annotated as Å. The proton `r1` handle calls `r1n_dnp(sys.magnet,inter.temperature,2.00230,1e-3,52,r(n),bet)`; electron `r1=1e3`, `r2_rates={200e3 50e3}`, relaxation retention is diagonal, and equilibrium is `dibari`. Magnet, temperature and rates are not assigned units in this source.

For each distance and each of six Gauss–Legendre B1 points spanning 10e6–20e6 Hz, it evaluates 201 microwave offsets from −100e6 to 100e6 Hz on `rep_2ang_800pts_sph`, using the unapproximated `sphten-liouv` basis, proton `Lz` detection, and `powder(...,@xixdnp_steady,...,'esr')`. The XiX protocol uses 48e-9 s pulses, 36 blocks, second-pulse phase `pi`, additional shift −13e6 (unit not annotated), and shot spacing `204e-6 - 2*nloops*pulse_dur`.

## Dependencies and output

Requires Spinach and `gaussleg`, `r1n_dnp`, `powder`, and `xixdnp_steady`. First the B1 results are quadrature-weighted; the distance average then applies the radial `r^2` factor with its quadrature weights and normalisation. It plots the real proton `Lz` signal against microwave offset in MHz and saves `xix_q_field_profile_ensemble_b1_r.fig` in the MATLAB current directory. No numerical profile is saved, and both ensembles use finite quadrature grids.
