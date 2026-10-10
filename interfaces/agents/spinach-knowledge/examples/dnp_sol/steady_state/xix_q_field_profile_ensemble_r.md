# examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_r.m

- Signature: `xix_q_field_profile_ensemble_r()` (no arguments)
- Source: [MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_r.m)
- Method: [steady-state XiX implementation](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/hyperpol/xixdnp_steady.m), which cites [10.1021/jacs.1c09900](https://doi.org/10.1021/jacs.1c09900).

## Purpose and protocol

Calculates a steady-state XiX DNP microwave-offset profile at one fixed electron nutation frequency, averaged over an electron–proton distance ensemble. Unlike the B1 ensemble variants, it fixes `parameters.irr_powers=18e6` Hz.

## Setup and scan

The zero-argument function creates an `E`/`1H` pair with `sys.magnet=1.2142`, spin temperature `80`, trityl electron Zeeman principal values `[2.00319 2.00319 2.00258]`, proton shift guess `[0 0 5]` ppm, and Euler inputs `(pi/180)*{[0 10 0],[0 0 10]}`. It samples four distances with `gaussleg(3.5,20,3)`; this source calls the variable a distance but does not annotate a unit. At each point the proton coordinate is set to z=`r(n)`. Relaxation uses `t1_t2`, electron `r1=1e3`, proton `r1n_dnp(sys.magnet,inter.temperature,2.00230,1e-3,52,r(n),bet)`, `r2_rates={200e3 50e3}`, diagonal retention, and `dibari` equilibrium. Units for magnet, temperature, and the rates are not specified in this source.

The script uses the unapproximated `sphten-liouv` basis, proton `Lz` detection, and `rep_2ang_800pts_sph` for powder averaging. It evaluates 201 microwave resonance offsets from −100e6 to 100e6 Hz with `powder(...,@xixdnp_steady,...,'esr')`. XiX settings are a 48e-9 s pulse duration, 36 blocks, inverted second-pulse phase `pi`, additional shift −13e6 (unit not annotated), and shot spacing `204e-6 - 2*nloops*pulse_dur`.

## Dependencies and output

Requires Spinach and `r1n_dnp`, `gaussleg`, `powder`, and `xixdnp_steady`. The three profiles are combined with distance quadrature weights and the radial `r^2` factor, then normalised by the corresponding weighted sum. The plot shows real proton `Lz` against offset in MHz and saves `xix_q_field_profile_ensemble_r.fig` in the MATLAB current directory. Only the figure is saved; the distance ensemble is finite and discretely sampled.
