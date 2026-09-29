# examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_b1.m

- Signature: `xix_q_field_profile_ensemble_b1()` (no arguments)
- Source: [MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_b1.m)
- Method: [steady-state XiX implementation](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/hyperpol/xixdnp_steady.m), which cites [10.1021/jacs.1c09900](https://doi.org/10.1021/jacs.1c09900).

## Purpose and protocol

Builds a steady-state XiX DNP microwave-offset profile and averages it over electron nutation-frequency (B1) quadrature points. This is the fixed-distance, B1-ensemble variant; the sibling files add distance averaging or vary the electron relaxation time. The example source estimates minutes of calculation time.

## Setup and scan

The script creates an `E`/`1H` pair, with `sys.magnet=1.2142`, spin temperature `80`, trityl electron Zeeman principal values `[2.00319 2.00319 2.00258]`, and the proton shift guess `[0 0 5]` ppm. Euler-angle inputs are `(pi/180)*{[0 10 0],[0 0 10]}`; the two coordinate rows place the proton at z=3.5 relative to the electron (no coordinate unit is stated in this source). Relaxation is `t1_t2`, with electron `r1=1e3`, proton `r1n_dnp(sys.magnet,inter.temperature,2.00230,1e-3,52,r_en,bet)`, `r2_rates={200e3 50e3}`, diagonal relaxation retention, and `dibari` equilibrium. The source does not annotate units for the magnet, temperature, or rate values.

It uses an unapproximated `sphten-liouv` basis, proton `Lz` detection, propagator chop tolerance `1e-12`, and powder grid `rep_2ang_800pts_sph`. The 201 microwave offsets span −100e6 to 100e6 Hz. Six Gauss–Legendre points sample B1 from 10e6 to 20e6 Hz; each point runs `powder(spin_system,@xixdnp_steady,parameters,'esr')`. XiX settings are 48e-9 s pulses, 36 blocks, second-pulse phase `pi`, additional shift −13e6 (the source does not annotate its unit), and shot spacing `204e-6 - 2*nloops*pulse_dur`.

## Dependencies and output

Requires Spinach setup/functions plus `gaussleg`, `r1n_dnp`, `powder`, and `xixdnp_steady`. B1 profiles are combined with the quadrature weights; the script plots `real(dnp)` (the real proton signal) against offset in MHz and saves `xix_q_field_profile_ensemble_b1.fig` in the MATLAB current directory. It saves no tabulated profile; the scan is a finite quadrature on the stated grid and offsets.
