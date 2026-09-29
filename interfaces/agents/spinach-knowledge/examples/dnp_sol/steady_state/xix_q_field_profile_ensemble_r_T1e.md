# examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_r_T1e.m

- Signature: `xix_q_field_profile_ensemble_r_T1e()` (no arguments)
- Source: [MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_r_T1e.m)
- Method: [steady-state XiX implementation](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/hyperpol/xixdnp_steady.m), which cites [10.1021/jacs.1c09900](https://doi.org/10.1021/jacs.1c09900).

## Purpose and protocol

Overlays five steady-state XiX DNP microwave-offset profiles for different electron longitudinal relaxation times, with each curve averaged over an electron–proton distance ensemble. It is the T1e-sweep variant of the fixed-B1 distance example: B1 remains 18e6 Hz and each distance/T1e pair runs the same XiX sequence. The source estimates minutes of calculation time.

## Setup and scan

The outer function scans `T1e=[10e-3,3.0e-3,1.0e-3,0.3e-3,0.1e-3]` seconds and calls its local profile routine for each value. Each run builds the `E`/`1H` system with `sys.magnet=1.2142`, spin temperature `80`, trityl electron Zeeman values `[2.00319 2.00319 2.00258]`, proton shift guess `[0 0 5]` ppm, and Euler inputs `(pi/180)*{[0 10 0],[0 0 10]}`. It samples four distances from `gaussleg(3.5,20,3)`; the source does not annotate their unit. Relaxation uses electron `r1=1/T1e`, proton `r1n_dnp(sys.magnet,inter.temperature,2.00230,T1e,52,r(n),bet)`, `r2_rates={200e3 50e3}`, diagonal retention, and `dibari` equilibrium. Magnet, temperature and rate units are not stated in this source.

Each run uses fixed electron nutation frequency 18e6 Hz, 201 microwave offsets spanning −100e6 to 100e6 Hz, the unapproximated `sphten-liouv` basis, proton `Lz` detection, and `rep_2ang_800pts_sph` powder averaging. The protocol uses 48e-9 s pulses, 36 XiX blocks, second-pulse phase `pi`, additional shift −13e6 (unit not annotated), and shot spacing `204e-6 - 2*nloops*pulse_dur`. Distance profiles are combined using quadrature weights and the radial `r^2` factor.

## Dependencies and output

Requires Spinach and `gaussleg`, `r1n_dnp`, `powder`, and `xixdnp_steady`; the five-curve wrapper calls its own local profile helper. It plots real proton `Lz` versus offset in MHz, fixes y limits to [−3e−3,3e−3], adds a T1e legend, and saves `xix_q_field_profile_ensemble_r_T1e.fig` in the MATLAB current directory. It saves no tabulated profiles; both distance and T1e scans are finite.
