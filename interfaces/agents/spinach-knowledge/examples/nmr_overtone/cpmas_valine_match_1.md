# examples/nmr_overtone/cpmas_valine_match_1.m

- MATLAB implementation: [examples/nmr_overtone/cpmas_valine_match_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_overtone/cpmas_valine_match_1.m)

## Purpose

This example calculates proton-to-14N overtone cross-polarisation in N-acetylvaline under magic-angle spinning (MAS), using Spinach's Fokker-Planck formalism. It evaluates a Hartmann-Hahn profile versus proton RF power on the source's rough powder grid. The source estimates hours of calculation time.

The MATLAB source attributes the valine quadrupolar tensor data to https://doi.org/10.1039/c4cp03994g. That is literature provenance for the model input; the spectra described here are simulations, not measurements or a claim of experimental reproduction.

## Spin model and settings

The source defines two spins, 14N and 1H, with sys.magnet set to 14.10220742. The nitrogen coupling is supplied as eeqq2nqi(3.21e6,0.27,1,[0 0 0]); its Zeeman-eigenvalue input is [57.5 81.0 227.0], with Euler angles [-90 -90 -17]*(pi/180). The proton Zeeman-eigenvalue input is [0 0 0] and its Euler angles are zero. Coordinates are [0.00 0.00 0.00] and [1.00 0.00 0.00], as coded. The source does not annotate units for these tensor and coordinate entries, so none are added here.

Relaxation is damp, with rlx_keep='diagonal', zero equilibrium, and damp_rate=2000. The basis is sphten-liouv with approximation none. The magic-angle parameter is atan(sqrt(2)); the spectrum axis input is [sqrt(2/3) 0 sqrt(1/3)]. The overtone irradiation spin is 14N.

## Rotor, RF, and spectrum workflow

The simulation sets max_rank=9, rotor-rate input -19840, and grid rep_2ang_800pts_sph. The spectrum sweep is [70.0e3 105.0e3] with 256 points and 256-point zero-fill; axis_units is kHz. Initial state and proton RF operator use the magic-angle-weighted Lz and Lx components. The receiver and nitrogen RF operator are weighted the same way for 14N. The 14N overtone uses average treatment, rf_frq=86.30e3, and rf_dur=7e-5.

Fifteen proton-power settings span linspace(26e3,40e3,15). At each setting, the source assigns rf_pwr=2*pi*[55e3 rf_powers(n)]/sin(theta) and calls singlerot with overtone_cp and the qnmr algorithm.

## Output and scope

Each simulated spectrum is plotted as its real part against the configured sweep, with the displayed x limits set to 70-105 and y limits to -8.7843e-4 through 8.7843e-4. The loop labels panels with the scanned proton power in kHz. This is a one-dimensional simulated spectrum for each power setting; the source does not report a measured match or a fit statistic.
