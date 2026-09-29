# examples/nmr_overtone/cpmas_valine_simple.m

- MATLAB implementation: [examples/nmr_overtone/cpmas_valine_simple.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_overtone/cpmas_valine_simple.m)

## Purpose

This script calculates a single proton-to-14N overtone cross-polarisation spectrum for N-acetylvaline under MAS using the Fokker-Planck formalism. The source estimates hours of calculation time. It attributes the valine quadrupolar tensor data to https://doi.org/10.1039/c4cp03994g; the calculated spectrum is simulation output, not a reported measurement or evidence of experimental reproduction.

## Spin model

The two-spin system is 14N and 1H at sys.magnet=14.10220742. The nitrogen coupling is supplied as eeqq2nqi(3.21e6,0.27,1,[0 0 0]); its Zeeman-eigenvalue input is [57.5 81.0 227.0] and its Euler angles are [-90 -90 -17]*(pi/180). The proton Zeeman-eigenvalue input and Euler angles are zero. The two coordinate entries are [0.00 0.00 0.00] and [1.00 0.00 0.00]. These arrays are transcribed as source inputs; the source does not give inline units for them.

The relaxation model is damp, with diagonal retention, zero equilibrium, and damp_rate=2000. The basis is sphten-liouv with approximation none. Algorithmic options disable krylov and trajlevel.

## MAS, RF, and simulated spectrum

The magic-angle parameter is atan(sqrt(2)); the spectrum axis input is [sqrt(2/3) 0 sqrt(1/3)]. The overtone irradiation spin is 14N. The spectrum uses max_rank=9, rotor-rate input -19840, grid rep_2ang_6400pts_sph, sweep [70.0e3 105.0e3], 256 points, and 256-point zero-fill. The source sets axis_units=kHz. Initial state, receiver, and RF operators use the stated magic-angle combinations of Lz and Lx for the relevant spin.

This is one parameter setting, not a power or spinning-rate scan: the proton RF-power input is 2*pi*[55.0e3 35.1e3]/sin(theta), the 14N overtone RF-frequency input is 86.30e3, and rf_dur=1e-4. The source calculates the spectrum with singlerot and overtone_cp, then plots its real part with plot_1d.
