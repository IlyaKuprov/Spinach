# examples/nmr_overtone/cpmas_valine_match_2.m

- MATLAB implementation: [examples/nmr_overtone/cpmas_valine_match_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_overtone/cpmas_valine_match_2.m)

## Purpose

This example calculates a proton-to-14N overtone cross-polarisation Hartmann-Hahn profile for N-acetylvaline under MAS with the Fokker-Planck formalism. Unlike the one-dimensional power scan in the companion example, this script scans both proton RF power and sample spinning rate and plots a two-dimensional simulation-derived intensity map. The source estimates hours of calculation time.

The MATLAB source attributes the valine quadrupolar tensor data to https://doi.org/10.1039/c4cp03994g. This is the model's cited literature provenance; the map is calculated output, not experimental data or a claim of agreement with an experiment.

## Spin model and settings

The source defines 14N and 1H at sys.magnet=14.10220742. The nitrogen coupling is eeqq2nqi(3.21e6,0.27,1,[0 0 0]); its Zeeman-eigenvalue input is [57.5 81.0 227.0], with Euler angles [-90 -90 -17]*(pi/180). The proton Zeeman-eigenvalue input and Euler angles are zero. Coordinates are [0.00 0.00 0.00] and [1.00 0.00 0.00], as coded. No units are added to tensor or coordinate inputs because the source does not annotate them inline.

Relaxation is damp, retaining the diagonal terms with zero equilibrium and damp_rate=10000. The basis is sphten-liouv with no approximation. The magic angle is atan(sqrt(2)); the spectrum axis input is [sqrt(2/3) 0 sqrt(1/3)]. The overtone irradiation spin is 14N.

## Rotor and RF scan

The calculation uses max_rank=5, grid rep_2ang_200pts_oct, and axis_units=kHz. It builds 50 proton-power values from linspace(10e3,200e3,50) and 50 spinning-rate values from linspace(20e3,90e3,50). At each pair, the RF power is assigned as 2*pi*[55e3 rf_powers(n)]/sin(theta) and the local rotor-rate input is assigned -spin_rates(k). The nitrogen RF frequency input is 46.30e3-2*localpar.rate, and the local sweep is that frequency plus or minus 15e3. The spectrum uses 256 points and 256-point zero-fill; RF duration is 1e-4. The source does not state additional units for the rotor-rate or RF-frequency inputs in these assignments.

For each grid point, singlerot runs overtone_cp with qnmr; the stored value is sum(real(spectrum)). The calculation uses a parfor loop over the power-rate pairs.

## Output and scope

The plotted image has proton RF power in Hz on the x axis and sample spinning rate in Hz on the y axis, with the code setting the y direction to normal. Pixel values are the summed real parts of the simulated spectra, not fitted parameters or experimental intensities. In particular, this output is a coarse scan on the specified powder grid rather than a plotted collection of full spectra.
