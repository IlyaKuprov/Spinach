# examples/optimal_control/distortions/distortions_figure_3.m

- Signature: `distortions_figure_3()`

## Purpose

Formulate a broadband 13C optimal-control pulse over a chemical-shift range and compare simulated spectra for an undistorted GRAPE pulse, its RLC-distorted form, an RLC-aware GRAPE design, and a conventional hard pulse. The source specifies calculations and plots; this note makes no claim about their numerical outcomes.

## Physical / mathematical content

The spin system contains 100 noninteracting 13C spins at equally spaced offsets from -100 to +100 ppm in a 28.18 T field. Its NMR drift Hamiltonian is used with two quadrature controls, Lx and Ly, on the 13C channel. The three initial operators are Sx, Sy, and Sz; their target operators are -Sz, Sy, and Sx, respectively. This states the component mapping used as the common control objective across the offsets. The basis is `sphten-liouv` with `IK-2`; the source comment specifies a complete basis for each spin while omitting multi-spin orders.

## Numerical / algorithmic content

The GRAPE control is piecewise constant over 45 slices of 1 microsecond each. Both quadratures share the 13C channel; ten power levels span `2*pi*linspace(50e3,70e3,10)` rad/s. The `fmaxnewton` call uses `@grape_xy`, `lbfgs`, up to 200 iterations, `NS` and `SNS` penalties with weights 0.01 and 10, and freezes slices 41-45. The initial guess is 0.25 in both channels on slices 1-40 and zero on the final five.

For the RLC model, the code sets Q to 1000 and forms two single-pole coefficients as `exp(-abs(omega)*dt/(2*Q))`, where `omega=-sys.magnet*spin('13C')` and `dt` is one slice. Those two `spf` filters are included as distortion functions in a second GRAPE optimisation; the optimised profile is also filtered for the comparison spectrum. The acquisition calculations use a 13C FID, a 70000 Hz sweep, 2048 points, zero-filling to 16384 points, Gaussian apodisation with parameter 10, and a Fourier transform. The comparison hard pulse is specified at `2*pi*60e3` rad/s for 4.2 microseconds, with phase `pi/2` and offset frequency zero.

## Syntax

Call `distortions_figure_3()` with no arguments. Source: [examples/optimal_control/distortions/distortions_figure_3.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/distortions/distortions_figure_3.m).

## Parameters / inputs

The field, offset ensemble, basis settings, control discretisation, power levels, penalties, filter model, and acquisition settings are defined in the function. There are no function inputs.

## Outputs

The function constructs and plots four simulated 13C spectra, labelled for the hard pulse, RLC-distorted GRAPE pulse, undistorted GRAPE pulse, and RLC-distorted RAW-GRAPE pulse. The spectra are based on simulated FIDs and plotted with intensity in arbitrary units; the source contains no experimental measurement or numerical outcome table.

## Header notes

The source labels this as Figure 3 from a paper by Rasulov and Kuprov and links [arXiv:2502.02198](https://arxiv.org/abs/2502.02198).
