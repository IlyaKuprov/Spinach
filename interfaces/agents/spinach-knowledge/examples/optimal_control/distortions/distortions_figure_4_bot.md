# examples/optimal_control/distortions/distortions_figure_4_bot.m

- MATLAB implementation: [examples/optimal_control/distortions/distortions_figure_4_bot.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/distortions/distortions_figure_4_bot.m)

- Signature: distortions_figure_4_bot()
- Source: [examples/optimal_control/distortions/distortions_figure_4_bot.m](../../../../../../examples/optimal_control/distortions/distortions_figure_4_bot.m)
- Paper: [Rasulov and Kuprov, arXiv:2502.02198](https://arxiv.org/abs/2502.02198)

## Purpose

Reconstructs the amplifier-saturation robustness calculation for the bottom panel of Figure 4. It optimises a broadband carbon XY pulse against a small ensemble of tanh amplifier-saturation parameters, then simulates the pulse on a wider saturation-versus-RF-amplitude grid and plots state-transfer infidelity.

## Spin system and pulse design

The model contains 100 non-interacting 13C spins at evenly spaced offsets from −100 to +100 ppm in a 28.18 T field. It uses Spinach's sphten-liouv formalism with IK-2 and projection level 1. The three normalised starting operators Sx, Sy and Sz target −Sz, Sy and Sx. Lx and Ly are the two quadrature controls on the 13C channel, with a common drift Hamiltonian.

The waveform has 125 piecewise-constant 1 μs intervals. Its last five slices are frozen as dead time; those guess samples are set to 1e−3. Nominal power levels are 50–70 kHz nutation frequency per channel (represented internally as angular frequency, 2π times Hz). The starting XY profile is read from guess.mat. optimcon, fmaxnewton and grape_xy run LBFGS-GRAPE with NS and SNS penalties weighted 0.01 and 0.10, up to 50 iterations.

## Saturation model and benchmark

During optimisation, five amplifier saturation parameters span 0.9–1.1 times the mean nominal angular-frequency level. Each ensemble member supplies an amp_tanh transformation for the XY control. The example calls that helper but does not define its mathematical form here, so no further transfer-function details are implied.

The post-optimisation benchmark sweeps 41 saturation values from 0.5 to 1.5 times the same mean and 41 nominal nutation frequencies from 40 to 80 kHz. At each point it rescales the optimised waveform, applies amp_tanh, splits the quadratures, and propagates the three starting states with shaped_pulse_xy (expv-pwc). The plotted heat map is log10 of the mean three-state infidelity, calculated as one minus the real target-overlap average. White guide lines mark 50 and 70 kHz and relative saturation factors 0.9 and 1.1. Optimisation also requests XY-control, robustness and spectrogram plots.

This is a numerical amplifier model and spin-dynamics simulation; the plotted map is not a measurement of an amplifier or a claim of hardware validation. The initial guess file and amp_tanh implementation are external dependencies of the example.
