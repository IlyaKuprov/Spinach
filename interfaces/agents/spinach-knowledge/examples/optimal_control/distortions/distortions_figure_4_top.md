# examples/optimal_control/distortions/distortions_figure_4_top.m

- MATLAB implementation: [examples/optimal_control/distortions/distortions_figure_4_top.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/distortions/distortions_figure_4_top.m)

- Signature: distortions_figure_4_top()
- Source: [examples/optimal_control/distortions/distortions_figure_4_top.m](../../../../../../examples/optimal_control/distortions/distortions_figure_4_top.m)
- Paper: [Rasulov and Kuprov, arXiv:2502.02198](https://arxiv.org/abs/2502.02198)

## Purpose

Reconstructs the RLC-distortion robustness calculation for the top panel of Figure 4. The script first optimises a broadband carbon XY pulse under an ensemble of two-stage RLC-like slice filters, then independently simulates the pulse over a wider Q-factor and RF-amplitude grid and plots the resulting state-transfer infidelity.

## Model and pulse design

The spin system contains 100 uncoupled 13C spins with evenly spaced offsets from −100 to +100 ppm at a 28.18 T field. Spinach uses the sphten-liouv formalism with the IK-2 approximation and projection level 1. Normalised Cartesian operators Sx, Sy and Sz are prepared as the three initial states; the desired outputs are −Sz, Sy and Sx. The same drift Hamiltonian is used for each state, while Lx and Ly drive the two quadratures of one 13C channel.

The control has 125 piecewise-constant intervals of 1 μs (125 μs total); the final five are frozen as dead time. Nominal channel amplitude levels span 50–70 kHz in nutation-frequency units, converted to angular frequency by multiplying by 2π. The starting waveform is loaded from guess.mat (xy_profile), with the final five samples set to zero. optimcon configures an LBFGS-GRAPE calculation via fmaxnewton and grape_xy; the run permits 25 iterations and applies the NS and SNS penalties with weights 0.01 and 0.10.

## RLC distortion and reported observable

For the optimisation ensemble, five Q values span 560–640. The source forms the per-slice attenuation parameter as exp(-abs(omega)*dt/(2*Q)), with omega=-sys.magnet*spin('13C') and dt=1 μs, and applies spf twice in succession to each quadrature. Thus the modeled response is a paired slice-filter operation parameterised by Q; it is not a measured probe response or a full hardware calibration.

After optimisation, the benchmark varies Q over 200–1000 and nominal RF nutation frequency over 40–80 kHz (41 points on each axis). For each grid point, the pulse is rescaled, passed through the two Q-dependent spf stages, and propagated with shaped_pulse_xy using expv-pwc. The plotted quantity is log10 of one minus the mean of the three real target-state overlaps. The heat map has RF nutation frequency in kHz and Q factor on its axes, with reference lines at 50/70 kHz and Q=550/650; the script also requests control, robustness and spectrogram plots during optimisation.

These are simulated responses of the stated spin model and filter, not experimental observations. The waveform initial guess and Spinach helper implementations are dependencies; their contents are not defined by this example.
