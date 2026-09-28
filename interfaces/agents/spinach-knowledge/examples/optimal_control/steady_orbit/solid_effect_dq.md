# examples/optimal_control/steady_orbit/solid_effect_dq.m

- Signature: `solid_effect_dq()`

## Purpose

Panoramic phase optimisation for stroboscopic steady-state DNP, using timing and power settings matching the XiX experiment. The source notes a calculation time of days on a large parallel cluster.

## Physical / mathematical content

- Models an electron (`E`) and proton (`1H`) at 3.500 Å separation in a 3.35316 T W-band magnet (HiPER at St Andrews), at 80 K. The trityl electron g-tensor principal values are `[2.00319 2.00319 2.00258]`; the proton Zeeman values are `[0 0 5]` ppm.
- Uses `t1_t2` relaxation, including an orientation- and distance-dependent nuclear R1 from `r1n_dnp`; R1 rates are `{1e3, r1n_rate}` and R2 rates are `{200e3, 50e3}`. The target is proton `Lz` magnetisation normalised to its thermal-equilibrium value.
- Computes ESR drift Liouvillians over the `rep_2ang_800pts_sph` powder grid, with the transmitter set to 94.0 GHz.

## Numerical / algorithmic content

- Builds a `sphten-liouv` basis without approximation and uses 240 parallel processes, `prop_chop=1e-14`, and `stst_tol=1e-10`.
- Optimises electron `Lx` and `Ly` controls with `grape_phase` and `fmaxnewton`. The settings specify `rbfgs`, steady-state optimisation, at most 10,000 iterations, and a budget of 500. Microwave power levels are `2π × linspace(5,25,20) × 10^6` rad/s; microwave offsets are −2, −1, 0, +1, and +2 MHz.
- The sequence comprises 720 pulse samples and 20 ringdown samples at 0.5 ns each, followed by a 167 µs delay. Only pulse phases are unfrozen; amplitudes are one during the pulse and zero thereafter. The initial phase guess is a wrapped, negative 140 MHz phase ramp across the 720 pulse samples.
- Applies a 16-tap FIR filter loaded from `hiper_kernel_trans.mat`, normalised to unit absolute DC gain, for optimisation and plotting. Plotting requests robustness and spectrogram views.

## Implementation structure

The function creates the spin system and control operators, obtains powder-orientation drifts, configures the control problem with `optimcon`, and runs the optimisation. The resulting `pulse_profile` is assigned locally; the source does not explicitly save the workspace.

## Authors

- guinevere.mathies@uni-konstanz.de
- shebha-anandhi.jegadeesan@uni-konstanz.de
- ilya.kuprov@weizmann.ac.il