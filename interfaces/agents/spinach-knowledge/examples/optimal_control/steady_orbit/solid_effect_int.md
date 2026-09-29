# examples/optimal_control/steady_orbit/solid_effect_int.m

Source: [MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/steady_orbit/solid_effect_int.m) · Function: `solid_effect_int()`

The example defines a phase-controlled, stroboscopic steady-state DNP design for an electron–proton pair. It builds the spin system, declares the proton longitudinal-magnetisation state as the destination, configures powder-dependent ESR drift and control constraints, and calls the GRAPE phase objective through `fmaxnewton`. This describes the optimisation setup, not a reported optimised pulse or measured enhancement.

The model uses `E` and `1H` at 3.35316 T and 80 K. The trityl electron Zeeman principal values are [2.00319, 2.00319, 2.00258], with Euler angles [0, 10, 0] degrees; the proton shift is [0, 0, 5] ppm, with Euler angles [0, 0, 10] degrees. Their coordinates are [0, 0, 0] and [0, 0, 3.500]; the script does not label the coordinate unit. The T1/T2 relaxation configuration uses a distance- and orientation-dependent proton R1 callback from `r1n_dnp`, R1 entries {1e3, callback}, R2 entries {200e3, 50e3}, diagonal relaxation retention, and `dibari` equilibrium. The basis is `sphten-liouv` with no approximation.

The target is proton `Lz`, normalised by its overlap with the equilibrium state. The code also constructs an equilibrium initial state, while its comment says the steady-state module ignores that initial state. Powder drifts use spins {`E`, `1H`}, grid `rep_2ang_800pts_sph`, and an ESR transmitter reference at 94.0 GHz. The controls are electron `Lx/Ly`, with `Lz` as the offset operator. The listed microwave control levels span 2π·5×10⁶ to 2π·25×10⁶ rad/s (20 levels), and the five offsets are −2, −1, 0, +1, and +2 MHz.

The control vector has 720 pulse samples of 0.5 ns each (360 ns), 20 frozen 0.5 ns ringdown samples, then a frozen 167 μs sequence delay. Amplitude is one during the pulse and zero in the ringdown and delay. The configured method is `rbfgs`, with maximum 10,000 iterations, `steady=true`, and budget 500; robustness and spectrogram plots are requested. A 16-coefficient prefix of `hiper_kernel_trans.mat` is normalised to unit absolute DC gain and passed through `firf` for optimisation distortion and plotting. The system requests 240 processes and sets propagation chop tolerance 1e-14 and steady-state tolerance 1e-10; the source comments that the calculation may take days on a large parallel cluster, which is not a run result.

The source does not specify the scalar fidelity/loss formula, a converged waveform, numerical robustness or spectrogram values, or a measured DNP enhancement. The requested robustness/spectrogram displays are outputs to inspect, not results reported here; this is an example setup, not a kernel test.

This variant supplies a sinusoidal phase starting guess: the 720-point phase array is `wrapTo2Pi(4π sin(−4·linspace(−π, π, 720)))`, followed by 20 ringdown zeros and one sequence-delay zero. That is an initialisation choice; it does not assert that the final optimised phase remains sinusoidal.

The source comments identify Guinevere Mathies, Shebha-Anandhi Jegadeesan, and Ilya Kuprov as contacts: guinevere.mathies@uni-konstanz.de, shebha-anandhi.jegadeesan@uni-konstanz.de, and ilya.kuprov@weizmann.ac.il.
