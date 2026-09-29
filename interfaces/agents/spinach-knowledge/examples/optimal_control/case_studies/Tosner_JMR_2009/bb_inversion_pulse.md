# examples/optimal_control/case_studies/Tosner_JMR_2009/bb_inversion_pulse.m

- Signature: `bb_inversion_pulse()`
- Source: [examples/optimal_control/case_studies/Tosner_JMR_2009/bb_inversion_pulse.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/case_studies/Tosner_JMR_2009/bb_inversion_pulse.m)

## Purpose

This liquid-state NMR example formulates robust broadband inversion of a single proton over transmitter offsets. Its target operation is longitudinal magnetisation inversion, +Iz to -Iz. The source describes it as a Spinach reproduction of the second example associated with the Journal of Magnetic Resonance DOI [10.1016/j.jmr.2008.11.020](https://doi.org/10.1016/j.jmr.2008.11.020).

## Model and optimisation

The model is one 1H spin at 14.1 T, with zero scalar chemical shift, in the NMR rotating frame. The initial and target states are normalised +Sz and -Sz. The pulse uses Cartesian Lx and Ly controls and an Lz offset operator, with 101 equally spaced offsets from -50 to +50 kHz in the design ensemble.

The design has 600 time slices of 1 microsecond each, for a total duration of 600 microseconds. The control power level is 2*pi*10 kHz in angular-frequency units. Starting from a random 2-by-600 guess scaled by 1/10, the example configures the L-BFGS method, the NS and SNSA penalties with weights 0.01 and 10, respectively, and a maximum of 200 iterations. It calls fmaxnewton with the GRAPE XY objective.

## Verification observable

After optimisation, the script rescales the two waveform channels to physical angular-frequency controls and simulates the pulse at 201 offsets from -100 to +100 kHz. At each offset it evaluates -real(Sz' * rho_f), the inversion transfer relative to the normalised initial state, and plots that profile. The wider verification grid is distinct from the +/-50 kHz optimisation ensemble. The source defines this evaluation but supplies no numerical profile values in the page; no performance result or experimental validation is asserted here.
