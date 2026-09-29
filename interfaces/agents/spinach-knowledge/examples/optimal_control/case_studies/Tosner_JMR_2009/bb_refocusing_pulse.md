# examples/optimal_control/case_studies/Tosner_JMR_2009/bb_refocusing_pulse.m

- Signature: `bb_refocusing_pulse()`
- Source: [examples/optimal_control/case_studies/Tosner_JMR_2009/bb_refocusing_pulse.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/case_studies/Tosner_JMR_2009/bb_refocusing_pulse.m)

## Purpose

This example designs a broadband x-phase pi refocusing operation for liquid-state NMR. The specified Cartesian-state map is Sx to Sx, Sy to -Sy, and Sz to -Sz over the offset ensemble. The source identifies it as the broadband refocusing example associated with the Journal of Magnetic Resonance DOI [10.1016/j.jmr.2008.11.020](https://doi.org/10.1016/j.jmr.2008.11.020).

## Model and optimisation

The model is one 1H spin at 14.1 T, with zero scalar chemical shift, in the NMR rotating frame. It uses normalised Sx, Sy, and Sz states as its three initial states, with the corresponding target states Sx, -Sy, and -Sz. Cartesian Lx and Ly controls act in the presence of the Lz offset operator. The design ensemble contains 101 equally spaced offsets from -12.5 to +12.5 kHz.

The waveform has 600 equal slices over 200 microseconds, so each slice is 200/600 microseconds. The control power scale is 2*pi*30 kHz in angular-frequency units. From a random 2-by-600 initial guess, the example configures L-BFGS with the NS and SNS penalties weighted 0.01 and 100, respectively, and a maximum of 200 iterations, then calls fmaxnewton with the GRAPE XY objective.

## Fidelity profile

The optimised channels are rescaled by the control power level and tested at 201 offsets from -25 to +25 kHz. For each offset, the script applies the pulse to all three Cartesian initial states and computes the real trace overlap of the three specified targets with their final states, averaged by division by three. It plots this fidelity profile. The wider test range is not the design range, and the source provides no numerical fidelity values here; no experimental validation is claimed.
