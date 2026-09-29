# examples/optimal_control/features_multitarget.m

- Signature: `features_multitarget()`
- Source: [examples/optimal_control/features_multitarget.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/features_multitarget.m)

## Purpose and spin model

This singlet-state NMR simulation optimises one pulse for two transfers: carbon-triplet/proton-triplet to carbon-singlet/proton-singlet (TT to SS), and carbon-triplet/proton-singlet to carbon-singlet/proton-triplet (TS to ST). The four sites are 1H, 13C, 13C, and 1H at 14.1 T, with zero chemical shifts. The fixed couplings are 15 Hz for pairs 1-2 and 3-4, 3 Hz for 1-3 and 2-4, 150 Hz for 2-3, and 8 Hz for 1-4. The pulse is designed in a model built by the script; no imported measurements are used.

## Robust control calculation

The TT/TS source and SS/ST target operators are built from singlet/triplet terms on the proton pair (sites 1 and 4) and carbon pair (sites 2 and 3). The 1H-13C coupling for sites 1-2 is varied over 11 values from 13 to 16 Hz, producing an ensemble of drift Hamiltonians. Four x/y controls address the proton and carbon channels. The pulse has 275 intervals of 500 microseconds (137.5 ms) and uses a single power of `2*pi*500` rad/s, the SNS penalty (weight 100), `lbfgs`, and a 150-iteration limit. A random 4-by-275 guess is optimised through `fmaxnewton` with `@grape_xy`; robustness and spectrogram plots are among the enabled diagnostics. The source comments estimate hours of calculation time; this is an estimate, not a measured runtime here.

## Output and limits

The optimisation is configured with all 11 drift Hamiltonians, but the final `shaped_pulse_xy` check propagates only with the sixth drift, corresponding to a 1-2 coupling of 14.5 Hz. It reports two real final overlaps, one for each target. Those two values are single-ensemble-member simulation outputs, not an aggregate robustness score. No hardware validation is performed.
