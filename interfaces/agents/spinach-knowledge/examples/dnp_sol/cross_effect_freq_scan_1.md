# examples/dnp_sol/cross_effect_freq_scan_1.m

- Signature: `cross_effect_freq_scan_1()`

## Purpose

Simulates a TOTAPOL-based cross-effect dynamic nuclear polarization (DNP) system and plots the proton (S_z) expectation value as a function of microwave-frequency offset. The example is intended to reproduce Figure 2c of [the cited Journal of Magnetic Resonance paper](http://dx.doi.org/10.1016/j.jmr.2011.09.047). The source notes that differences in intensity arise from its relaxation model and from minor inconsistencies between the geometry stated in the paper and the interaction amplitudes used there.

The calculation uses electron rotating-frame dynamics with Nottingham DNP relaxation theory, as described in [the cited Applied Magnetic Resonance paper](http://dx.doi.org/10.1007/s00723-012-0367-0). The source estimates the calculation time as seconds.

## Physical / mathematical content

- The spin system contains two electrons and one proton in a 3.4 T field. The electron Zeeman scalars are 2.0023193 and 2.0021091; the proton scalar is 0.
- The specified coordinates (in the source's coordinate units) are ([0,0,0]), ([12.80,0,0]), and ([-3.12,0,3.12]) for the two electrons and proton, respectively.
- Relaxation is `nottingham`, with secular terms retained, zero equilibrium, temperature 10, and the source's listed electron and nuclear (T_1/T_2) parameters.
- The output is the proton (S_z) expectation value across 50,000 microwave offsets from -350 to 350 MHz, for a 100 kHz microwave power parameter and the specified static orientation.

## Numerical / algorithmic content

Builds the Spinach system and a full `sphten-liouv` basis with `approximation='none'`, then calls `crystal` with the `dnp_freq_scan` callback and ESR mode. It plots the real part of the returned frequency-scan signal.

## Implementation structure

The function sets the field, isotopes, Zeeman scalars, coordinates, basis, and Nottingham relaxation parameters; creates and bases the Spinach system; configures the electron microwave operators, proton detection state, frequency offsets, and orientation; evaluates the scan; and labels the plotted proton signal.
