# examples/optimal_control/case_studies/Tosner_JMR_2009/bb_refocusing_pulse.m

- Signature: `bb_refocusing_pulse()`

## Purpose

Spinach implementation of the broadband refocusing example from GRAPE is used to design a 200 µs broadband x-phase π pulse: {Sx -> Sx, Sy -> -Sy, Sz -> -Sz} over an offset range of ±12.5 kHz.

## Implementation

- Implements the refocusing example from [Tosner et al., JMR (2009)](http://dx.doi.org/10.1016/j.jmr.2008.11.020). Builds a 14.1 T single-proton Liouville-space model and uses Lz as the offset operator on 101 design offsets from −12.5 to +12.5 kHz.
- Optimizes Cartesian Lx/Ly controls for 600 time steps over 200 µs, with a 2π·30 kHz power level, L-BFGS, and at most 200 iterations. The objective implements the stated x-phase pi-refocusing map.
- Tests the pulse at 201 offsets from −25 to +25 kHz with shaped_pulse_xy and plots the fidelity profile.
