# examples/esr_sol_pulsed/hard_3_pulse_deer_gd_2.m

- Signature: `hard_3_pulse_deer_gd_2()`

## Purpose

Simulates a three-pulse DEER experiment on a pair of Gd(III) centres. The probe pulse excites an outer ZFS transition and the pump pulse excites the central transition; both are treated as hard pulses. The calculation uses brute-force time propagation and powder averaging.

## Model and calculation

- Two `E8` spins are placed 30 Å apart at 3.5 T. Each has an axial ZFS tensor with principal values proportional to `[1, 1, -2]`; their orientations differ.
- The script constructs the probe and pump transition operators, sets a 100-step trace with 20 ns steps, and averages over the 1,600-point spherical grid.
- It apodises the hard-pulse, probe, and pump signals, Fourier transforms them for excitation-profile plots, and plots the DEER time trace.

The source notes that the Gd spin echo can be very sharp; the simulation omits the ZFS distributions present in experimental systems. Estimated calculation time: minutes.
