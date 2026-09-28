# examples/optimal_control/static_powder_control.m

- Signature: `static_powder_control()`

## Purpose

Optimise a pulse that prepares deuterium magnetisation in alanine's -CD3 group for rephasing 100 microseconds after the pulse. The optimisation covers a 100-orientation powder, RF power levels from 46 to 54 kHz, and transmitter offsets from -1 to 1 kHz. It uses Goodwin's GRAPE Hessian algorithm because the propagator dimensions are small, yielding a spin echo.

## Physical / mathematical content

- A 600 MHz magnet and an alanine -CD3 deuterium nuclear quadrupole interaction define the spin system.
- The initial deuterium `Lz` state and target `Lx` state are normalised. Deuterium `Lx` and `Ly` operators control the pulse; `Lz` supplies the offset operator.
- The optimisation includes powder orientations, five RF power levels, and five transmitter offsets.

## Numerical / algorithmic content

- The pulse has 100 slices of 2 microseconds each. `fmaxnewton` optimises a random two-channel waveform guess using `grape_xy`, the `goodwin` method, a 100-iteration limit, a 100-microsecond dead time, and `NS` and `SNS` penalties.
- A parallel loop tests the optimised pulse across the powder drifts and compares its echo with an ideal free induction decay. The signals are exponentially apodised and Fourier transformed for comparison.

## Implementation structure

- Create the deuterium spin system and powder drift Liouvillians, then configure and run the pulse optimisation.
- Apply the resulting pulse in a test calculation, plot the time-domain echo, and compare the Fourier transforms of the optimised half-echo and ideal free induction decay.
- Calculation time: minutes.
