# examples/fundamentals/derivative_tests/difdiff_bs_rect.m

- Signature: `difdiff_bs_rect()`

## Purpose

Checks directional derivatives of Cartesian GRAPE with Bloch–Siegert corrections by comparing the analytical gradient with central finite differences at the first, last, and midpoint waveform samples.

## Physical / mathematical content

The test system has `sys.magnet=10.2` and isotopes `{'1H','1H','13C','13C'}`. The scalar Zeeman values are `{1.5,2.0,30.0,40.0}`; the nonzero scalar couplings are 1–2: `7.0`, 1–3: `150`, 2–4: `150`, and 3–4: `50`. It starts from the normalized singlet of spins 1 and 2 and targets normalized `Lz` on spin 4. Cartesian `Lx/Ly` controls use channels `[1,1,2,2]`, with transmitter offsets `{1050,5285}`; the rectangular integrator uses 100 pulse intervals of `1.5e-4` each, power levels `2*pi*500`, and `max_iter=1000`. Bloch–Siegert corrections are enabled for `{'1H','13C'}.

## Numerical / algorithmic content

A random control waveform is evaluated with `grape_xy`. Its analytical gradient is compared at the left edge, right edge, and midpoint with central differences using `h=1e-5`. Each relative error must be below `5e-6`; the script reports individual test outcomes and errors if any check fails.

## Implementation structure

The script creates and optimizes the spin system and control configuration, forms the random waveform, and computes the analytical gradient once. It perturbs one waveform sample at a time in both directions, reevaluates fidelity with `grape_xy`, and checks the three selected positions.
