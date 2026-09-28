# examples/nmr_solids/cp_contact_mas_nhh.m

- Signature: `cp_contact_mas_nhh()`

## Purpose

Simulates spinning-powder ¹H→¹⁵N cross-polarisation in the doubly rotating frame for a ¹⁵N coupled to eight protons distributed around it. The calculation uses a restricted Liouville space retaining correlations through three spins. The source estimates minutes on an NVIDIA Tesla A100 and much longer on a CPU.

## Physical / mathematical content

The nine-spin system is specified by isotropic shifts and explicit proton/nitrogen coordinates, with the proton environment arranged on an approximately 2 Å sphere around ¹⁵N; the temperature is 298 K. The experiment applies proton and nitrogen spin-lock fields along the MAS axis `[sqrt(2/3) 0 sqrt(1/3)]` at a 10 kHz rotor rate and observes the ¹⁵N signal.

## Numerical / algorithmic content

The basis uses `sphten-liouv`, approximation `IK-0`, and `inter_level=3`; interaction and proximity cutoffs are 5.0 and 4.0, respectively. The source disables `trajlevel`, enables `greedy`, selects the `rep_2ang_100pts_sph` powder grid, and runs `singlerot` with `@cp_contact_hard`. It uses `max_rank=3`, 100 steps of 10 μs, RF powers of 50 kHz on ¹H and 40 kHz on ¹⁵N, and requests `iso_eq`.

## Implementation structure

- Defines the ¹H₈–¹⁵N system, shifts, coordinates, and temperature.
- Constructs the restricted basis and Spinach system, then obtains the CP irradiation and detection operators.
- Sets MAS, powder grid, RF powers, equilibrium requirement, and time steps.
- Simulates the FID and plots its real part against elapsed time.
