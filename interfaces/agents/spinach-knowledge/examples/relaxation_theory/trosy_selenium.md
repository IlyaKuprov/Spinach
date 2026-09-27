# examples/relaxation_theory/trosy_selenium.m

- Signature: `trosy_selenium()`

## Purpose

Calculate transverse relaxation matrix elements for `77Se` and its directly bonded `13C` in ethylselenol as a function of applied magnetic field. Calculation time: minutes.

## Physical / mathematical content

- The two-spin model uses coordinates and Zeeman shielding tensors extracted from `../standard_systems/ethylselenol.out`.
- Relaxation uses the `redfield` model with `rlx_keep='labframe'`, `equilibrium='zero'`, and correlation time `tau_c=25e-9` s.
- For each field, the code evaluates `-v'*R*v` for normalized single-spin `L+` states and for normalized combinations `Se+ ± 2 Se+ Cz` and `C+ ± 2 C+ Sez`, where `R` is the relaxation superoperator.

## Numerical / algorithmic content

- The basis uses `sphten-liouv` formalism with `approximation='none'`; startup hygiene checks are disabled.
- A grid of 20 proton Larmor frequencies from 200 to 800 MHz is converted to magnetic fields using `B0=2*pi*lin_freq*1e6/spin('1H')`. At each field, the spin system and basis are created, `R` is calculated, and six relaxation matrix elements are evaluated.
- Separate plots show the three selenium and three carbon matrix elements against proton Larmor frequency in MHz; the vertical axis is relaxation matrix element in Hz.