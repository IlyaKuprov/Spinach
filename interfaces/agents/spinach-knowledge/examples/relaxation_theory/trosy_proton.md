# examples/relaxation_theory/trosy_proton.m

- Signature: `trosy_proton()`

## Purpose

Calculate transverse relaxation rates as a function of applied magnetic field for the C-H group at position 3 of the tyrosine aromatic ring. Calculation time: minutes.

## Physical / mathematical content

- The example reads a 3-fluorotyrosine DFT calculation from `../standard_systems/amino_acids/tyr.log`, maps carbon and hydrogen to `13C` and `1H`, and extracts the selected proton and carbon coordinates and shielding tensors for a two-spin system.
- Relaxation uses the `redfield` model with `labframe` terms retained, `zero` equilibrium, and a correlation time of `25e-9` s.
- For each field, the script evaluates relaxation matrix elements as `-v'*R*v` for normalized proton and carbon raising-operator states and for their respective two-spin combinations with `Lz` on the other spin.

## Numerical / algorithmic content

- The basis uses `sphten-liouv` formalism with `none` approximation; startup `hygiene` checks are disabled.
- A grid of 20 proton Larmor frequencies from 200 to 800 MHz is converted to magnetic fields using `2*pi*lin_freq*1e6/spin('1H')`. At each field, the script creates the spin system, constructs its basis, and computes the relaxation superoperator.
- Separate plots show the three proton and three carbon relaxation matrix elements against proton Larmor frequency. The proton plot uses a 0–160 Hz vertical range; the carbon plot uses 0–400 Hz.