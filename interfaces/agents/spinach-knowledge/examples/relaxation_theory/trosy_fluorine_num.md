# examples/relaxation_theory/trosy_fluorine_num.m

- Signature: `trosy_fluorine_num()`

## Purpose

Calculate transverse relaxation rates as a function of applied magnetic field for the fluorine atom and its directly bonded carbon in a 3-fluorotyrosine-labelled protein. Calculation time: minutes.

## Physical / mathematical content

- The two-spin system contains `19F` and `13C`. Their coordinates and Zeeman shielding matrices are extracted from a 3-fluorotyrosine DFT calculation.
- Relaxation uses the `redfield` model with `labframe` terms, `zero` equilibrium, and a correlation time of `25e-9` s.
- At each field, rates are evaluated as `-v'*R*v` for normalized single-spin transverse states and for the corresponding states combined with longitudinal order on the other spin.

## Numerical / algorithmic content

- Read `../standard_systems/3_fluoro_tyr.log` using `gparse` and `g2spinach`, with isotope substitutions `C` to `13C` and `F` to `19F` and reference values `[186.38 192.97]`.
- Use the `sphten-liouv` formalism with `none` basis approximation; disable `hygiene` startup checks.
- Sweep 20 proton Larmor frequencies from 200 to 800 MHz, converting each frequency to a magnetic field with `B0=2*pi*lin_freq*1e6/spin('1H')`. Create the spin system and basis and calculate the relaxation superoperator at each field.

## Implementation structure

- Extract the fluorine and carbon shielding matrices and coordinates from DFT entries 8 and 7, respectively.
- Construct and normalize the `L+` state for each spin and the combinations `F+ - 2 F+ Cz`, `F+ + 2 F+ Cz`, `C+ - 2 C+ Fz`, and `C+ + 2 C+ Fz`.
- Calculate six relaxation matrix elements at each field: three for fluorine states and three for carbon states.
- Produce separate fluorine and carbon plots against proton Larmor frequency (MHz), with relaxation matrix elements in Hz.