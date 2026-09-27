# examples/relaxation_theory/trosy_nh.m

- Signature: `trosy_nh()`

## Purpose

Calculate transverse relaxation matrix elements as a function of magnetic field for a typical protein amide `1H`–`15N` group. The rotational correlation time is 25 ns. The nitrogen CSA parameters are from [the cited study](http://dx.doi.org/10.1021/ja0016194); the nitrogen–proton bond length is from DFT. The source estimates a calculation time of minutes.

## Physical / mathematical content

- The two-spin system uses proton and nitrogen shielding-tensor eigenvalues `[6, 0, -6]` and `[-108, 62, 46]`, respectively. Their Euler angles are `[0, 0, 0]` and `[0, 0, -19]` degrees; the spin coordinates are `[1.04, 0, 0]` and `[0, 0, 0]`.
- Relaxation is calculated with the `redfield` model, `labframe` relaxation terms, `zero` equilibrium, and `tau_c={25e-9}`.
- Six rates are evaluated as negative diagonal matrix elements of the relaxation superoperator: single-spin `H+` and `N+` coherences, plus the normalized `H+ - 2 H+ Nz`, `H+ + 2 H+ Nz`, `N+ - 2 N+ Hz`, and `N+ + 2 N+ Hz` states.

## Numerical / algorithmic content

- The calculation uses the `sphten-liouv` formalism with no basis approximation. Startup hygiene checks are disabled.
- A 30-point grid spans proton Larmor frequencies from 200 to 1500 MHz. Each frequency is converted to a magnetic field using `spin('1H')`; the spin system and basis are then built and the relaxation superoperator is calculated at that field.
- Two figures plot the three proton and three nitrogen relaxation matrix elements against proton Larmor frequency (MHz), with relaxation matrix elements labelled in Hz.