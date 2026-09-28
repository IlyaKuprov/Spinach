# tests/kernel/test_zeeman_hamiltonian.m

- Signature: `result=test_zeeman_hamiltonian()`

## Purpose

Checks the sign and units of a one-proton Zeeman Hamiltonian in the Spinach NMR rotating-frame convention.

## Physical / mathematical content

- For a positive 1 ppm chemical shift, the rotating-frame Hamiltonian contribution is `-2*pi*nu*Lz`.

## Numerical / algorithmic content

- Constructs a one-proton system at `14.1` T, computes `nu` with `ppm2hz(1,sys.magnet,'1H')`, and compares `hamiltonian(assume(spin_system,'nmr'))` with `-2*pi*nu*operator(spin_system,'Lz',1)` using absolute and relative tolerances `1e-6` and `1e-12`.

## Outputs

- `result` - regression test result with explanatory messages.

## Implementation structure

- Defines the proton isotope, scalar shift, `zeeman-hilb` formalism, and `none` approximation; builds the test spin system, forms the observed and reference Hamiltonians, and records their comparison with `test_close`.
