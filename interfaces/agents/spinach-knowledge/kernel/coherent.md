# kernel/coherent.m

- Signature: `rho=coherent(spin_system,mode,alpha)`

## Purpose

Builds the normalised, Fock-space-truncated coherent state with amplitude `alpha` on the specified bosonic mode, with unit operators on all other particles of the system.

## Physical / mathematical content

- The Fock state amplitudes up to the truncation level are proportional to `alpha^n/sqrt(n!)`.
- Fock space truncation removes the tail of the Poisson distribution. The lost weight is reported, and the truncated state is renormalised.

## Numerical / algorithmic content

- Constructs the mode density matrix from the normalised amplitudes and takes its Kronecker product with identity operators on the other particles.
- Returns a full density matrix in `zeeman-hilb` formalism or its vectorisation in `zeeman-liouv` formalism.

## Parameters / inputs

- `mode` - index of a bosonic mode in `sys.isotopes`.
- `alpha` - coherent state amplitude, a complex scalar.

## Outputs

- `rho` - coherent state density matrix (`zeeman-hilb`) or its vectorisation (`zeeman-liouv`).

## Implementation structure

- Checks that basis information is present, `mode` indexes a bosonic particle, and `alpha` is a finite numeric scalar.
- Computes and normalises the truncated amplitudes, reports the lost Poisson-distribution weight, builds the composite density matrix, and converts it to the current Zeeman formalism. Other formalisms produce an error.

<https://spindynamics.org/wiki/index.php?title=coherent.m>
