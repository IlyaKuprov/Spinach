# kernel/utilities/sorensen.m

## Purpose

Computes the Sorensen bound for the maximum transfer efficiency between two states under arbitrary control operators, implementing Equation 186 of Sorensen's review ([DOI: 10.1016/0079-6565(89)80006-8](https://doi.org/10.1016/0079-6565(89)80006-8)).

## Behaviour

- Syntax: `b=sorensen(rho_init,rho_targ)`.
- Validates the inputs through an internal `grumble` subfunction, which errors with `'the inputs must be Hermitian matrices of the same size.'` unless both arguments are numeric, Hermitian, and of equal element count.
- Diagonalises both density matrices with `eig(full(...))`, extracts the eigenvalues, and sorts each set in ascending order.
- Computes the bound as `b=(sigma_init'*sigma_targ)/trace(rho_targ^2)`, i.e. the inner product of the sorted eigenvalue vectors divided by the trace of the square of the target density matrix.
- The header notes this is an exact unitary bound; the amount reachable with realistically available instrumental controls may be smaller, with a detailed analysis at [DOI: 10.1080/00268979909483117](https://doi.org/10.1080/00268979909483117).

## Inputs and outputs

**Inputs**

- `rho_init` — initial density matrix, Hilbert space.
- `rho_targ` — target density matrix, Hilbert space.

**Outputs**

- `b` — Sorensen bound.

## References

- Source: [kernel/utilities/sorensen.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/sorensen.m)
- [Spinach Wiki: sorensen.m](https://spindynamics.org/wiki/index.php?title=sorensen.m)
- O. W. Sørensen, *Polarization transfer experiments in two-spin systems*, J. Magn. Reson. — [DOI: 10.1016/0079-6565(89)80006-8](https://doi.org/10.1016/0079-6565(89)80006-8)
- [DOI: 10.1080/00268979909483117](https://doi.org/10.1080/00268979909483117)
