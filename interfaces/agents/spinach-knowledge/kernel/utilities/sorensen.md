# kernel/utilities/sorensen.m

- Signature: `b=sorensen(rho_init,rho_targ)`

## Purpose

Sorensen bound for the maximum transfer efficiency between two states under arbitrary control operators. Equation 186 from https://doi.org/10.1016/0079-6565(89)80006-8. Syntax: b=sorensen(rho_init,rho_targ)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Parameters / inputs

- rho_init -initia ldensity matrix, Hilbert space
- rho_targ -target density matrix, Hilbert space
- Output:
- b -Sorensen bound
- Note: this is an exact unitary bound; the amount reachable
- with realistically available instrumental controls
- may be smaller, see the detailed analysis here:

## Implementation structure

- Sorensen bound for the maximum transfer efficiency between
- two states under arbitrary control operators. Equation 186
- from https://doi.org/10.1016/0079-6565(89)80006-8. Syntax:
- b=sorensen(rho_init,rho_targ)
- rho_init -initia ldensity matrix, Hilbert space
- rho_targ -target density matrix, Hilbert space
- Output:
- b -Sorensen bound
- Note: this is an exact unitary bound; the amount reachable
- with realistically available instrumental controls
- may be smaller, see the detailed analysis here:
- Check consistency
