# kernel/utilities/rspt_eig.m

- Signature: `[E,V,dE,T,LP]=rspt_eig(spin_system,parameters,Hz,Hc,Hmw,B)`

## Purpose

Eigensystem of sparse Hamiltonians to user-specified order in RSPT with careful handling of diagonal dominance and an opti- on to do exact diagonalisation (expensive). The function also returns eigenvalue derivatives and transition moments between eigenvectors under a user-specified operator. Parametrisation matches use cases in field-swept EPR spectroscopy. Syntax: [E,V,dE,T,LP]=rspt_eig(spin_system,parameters,Hz,Hc

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.

## Parameters / inputs

- Hz -laboratory frame Hamiltonian, containing only
- Zeeman terms at 1 Tesla
- Hc -laboratory frame Hamiltonian, containing all
- spin-spin couplings, but no Zeeman terms
- Hmw -observable operator without the amplitude pre-
- factor (2-norm should be around 1)
- B -magnetic field, Tesla
- parameters.rspt_order -perturbation theory order to use
- to account for the off-diagonal
- part of the Hamiltonian, Inf for
- exact diagonalisation
- parameters.rho0 -[optional] when a matrix, sets a user-
- specified thermal equilibrium state;
- when a function handle f(B,alp,bet,gam)
- sets the function to call to obtain
- the thermal equilibrium at each orien-
- tation and magnetic field; if not pro-
- vided, the thermal equilibrium is com-
- puted at the current temperature, ori-
- entation, and magnetic field

## Outputs

- E -a column vector of energies, sorted in ascen-
- ding order (rad/s)
- V -a matrix with eigenvectors in columns, sorted
- left to right in the same order as the energies
- dE -a column vector of dE/dB derivatives, sorted in
- the same order as the energies
- T -a matrix of transition moments under Hmw
- LP -a column vector of energy level populations, sor-
- ted in the same order as the energies

## Implementation structure

- Eigensystem of sparse Hamiltonians to user-specified order in
- RSPT with careful handling of diagonal dominance and an opti-
- on to do exact diagonalisation (expensive). The function also
- returns eigenvalue derivatives and transition moments between
- eigenvectors under a user-specified operator. Parametrisation
- matches use cases in field-swept EPR spectroscopy. Syntax:
- [E,V,dE,T,LP]=rspt_eig(spin_system,parameters,Hz,Hc,Hmw,B)
- Hz - laboratory frame Hamiltonian, containing only
- Zeeman terms at 1 Tesla
- Hc - laboratory frame Hamiltonian, containing all
- spin-spin couplings, but no Zeeman terms
- Hmw - observable operator without the amplitude pre-
