# kernel/utilities/rspt_eig.m

- Signature: `[E,V,dE,T,LP]=rspt_eig(spin_system,parameters,Hz,Hc,Hmw,B)`

## Purpose

Computes the eigensystem of a sparse Hamiltonian to a specified order of Rayleigh–Schrödinger perturbation theory (RSPT), with handling of diagonal dominance, or by exact diagonalisation. The parameterisation supports field-swept EPR spectroscopy. Optional outputs provide energy derivatives, transition moments and level populations.

## Parameters / inputs

- `Hz` — laboratory-frame Hamiltonian containing only Zeeman terms at 1 Tesla.
- `Hc` — laboratory-frame Hamiltonian containing spin–spin couplings but no Zeeman terms.
- `Hmw` — observable operator without its amplitude prefactor; its 2-norm should be around 1.
- `B` — magnetic field in Tesla.
- `parameters.rspt_order` — order used to account for the off-diagonal Hamiltonian. Supported perturbative orders are `1`, `2`, `3` and `4`; `Inf` selects exact diagonalisation, which is expensive.
- `parameters.rho0` — optional equilibrium-state matrix, or a function handle `f(B,alp,bet,gam)` evaluated at the magnetic field and the three entries of `parameters.orientation`. If omitted, equilibrium is computed for the current temperature, orientation and magnetic field.

## Outputs

- `E` — column vector of energies in ascending order, in rad/s.
- `V` — eigenvectors in columns, ordered to match `E`.
- `dE` — column vector of `dE/dB` derivatives, ordered to match `E`.
- `T` — matrix of squared magnitudes of transition moments under `Hmw`: `abs(V'*Hmw*V).^2`.
- `LP` — column vector of level populations, ordered to match `E`.

## Numerical / algorithmic details

The Hamiltonian is `B*Hz+Hc` and is symmetrized before calculation. For perturbative orders, its diagonal and off-diagonal parts are passed separately to `rspert`; exact diagonalisation uses `eig`. If irreducible-representation projectors are present, both Hamiltonian components and the observable are projected and symmetrized before recursive eigensystem calculations. The resulting eigenvectors are projected back to the original basis, and all energies and eigenvectors are sorted together. Optional `dE`, `T` and `LP` outputs are then calculated in the original basis. Derivatives use the Hellmann–Feynman expression `real(diag(V'*Hz*V))`; populations use `real(diag(V'*rho0*V))`.

`parameters.rspt_order` is required. The Hamiltonian components and observable must be numeric square matrices, and `B` must be a real scalar.

Contact: ilya.kuprov@weizmann.ac.il

<https://spindynamics.org/wiki/index.php?title=rspt_eig.m>