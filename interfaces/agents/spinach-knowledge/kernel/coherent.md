# kernel/coherent.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/coherent.m) · [Spin Dynamics Wiki: coherent.m](https://spindynamics.org/wiki/index.php?title=coherent.m)

## Contract

`rho=coherent(spin_system,mode,alpha)` prepares the density operator of a coherent state in one truncated bosonic mode. The amplitude `alpha` is a finite complex scalar; `mode` is the particle index in the system and must identify a bosonic particle of type `C`, `V`, or `T`.

## State and dimensions

Let `N=spin_system.comp.mults(mode)` be the mode's Fock truncation. Before normalisation, the state-vector coefficients for occupation `n=0,...,N-1` are `alpha^n/sqrt(n!)`. The cutoff omits the upper Poisson tail with mean `abs(alpha)^2`; the function reports its lost probability weight, then normalises the retained coefficients. Thus the constructed single-mode state has unit trace after truncation, not the untruncated infinite-dimensional coherent state.

The function forms the mode projector `|alpha><alpha|` and Kronecker-products it with identity operators for every other particle, in the particle order of the spin system. If `D=prod(spin_system.comp.mults)`, the Hilbert-space density matrix is `D-by-D`; the selected mode contributes an `N-by-N` block. This is a state operator, not a state ket. For example, `alpha=0` gives the truncated vacuum in the selected mode and leaves all other particles as identities.

## Formalism and basis

In `zeeman-hilb`, the function returns the density matrix. In `zeeman-liouv`, it returns the column vectorisation of that matrix, with `D^2` entries. Other formalisms are rejected. The construction uses the system's Zeeman product basis and its particle ordering; it does not rotate, average, or otherwise transform the state.

## Inputs

- `mode`: particle index, checked against the system's particle count and bosonic type.
- `alpha`: finite numeric scalar, with real or complex value; it is the dimensionless coherent-state amplitude.
- The mode cutoff is read from `spin_system.comp.mults(mode)`; it is not an extra function argument.
