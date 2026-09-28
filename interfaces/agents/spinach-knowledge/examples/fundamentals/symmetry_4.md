# examples/fundamentals/symmetry_4.m

- Signature: `symmetry_4()`

## Purpose

Builds the Hamiltonian for a radical pair with four equivalent proton nuclei, imposes S4 permutation symmetry on those nuclei, and visualises the Hamiltonian sparsity before and after symmetry factorisation.

## Physical / mathematical content

- The spin system contains two electrons and four protons; the four proton spins are grouped under the S4 symmetry group.
- The code sets zero magnetic field, equal electron Zeeman scalars of 2.002, and electron-proton scalar coupling entries of 0.295 (converted with `mt2hz`).
- It forms the Hamiltonian and concatenates the irreducible-representation projectors returned by the symmetry-adapted basis. The transformed matrix `S'*H*S` is compared with `H`.

## Numerical / algorithmic content

- `bas.sym_spins={[3 4 5 6]}` and `bas.sym_group={'S4'}` request the permutation-symmetry basis for the four equivalent nuclei.
- The final two-panel plot shows entries for which `abs(H)>1e3` and `abs(S'*H*S)>1e3`; guide lines are drawn at indices 20 and 28.

## Implementation structure

- Define the spin system, Zeeman and scalar-coupling interactions, and S4 basis.
- Create the Spinach system, build the basis under the lab-frame assumption, and construct the Hamiltonian.
- Concatenate the irrep projectors and plot the original and transformed Hamiltonian sparsity patterns.
