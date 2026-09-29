# examples/fundamentals/symmetry_4.m

- MATLAB implementation: [examples/fundamentals/symmetry_4.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/symmetry_4.m)

- Signature: `symmetry_4()`
- Source: [examples/fundamentals/symmetry_4.m](../../../../../examples/fundamentals/symmetry_4.m)

## Purpose

Shows how S4 symmetry factorisation changes the sparsity pattern of the Hamiltonian for a zero-field radical pair with four equivalent protons.

## Spin model and basis

The source sets `sys.magnet=0` for two electrons and four protons, groups protons 3–6 under S4, and uses the Zeeman-Hilbert formalism (`zeeman-hilb`) with approximation `none`. The electron Zeeman scalars are both 2.002 and the proton entries are zero. In the scalar-coupling matrix, electron 1 couples to each proton with input value 0.295 passed through `mt2hz`; the other listed entries are zero.

## Construction and display

The system and basis are built, the `labframe` assumption is applied, and the source calls `hamiltonian(spin_system)` and concatenates the irrep projectors into `S`. The source labels this section “Hamiltonian superoperator”, while the selected basis formalism is `zeeman-hilb`; this page reports both source details without resolving that terminology. The plotted masks use `abs(H)>1e3` and `abs(S'*H*S)>1e3`, with guide lines at indices 20 and 28.

## Scope

This is a matrix-sparsity visualisation, not a time-domain simulation. The source specifies plot thresholds and guide lines but does not report a numerical comparison or observed run result.
