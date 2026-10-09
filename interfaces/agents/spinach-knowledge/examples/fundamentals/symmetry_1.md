# examples/fundamentals/symmetry_1.m

- MATLAB implementation: [examples/fundamentals/symmetry_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/symmetry_1.m)

- Signature: `symmetry_1()`
- Source: [examples/fundamentals/symmetry_1.m](../../../../../examples/fundamentals/symmetry_1.m)

## Purpose

Demonstrates S4 symmetry factorisation of a zero-field radical-pair Liouvillian and compares its sparsity with the untransformed operator.

## Spin model and basis

The source defines two electron spins followed by four protons, with `sys.magnet=0`. Protons 3–6 form one `S4` symmetry group. It builds a full spherical-tensor Liouville basis (`sphten-liouv`, approximation `none`) and sets `sym_a1g_only=0`, so the basis is not limited to the A1g sector. The electron Zeeman scalars are both 2.002; the proton entries are zero. The scalar-coupling matrix couples electron 1 to each of protons 3–6 with input value 0.295, converted by `mt2hz`; the other listed entries are zero.

## Construction and display

After creating the system and basis, the example applies the `labframe` assumption and obtains `H=hamiltonian(spin_system)`. It concatenates `bas.sym_fact(1).irr_projectors` as `S` and displays sparsity masks for `abs(H)>1e3` and `abs(S'*H*S)>1e3`. Guide lines mark indices 560, 1216, and 1936. These are plot settings, not reported dimensions or measured results.

## Scope

This is an operator-structure demonstration: the source constructs and plots the matrices, but contains no time evolution, acquisition, or numerical comparison of dynamics.
