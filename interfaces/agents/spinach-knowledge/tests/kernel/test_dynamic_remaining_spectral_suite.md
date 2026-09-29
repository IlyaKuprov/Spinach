# tests/kernel/test_dynamic_remaining_spectral_suite.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_remaining_spectral_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_remaining_spectral_suite.m)

## Purpose

Regression test suite for the remaining spectral, symmetry, and Fokker-Planck utilities in the Spinach kernel. It verifies field-swept eigensystem helpers, rotor-stack assembly, permutation-symmetry projectors, and one-dimensional Fokker-Planck operators against compact analytical matrix references.

## Behaviour

The suite runs the following checks:

1. **`rspt_eig` exact diagonalisation path** — with `rspt_order` set to `Inf`, compares energies, Hellmann-Feynman derivatives `dE/dB` (Zeeman operator expectation values), transition moments (squared microwave-operator matrix elements), and level populations (equilibrium density-matrix expectations) against dense Hermitian diagonalisation via `eig()`, at tolerances of `1e-14`.
2. **`eigenfields` Liouville-space resonance-field extraction** — for a diagonal pencil with `mw_freq = 100`, `window = [9 11]`, `fwhm = 0.01`, and `Iz = 2*pi*diag([10 20])`, checks that a 100 Hz transition under a 10 Hz/T Liouville pencil occurs at 10 T, the selected normalised dyadic has unit transition moment, the returned width equals the requested FWHM, population differences are unit placeholders, and the field-sweep Jacobian is scaled by the free-electron angular gyromagnetic ratio `abs(spin('E'))/(2*pi*10)`.
3. **`eigenfields` two-root Hilbert-space resonance extraction** — in a curved level gap with `mw_freq = 1/(2*pi)`, `window = [-0.8 0.8]`, and `pp_tol = 1e-20`, checks that a symmetric avoided crossing produces both resonance fields at `±sqrt(0.25-0.01)`, that multiple roots of the same level pair receive distinct branch identities (`[1 2 1;1 2 2]`), and that symmetric roots carry the same absolute field-sweep Jacobian `abs(spin('E'))./(4*root_field)`.
4. **`g2fplanck` one-dimensional gradient operator** — with `dims = 2`, `npts = 3`, checks that the 1D gradient operator is the centred spatial coordinate Kronecker product with the Zeeman Hamiltonian per tesla, and that inactive spatial dimensions return empty gradient operators.
5. **`v2fplanck` one-dimensional velocity generator** — with `dims = 1`, `npts = 10`, `deriv = {'period',3}`, and a uniform velocity field `u = 0.2`, checks that the assembled flow generator matches the finite-difference derivative construction.
6. **`rotor_stack` rank-zero assembly** — on a 14.1 T one-proton system with `max_rank = 0`, `offset = 0`, `rframes = {}`, and `masframe = 'rotor'`, checks that the returned Hamiltonian equals the direct NMR Hamiltonian (with frequency offset applied) and that the rotor phase grid is a single zero phase.
7. **`symmetry` S2 fully symmetric projector** — on a four-state two-proton product basis with `sym_group = {'S2'}`, `sym_spins = {[1 2]}`, and `sym_a1g_only = true`, checks projector orthonormality (`projector'*projector = eye(3)`), that the mixed S2 orbit symmetrises the two exchanged basis states to `[0;1;1;0]/sqrt(2)`, and that the A1g irrep dimension is 3.

## Inputs and outputs

**Syntax:**

```matlab
result = test_dynamic_remaining_spectral_suite()
```

**Outputs:**

- `result` — regression test result structure with explanatory messages, populated via `new_test_result`, `test_close`, and `test_true`.

The function takes no inputs.

## References

1. Spinach source: `tests/kernel/test_dynamic_remaining_spectral_suite.m` ([GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_remaining_spectral_suite.m))
2. Functions exercised: `rspt_eig`, `eigenfields`, `g2fplanck`, `v2fplanck`, `rotor_stack`, `symmetry`, `orientation`, `spin`, `hamiltonian`, `assume`, `frqoffset`, `clean_up`, `fdmat`, `inflate`, `spdiags`, `kron`, `parpool`, `gcp`.
