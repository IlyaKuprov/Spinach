# VO(TPP) energy levels

`[energies,fields]=gcce_votpp()` calculates the 16 energy levels of the electron–vanadium-51 register over 0–0.15 T, using Eq. (10) and Table I of [Ryan et al.](https://doi.org/10.1063/5.0278659) to reconstruct Figure 2. The energy array has one row per entry in the field row vector. Energies are in GHz, and fields are in tesla.

The applied field is transverse to the unique axis of the electronic g, electron–vanadium hyperfine, and vanadium quadrupolar tensors. The interaction (p I_{\parallel}^2) is represented by its traceless part, which changes no transition frequency. The field dependence is linear at the Hamiltonian level, not at the level of individual eigenvalues.

This is the observed (empty-cluster) Hamiltonian used in gCCE. It does not calculate a Hahn echo, include a nuclear or electronic bath, or reproduce the crystal decoherence panels in Figure 3. The separate [source-input regression](../../../../../tests/kernel/test_cce_lunghi_fixtures.m) checks a reduced finite-bath echo and distinguishes Spinach's unitary pulse from the pulse in the archived PyCCE fixture.

[MATLAB source](../../../../../examples/cce/gcce_votpp.m) · [Independent spectrum regression](../../../../../tests/kernel/test_cce_examples.m)
