# VO(TPP) energy levels and proton-limited coherence

`[energies,fields,echo,times]=gcce_votpp()` calculates the 16 energy levels of the electron–vanadium-51 register over 0–0.15 T, using Eq. (10) and Table I of [Ryan et al.](https://doi.org/10.1063/5.0278659) to reconstruct Figure 2. The energy array has one row per entry in the field row vector. Energies are in GHz, and fields are in tesla.

The applied field is transverse to the unique axis of the electronic g, electron–vanadium hyperfine, and vanadium quadrupolar tensors. The interaction (p I_{\parallel}^2) is represented by its traceless part, which changes no transition frequency. The field dependence is linear at the Hamiltonian level, not at the level of individual eigenvalues.

The middle panel of Figure 3 is calculated with projected order-two CCE at 0.33 T. The conditional hyperfine fields are computed from the 4th and 12th register eigenstates, the proton bath is unpolarised, and the ideal Hahn pulse interchanges the two levels. The full dipolar proton–proton tensors are retained for 28,192 pairs within 8 Å among 1,064 protons within 20 Å. The returned time row spans 0–30 µs in seconds, and the echo is complex. Register-changing bath matrix elements are neglected in this high-field projection; the paper used generalised CCE. The electron and deuteron bath curves are not included. The separate [source-input regression](../../../../../tests/kernel/test_cce_lunghi_fixtures.m) distinguishes Spinach’s unitary pulse from the pulse in the archived PyCCE fixture.

[MATLAB source](../../../../../examples/cce/gcce_votpp.m) · [Independent spectrum regression](../../../../../tests/kernel/test_cce_examples.m)
