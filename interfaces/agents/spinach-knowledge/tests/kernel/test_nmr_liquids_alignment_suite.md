# tests/kernel/test_nmr_liquids_alignment_suite.m

## Purpose

Regression test that verifies compact literature-alignment probes for liquid-state NMR pulse sequences in Spinach. It exercises the gCOSY, HMBC, HSQC, CT-HSQC, NOESY-HSQC, and TOCSY pulse sequence paths that were updated during the `nmr_liquids` literature-alignment pass, checking that each produces finite, non-zero compact FIDs and that pathway selection behaves as documented.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_nmr_liquids_alignment_suite.m>

## Behaviour

The function announces the test target with `fprintf`, then creates a test result object via `new_test_result` under the identifier `kernel/nmr_liquids_alignment_suite`, with the description "NMR liquids literature-alignment probes" and the criterion that updated liquid-state pulse sequence paths must produce finite non-zero compact FIDs.

For gCOSY, it builds a two-spin `1H`/`1H` system (9.4 T magnet, scalar coupling of 8 between spins 1 and 2, `sphten-liouv` formalism, `none` approximation) using `test_spin_system`. With sweep width 100, `npoints` `[4 4]`, `1H` spins, pulse angle `pi/2`, gradient amplitude 3, gradient duration `1e-3`, gradient stabilisation delay 0, and spin-lock length 1.5, it runs `liquid(spin_system,@gcosy,parameters,'nmr')` three times with `parameters.pathway` set to `'P'`, `'N'`, and `'P+N'`. It checks that the P and N FIDs are finite, non-zero, and distinct from each other; that the P+N result is a struct with `pos` and `neg` fields; and that the `pos` and `neg` branches match the direct P and N acquisitions to tolerances of `1e-10` (absolute and relative). It then recombines the pathways for phase-sensitive processing by applying `fft` along dimension 1 with `fftshift` to `fid_pn.pos` and `fid_pn.neg`, forming `fid_states = f1_pos + conj(f1_neg)`, and transforming along dimension 2 to obtain `spec_pn`, which must be finite with non-zero real part.

For HMBC and HSQC, it builds a `1H`/`15N` heteronuclear system with scalar coupling 90. HMBC uses sweep `[100 100]`, `npoints` `[4 4]`, spins `{'15N','1H'}`, `J` of 90, and `delta_b` of `5e-3`; the test checks the FID is finite and non-zero, with the message that HMBC coherence selection must follow `parameters.spins{1}` (i.e., a non-carbon heteronucleus). HSQC adds `decouple_f1` of `{'1H'}` and `decouple_f2` of `{'15N'}`; both the positive (`pos`) and negative (`neg`) pathway components must be finite and non-zero. CT-HSQC runs with the same parameters, and each pathway component must additionally have size `[parameters.npoints(2) parameters.npoints(1)]`.

For NOESY-HSQC, it sets sweep `[100 100 100]`, `npoints` `[2 2 2]`, spins `{'1H','15N','1H'}`, `J` of 90, and `tmix` of `1e-3`, then runs `liquid(spin_system,@noesyhsqc,parameters,'nmr')`. All four pathway components (`pos_pos`, `pos_neg`, `neg_pos`, `neg_neg`) must be finite.

For TOCSY, it builds a single-spin `1H` system, prepares the initial state with `state(spin_system,'Lz','1H','cheap')`, and constructs empty Hamiltonian `H`, relaxation superoperator `R = -speye(dim)`, and kinetics `K` matrices of the Liouville-space dimension. With sweep `[100 100]`, `npoints` `[4 4]`, `1H` spins, `lamp` of 1000, and `rho0` set from the prepared state, it calls `tocsy(spin_system,parameters,H,R,K)` twice: once with `tmix` of 0 (reference) and once with `tmix` of 0.2. The test asserts that the norm of the concatenated `cos` and `sin` components with relaxation is strictly less than without, verifying that relaxation acts during the spin-lock mixing interval.

All checks are recorded through `test_close`, which appends pass/fail results with explanatory messages to the returned result object.

## Inputs and outputs

```matlab
result = test_nmr_liquids_alignment_suite()
```

**Inputs**

None. The function takes no arguments; all spin systems and pulse sequence parameters are constructed internally.

**Outputs**

- `result` — regression test result object with explanatory messages, as produced by `new_test_result` and accumulated by `test_close`.

## References

- Source file: <https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_nmr_liquids_alignment_suite.m>
