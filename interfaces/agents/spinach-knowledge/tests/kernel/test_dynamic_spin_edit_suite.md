# tests/kernel/test_dynamic_spin_edit_suite.m

Source: [tests/kernel/test_dynamic_spin_edit_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_spin_edit_suite.m)

## Purpose

Regression test suite for the deterministic spin-system editing support utilities in the Spinach kernel: `kill_spin`, `dilute`, `dictum`, and `merge_inp`. The suite verifies that local spin-system editing helpers update dependent metadata without touching external state.

## Behavior

The test announces its target with `TESTING: Spin-system editing utilities`, registers a result under the identifier `kernel/dynamic_spin_edit_suite`, and builds a three-spin fixture (`1H`, `13C`, `13C` with labels `h`, `c1`, `c2`) via `local_edit_spin_system`.

For `kill_spin` with numeric index 2, the suite checks:

- `comp.nspins` decrements to 2 and the isotope and label lists become `{'1H','13C'}` and `{'h','c2'}`.
- `inter.coordinates` shrinks to 2 entries, and `inter.coupling.matrix` and `inter.proxmatrix` shrink to 2-by-2.
- `chem.parts` is renumbered to `{[1 2]}`.
- `comp.types` shrinks to `{'S','S'}`.
- `comp.iso_hash` is recomputed as `md5_hash({'1H','13C'})`.
- `rlx.srsk_sources` is renumbered to `[1 2]`.

A multi-part fixture with `chem.parts` equal to `{[1 2],3}` must renumber to `{[1 2],zeros(1,0)}` when spin 3 is removed. A stale-metadata fixture must have `bas`, `inter.conmatrix`, `comp.sym_group`, `inter.assumptions`, `inter.zeeman.strength`, `inter.giant.strength`, and `inter.coupling.strength` destroyed on spin removal. A logical mask `[false true false]` must remove the same spin as the numeric index.

For `dilute(spin_system,'13C',1)`, two labelled `13C` spins must produce two two-spin subsystems, each containing one `13C`.

For `dictum`, numeric calls `dictum(assumed,1,'strong')` and `dictum(assumed,[1 2],'weak')` must replace the Zeeman strength at spin 1 with `'strong'` and the coupling strength at position (1,2) with `'weak'`; an isotope call `dictum(assumed,{'13C'},'secular')` must replace the matching Zeeman strength at spin 2.

For `merge_inp`, the suite checks:

- `sys.magnet` stays 14.1 while isotopes concatenate to `{'1H','13C','15N'}` and labels to `{'h','c','n'}`.
- `inter.coordinates` concatenates in subsystem order to 3 entries.
- Identical non-extensive fields are kept: `temperature` 298, `relaxation` `{'redfield'}`, `sys.enable` `{'greedy'}`.
- Square fields block-merge with empty-cell or zero numeric padding: `coupling.scalar` becomes 3-by-3 with `{2,3}` equal to 5.0 and `{1,3}` empty; `srfk_mdepth{3,2}` equals 7.0 with `{1,2}` empty; `weiz_r1d` becomes `blkdiag(0.1,[0 0.2; 0.2 0])`.
- Per-spin rate arrays merge into columns: `r1_rates` `{0.1;0.2;0.3}` and `r2_rates` `{1.1;1.2;1.3}`.
- Spin and subsystem indices are offset by preceding spin counts: `srsk_sources` `[1 3]`, `ignore` `{[2 3]}`, `chem.parts` `{1,[2 3]}`, `chem.rates` `zeros(2)`, `chem.concs` `[1 1]`.
- Susceptibility centre lists concatenate: `suscept.chi` `{0.01*eye(3)}`, `suscept.xyz` `{[5 5 5]}`.
- Column-oriented `chem.parts` given as `{1;2}` merges into offset row lists `{1,2,3}`.
- Partless chemistry (no species split, with `rp_theory` `'haberkorn'`, `rp_electrons` 1, `rp_rates` `[1e6 2e6]`, and `tau_c` `{1e-9}`) keeps `tau_c` common and offsets electron indices to `[1 2]`.

Refusal cases verified through `local_merge_fails` (which reports success when `merge_inp` throws an error):

- Differing non-extensive values (`temperature` 300 in one subsystem).
- A nested group (`chem`) present in only some subsystems.
- An unknown `inter` subfield (`bogus`).
- An unknown subfield inside a nested group (`zeeman.bogus`).
- A `sys` structure without isotope lists.

## Inputs and outputs

```matlab
result = test_dynamic_spin_edit_suite()
```

- **Output**: `result` — regression test result with explanatory messages, built with `new_test_result` and extended by `test_true` assertions.
- The function takes no inputs.

## References

- Tested functions: `kill_spin`, `dilute`, `dictum`, `merge_inp`.
- Test infrastructure: `new_test_result`, `test_true`.
- Support utilities referenced by the fixture: `spin`, `md5_hash`.
