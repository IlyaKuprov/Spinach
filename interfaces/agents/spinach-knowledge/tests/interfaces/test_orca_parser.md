# tests/interfaces/test_orca_parser.m

Source: [tests/interfaces/test_orca_parser.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/interfaces/test_orca_parser.m)

## Purpose

Regression test for the ORCA log parser (`oparse`) that runs on the ORCA output files bundled with the Spinach examples. It verifies that magnetic parameters are read from ORCA logs and attached to the correct atoms.

## Behaviour

- Announces the test target with `TESTING: ORCA log parser`.
- Creates a test result object via `new_test_result` for `interfaces/orca_parser`.
- Locates the Spinach root directory three levels above the test file using `fileparts(fileparts(fileparts(mfilename('fullpath'))))`.
- Parses the methyl radical log (`examples/esr_liq_pulsed/data_import/orca_methyl_radical.out`, a vacuum DFT calculation with a g-tensor and hyperfines) and checks:
  - Version detection: `props.orca_version` equals `'4.1.1'`, confirming the parser branch is chosen from the ORCA version banner.
  - g-tensor: `props.g_tensor.matrix(1,1)` is close to `2.0022269` with absolute and relative tolerances `1e-7`, confirming the g-matrix is read from the electronic g-matrix section.
  - Hyperfine: `gauss2mhz(props.hfc.full.matrix{2})` is converted, and `trace(proton_hfc)/3` is close to `-65.1341` with tolerances `1e-4`, confirming the isotropic hyperfine coupling matches the A(iso) printed by ORCA.
- Parses the copper porphyrin log (`examples/visualisation/porphyrine.out`, where hyperfine tensors are printed for the protons only) and checks:
  - Skipped nucleus numbering: the atoms with non-empty `props.hfc.full.matrix` entries exactly match the atoms whose symbol is `'H'`, because ORCA numbers printed nuclei from zero and skips the ones it was not asked about.
  - Atom count: `props.natoms == 37`, confirming the coordinate table sets the length of every per-atom output.

## Inputs and outputs

Syntax:

```matlab
result = test_orca_parser()
```

- **Outputs**:
  - `result` — regression test result with explanatory messages.
- **Inputs**: none.

## References

- `oparse` — ORCA log parser under test.
- `new_test_result`, `test_true`, `test_close` — test harness utilities.
- `gauss2mhz` — unit conversion used for the hyperfine check.
