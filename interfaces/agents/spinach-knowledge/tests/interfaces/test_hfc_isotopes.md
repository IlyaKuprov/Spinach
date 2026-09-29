# tests/interfaces/test_hfc_isotopes.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/interfaces/test_hfc_isotopes.m)

## Purpose

Regression test for isotope-resolved hyperfine coupling imports from Gaussian and ORCA electronic structure logs. It verifies that hyperfine tensors follow the requested nuclear gyromagnetic ratio during `g2spinach` conversion, including anisotropic and off-diagonal components, and that source isotopes are taken from the shipped logs rather than assumed from natural abundance.

## Behavior

- Registers a test named `interfaces/hfc_isotopes` with description `EPR isotope conversion` and the physical statement `hyperfine tensors follow the requested nuclear gyromagnetic ratio`.
- Parses the Gaussian log `examples/standard_systems/nitroxide.log` with `gparse` and locates the nitrogen atom by symbol.
- Confirms the Gaussian source isotope is 14N, so the printed hyperfine belongs to 14N rather than the target 15N.
- Imports the same log twice via `g2spinach`, once requesting `14N` and once `15N`, and compares the full 3x3 coupling tensor against `1e6*gauss2mhz(props.hfc.full.matrix{atom}/2)` for the same-isotope case (tolerances 0, 0).
- Uses `isoswap(sys14,inter14,1,'15N')` and checks that the direct 15N import matches isotope replacement within `1e-8` absolute and `1e-14` relative tolerance, for both stored tensor halves (`matrix{1,2}` and `matrix{2,1}`).
- Scales Gaussian proton tensors to deuterons and verifies every proton tensor scales with the positive deuterium gamma ratio `spin('2H')/spin('1H')` within `1e-8`/`1e-14`.
- Sets `options.min_hfc` to `norm(source_hfc,'fro')*(1+abs(spin('15N')/spin('14N')))/2` with `options.purge='on'` and checks that the 15N coupling grows above the threshold (retained isotopes `{'15N','E'}`) while the source 14N stays below it (purged to `{'E'}`).
- Sets `options.min_hfc=norm(source_hfc,'fro')` and checks the strict less-than boundary: a tensor exactly at the existing threshold is retained (two isotopes).
- Verifies NMR conversion is independent of hyperfine isotope provenance: importing `{{'N','15N'}}` with scalar `0` couplings gives identical system and interaction structures whether or not `props.isotopes` is present (`isequaln` on both outputs).
- Confirms that a nonempty hyperfine tensor with the `isotopes` field removed raises an error containing `not implemented` and `source isotope`.
- Confirms that a zero entry in the parser isotope metadata for the nitrogen atom is likewise rejected with the same error message, since zero does not identify a source isotope.
- Builds eight invalid provenance variants (field removed; zero entry; NaN entry; truncated isotope list; `'14C'` string in a cell array; two-row `['14N';'15N']` entry; numeric `14` entry; scalar `'14N'` string) and, for each, exercises three modes (normal input, `std_geom` removed, `g_tensor` removed). Each case must be rejected by `grumble` with a message containing `explicit HFC source isotope`, with the rejection occurring before console output or parameter processing (checked via `evalc` capturing empty output).
- Verifies that a selected but empty hyperfine tensor (`props.hfc.full.matrix{atom}=[]`) imports without an isotope field, leaving `coupling.matrix{1,2}` empty.
- Verifies an electron-only import (`{{'E','E'}}` with scalar `0`) needs no nuclear provenance and yields isotopes `{'E'}`.
- Parses the ORCA log `examples/esr_liq_pulsed/data_import/orca_methyl_radical.out` with `oparse`, finds all hydrogen atoms, and for each proton:
  - checks `props.isotopes{atoms(n)}` equals `'1H'`, i.e. ORCA HFC provenance comes from `A:ISTP`, not `Q:ISTP`;
  - compares the same-isotope proton import against `1e6*gauss2mhz(props.hfc.full.matrix{atoms(n)}/2)` with zero tolerances;
  - compares the deuteron import against the source tensor scaled by `spin('2H')/spin('1H')` within `1e-8`/`1e-14`, covering the full non-diagonal tensor.
- Imports a direct zero-spin target `12C` and checks the resulting coupling `matrix{1,2}` equals `zeros(3)` with isotopes `{'12C','E'}`, i.e. a zero-gamma target produces a zero tensor without division by its gamma.
- Sets the carbon source isotope metadata to `'12C'` and confirms a `13C` request is rejected by `grumble` with a message containing `zero-gamma HFC source isotope` before warnings, rather than divided into.
- Sets `options.min_hfc` to `source_norm*(1+spin('2H')/spin('1H'))/2` where `source_norm` is the Frobenius norm of the first proton coupling: without purge, thresholding removes the deuterium coupling (`pruned.coupling.matrix{1,end}` empty) but keeps four isotopes; with `options.purge='on'`, all three deuterium tensors fall below the threshold and only `{'E'}` remains.
- Empties one ORCA isotope entry (`unknown.isotopes{atoms(1)}=[]`) and confirms a `2H` request raises the `not implemented` / `source isotope` error, since empty metadata is not a licence to assume 1H.
- Parses the partial ORCA output `examples/visualisation/porphyrine.out` and imports `{{'E','E'},{'N','15N'},{'H','2H'}}`: unprinted nitrogen tensors stay empty without fabricated provenance.
- Adds a declared antisymmetric orbital contribution `[0 1 -2;-1 0 3;2 -3 0]/10` to one proton hyperfine tensor in the parsed data and verifies the imported coupling equals `1e6*gauss2mhz(source_hfc/2)*(spin('2H')/spin('1H'))` within `1e-8`/`1e-14`, i.e. the antisymmetric part is scaled without symmetrisation.

## Inputs and outputs

**Syntax**

```matlab
result = test_hfc_isotopes()
```

The function takes no arguments.

**Outputs**

- `result` — regression check accumulator returned by `new_test_result` and progressively updated by `test_true` and `test_close`, covering tensor scaling, provenance, thresholding, purging, and unchanged NMR imports.

## References

- Uses `new_test_result`, `test_true`, `test_close`, `gparse`, `oparse`, `g2spinach`, `isoswap`, `gauss2mhz`, and `spin`.
