# tests/interfaces/test_hfc_isotopes.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/interfaces/test_hfc_isotopes.m)

## Purpose

Regression test for isotope-resolved hyperfine coupling imports from Gaussian and ORCA electronic structure logs. It verifies that hyperfine tensors follow the requested nuclear gyromagnetic ratio during `g2spinach` conversion, including anisotropic and off-diagonal components, and that source isotopes are taken from the shipped logs rather than assumed from natural abundance.

## What the suite checks

- Gaussian nitrogen hyperfine tensors must retain their source value for the printed isotope and rescale, including anisotropic components and the negative gyromagnetic-ratio sign, when a different target isotope is requested. Importing the same isotope must not alter the tensor. Hyperfine thresholding and purging must apply to the isotope-adjusted strength rather than silently changing the surviving interaction.
- Nonempty imported hyperfine tensors require identifiable source-isotope provenance: missing, malformed or unprinted isotope metadata cannot be replaced by a natural-abundance guess. Ordinary NMR-only conversion must remain usable without hyperfine isotope provenance.
- ORCA proton hyperfine tensors must use isotope metadata from the hyperfine records and rescale consistently for deuterium, including asymmetric components; invalid provenance is rejected without corrupting independent valid imports. The suite checks these import invariants, not a newly computed electronic-structure result.

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
