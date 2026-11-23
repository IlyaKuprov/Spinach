# tests/interfaces/test_hfc_isotopes.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/interfaces/test_hfc_isotopes.m)

## Purpose

Regression test for isotope-resolved hyperfine coupling imports from Gaussian and ORCA electronic structure logs. `g2spinach` selects elements and keeps the source isotopes and full tensors. A subsequent `isoswap` changes the selected isotope and rescales the tensor.

## What the suite checks

- In Gaussian nitroxide, the imported nitrogen-14 tensor agrees with the parsed tensor in Spinach units. Replacing nitrogen-14 by nitrogen-15 changes the sign and scales both stored pair directions. The threshold uses the imported norm: equality retains the tensor, a larger threshold removes it, and purging then removes the uncoupled nucleus.
- A missing label on a selected atom is an error. An electron-only selection still returns the electron g-tensor.
- In ORCA porphyrine, all imported nuclei have isotope labels, but only the protons have printed hyperfine tensors. Their hyperfine records identify them as hydrogen-1, not the isotope specified for the quadrupole calculation. Replacing the protons by deuterons scales every tensor component by the positive magnetogyric ratio.

## Inputs and outputs

**Syntax**

```matlab
result = test_hfc_isotopes()
```

The function takes no arguments.

**Outputs**

- `result` — regression check accumulator covering isotope labels, hyperfine tensors, thresholding, purging, and isotope substitution.

## References

- Uses `new_test_result`, `test_true`, `test_close`, `gparse`, `oparse`, `g2spinach`, `isoswap`, `gauss2mhz`, and `spin`.
