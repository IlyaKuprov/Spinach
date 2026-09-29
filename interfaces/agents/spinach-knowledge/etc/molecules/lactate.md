# etc/molecules/lactate.m

- MATLAB implementation: [etc/molecules/lactate.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/molecules/lactate.m)

**Call:** `[sys, inter] = lactate(spins)` (for example, `[sys,inter]=lactate({'1H','CA'})`).

**Input:** `spins` must be a cell array whose elements are MATLAB character arrays. Each entry is matched against the isotope list or atom labels; examples include `{'1H','13C'}`, `{'CO','CA','HA'}`, and mixtures such as `{'1H','CA'}`. Selection is membership-based: isotope entries select every spin of that isotope, labels select the named spin, and the returned subset remains in the source's canonical order rather than the order of the selector list. Unsupported strings are not rejected; they simply select no additional spins.

## Model and parameters

The template contains seven spins in this order: `CO, CA, CB, HA, HB1, HB2, HB3`, with isotopes `13C, 13C, 13C, 1H, 1H, 1H, 1H`. The source describes this as 13C-labelled lactate with OH protons assumed to exchange rapidly with water; there is no OH proton among the seven explicit spins. Approximate shifts, in that order, are 182.0, 68.4, 20.6, 4.0, 1.2, 1.2, and 1.2.

The scalar-coupling template is also explicitly marked “very approximate”. It assigns CO–CA and CA–CB as 40.0, CA–HA as 145.0, and CB–HB1/HB2/HB3 as 128.0; the remaining populated entries are small couplings: CO–HA 3.9, CO–HB1/HB2/HB3 3.7, CO–CB 3.0, CA–HB1/HB2/HB3 4.4, CB–HA 3.6, and HA–HB1/HB2/HB3 6.9. Values are reproduced as given by the source, not presented as a precision reference dataset.

## Outputs

- `sys`: the selected isotope and label lists.
- `inter`: the corresponding selected shift vector and coupling submatrix.

Only these two outputs are returned; despite the old page's purpose text showing a three-output call, the source signature has no `bas` output. An empty or entirely unmatched selector can produce an empty selected system because the implementation does not require at least one match.

**Source reference:** [Spinach Wiki: lactate.m](https://spindynamics.org/wiki/index.php?title=lactate.m).
