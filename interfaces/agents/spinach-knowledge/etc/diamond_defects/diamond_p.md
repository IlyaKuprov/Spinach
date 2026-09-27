# etc/diamond_defects/diamond_p.m

- Signature: [sys,inter]=diamond_p(parameters)

## Purpose

Builds spin systems for the phosphorus-related diamond centres MA1, NP1–NP6, NP8, and NP9. The tabulated magnetic parameters are from Nadolinny et al., *Crystals* **7**, 237 (2017), https://doi.org/10.3390/cryst7080237.

## Physical / mathematical content

The selected centre determines the electron g tensor and the set of nuclear hyperfine tensors. The model includes 31P nuclei and, for NP1–NP3, the specified 14N nuclei. MA1 may additionally include the reported 13C hyperfine coupling. The routine does not add zero-field splitting or nuclear quadrupole interactions.

## Numerical / algorithmic content

The tabulated principal values are converted from mT to frequency units and transformed from their principal-axis frame to the requested crystal orientation. MA1's 13C coupling is optional and disabled by default.

## Parameters / inputs

- parameters.centre: 'ma1', 'np1', 'np2', 'np3', 'np4', 'np5', 'np6', 'np8', or 'np9'.
- parameters.orientation: '111', '110', or '100'; the corresponding crystal-plane normal is aligned with the magnetic field.
- parameters.include_13c: scalar logical; include the reported 13C hyperfine coupling for MA1. Defaults to false; true is supported only for MA1.

## Outputs

- sys: Spinach system specification structure.
- inter: Spinach interaction specification structure.

## Implementation structure

The function selects the table entry, builds its nuclei and tensors, rotates them for the selected orientation, and returns the Zeeman and electron–nuclear coupling matrices.
