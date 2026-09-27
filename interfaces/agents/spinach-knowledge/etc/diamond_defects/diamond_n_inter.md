# etc/diamond_defects/diamond_n_inter.m

- Signature: [sys,inter]=diamond_n_inter(parameters)

## Purpose

Builds the electron–nitrogen spin system for a nitrogen interstitial in diamond, using the War9 or War10 centre parameters reported by Felton et al., *J. Phys.: Condens. Matter* **21**, 364212 (2009), https://doi.org/10.1088/0953-8984/21/36/364212.

## Physical / mathematical content

The electron g tensor and nitrogen hyperfine tensor are assembled in the crystal principal-axis frame and rotated so the requested crystal direction is aligned with the magnetic field. War9 uses the reported electron and nitrogen tensors. War10 uses the reported electron tensor and the nitrogen tensor defined in the source; the 14N coupling is obtained by scaling the 15N tensor by the nuclear gyromagnetic-ratio ratio. No zero-field splitting or nuclear quadrupole tensor is added.

## Numerical / algorithmic content

The routine constructs an orthonormal crystal frame, selects the War9 or War10 parameter set, builds the requested crystal-to-field rotation, and returns Zeeman and electron–nuclear coupling matrices in Spinach format. Unsupported centre, isotope, or orientation values raise an error.

## Parameters / inputs

- parameters.centre: 'war9' or 'war10'.
- parameters.orientation: '111', '110', or '100'; the corresponding crystal-plane normal is aligned with the magnetic field.
- parameters.nitrogen: '14N' or '15N'.

## Outputs

- sys: Spinach system specification structure.
- inter: Spinach interaction specification structure.

## Implementation structure

The function validates the single input structure, selects the centre-specific tensors, applies the isotope scaling when requested, rotates the tensors for the selected orientation, and populates the Zeeman and coupling matrices.
