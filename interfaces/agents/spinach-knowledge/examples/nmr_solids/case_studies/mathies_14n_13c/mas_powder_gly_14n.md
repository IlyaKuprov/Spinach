# examples/nmr_solids/case_studies/mathies_14n_13c/powder_gly_14n.m

- Signature: `powder_gly_14n()`

## Purpose

Calculates a static 14N powder spectrum for glycine, assuming 1H and 13C decoupling. CASTEP tensors provide the nitrogen shielding and quadrupolar interaction; the calculated spectrum is compared with the quadrupolar parameters reported by O'Dell and Schurko (PCCP 2009, DOI 10.1039/b906114b). The example reports a runtime of seconds. A numerical rotating-frame transformation is used because the 14N quadrupolar interaction is large.

## Model and calculation

- Reads `glycine.magres`, removes H, O, and C, and retains the single 14N site.
- Uses a 9.4 T field, sets the isotropic shift to 110.0 ppm, and obtains the CASTEP quadrupolar interaction with `castep2nqi(props.efg{1},20.44e-3,1)`.
- Builds an `sphten-liouv` basis with no approximation. The powder calculation uses `rep_2ang_12800pts_sph`; acquisition settings are a 3 MHz sweep, 256 points, and 1024-point zero filling.
- Applies exponential apodisation with parameter 6 and Fourier transforms the FID. A second powder simulation uses `eeqq2nqi(1.18e6,0.53,1,[0 0 0])`, the comparison parameters attributed in the plot legend to O'Dell (PCCP 2009).
