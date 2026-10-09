# examples/nmr_solids/case_studies/mathies_14n_13c/powder_gly_14n.m

- MATLAB implementation: [examples/nmr_solids/case_studies/mathies_14n_13c/powder_gly_14n.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/case_studies/mathies_14n_13c/powder_gly_14n.m)

- Signature: `powder_gly_14n()`

## Purpose

Computes a static ¹⁴N powder spectrum of glycine, assuming ¹H and ¹³C decoupling. The first simulated trace uses CASTEP shielding and electric-field-gradient tensors; a second simulated trace uses the ¹⁴N quadrupolar parameters measured by O'Dell and Schurko (PCCP 2009, DOI: [10.1039/b906114b](https://doi.org/10.1039/b906114b)). The source identifies the calculation as static powder, not MAS, and notes that it uses a numerical rotating-frame transformation because the ¹⁴N quadrupolar interaction is large. Its header reports a calculation time of seconds.

## Model and calculation

- Reads `glycine.magres` with `c2spinach`, discards H, O, and C atoms, and retains the single nitrogen site as ¹⁴N.
- Converts the CASTEP shielding tensor to a shift, then sets the isotropic shift to 110.0 ppm. The CASTEP quadrupolar interaction is formed with `castep2nqi(props.efg{1},20.44e-3,1)` and made traceless with `remtrace`.
- Uses a 9.4 T field and an `sphten-liouv` basis with no approximation. Powder averaging uses `rep_2ang_12800pts_sph`; acquisition uses a 3 MHz sweep, 256 points, 1024-point zero filling, zero offset, and MHz axis units. The initial state and receiver are both the ¹⁴N raising operator.
- Simulates with `powder(...,@acquire,...,'nmr')`, applies exponential apodisation with parameter 6, and Fourier-transforms the zero-filled FID.

## Comparison trace

The overlaid second trace is also simulated, not an experimental spectrum: the source replaces the quadrupolar matrix with `2*pi*eeqq2nqi(1.18e6,0.53,1,[0 0 0])`, using the O'Dell PCCP 2009 measured quadrupolar parameters, then repeats the powder simulation and processing. The plot legend labels the two simulations `CASTEP` and `O'Dell PCCP 2009`.
