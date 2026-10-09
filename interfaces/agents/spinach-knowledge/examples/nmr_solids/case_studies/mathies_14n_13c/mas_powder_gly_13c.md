# examples/nmr_solids/case_studies/mathies_14n_13c/mas_powder_gly_13c.m

Source: [mas_powder_gly_13c.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/case_studies/mathies_14n_13c/mas_powder_gly_13c.m)

## Model and input data

The source describes a glycine-powder ¹³C MAS spectrum assuming ¹H decoupling. It uses the Fokker–Planck MAS formalism and a spherical orientation grid to show the field dependence of the ¹³Cα line shape in the presence of quadrupolar ¹⁴N. The calculation uses a rotating frame for ¹³C and the laboratory frame for ¹⁴N.

The script reads CASTEP data from glycine.magres, removes H and O atoms, and retains ¹³Cα and ¹⁴N. It converts their shielding tensors into shifts and sets the isotropic shifts to 43.6 and 110.0, which the source explicitly labels experimental values; it does not state their units. The ¹⁴N quadrupolar interaction is formed from the ¹⁴N electric-field-gradient tensor with castep2nqi using 20.44 × 10⁻³ and spin 1, then passed through remtrace. Thus the calculation has structural/tensor input and two stated experimental isotropic shifts, but no experimental spectrum or measured FID as input.

## Field series and detection

The script computes fields of 4.7, 9.4, and 14.1 Tesla, as stated in its plot legend. It sets the MAS rate parameter to 10000 without stating a unit, uses the rep_2ang_200pts_sph grid and maximum rank 5, and configures 128 acquired points with zero-fill to 1,024. The sweep is 200 times the field; the source comments that this is unchanged in ppm. Numerical rotating-frame transforms are set for ¹³C harmonic 1 and ¹⁴N harmonic 2. No RF pulse or CP-contact sequence is specified; ¹H decoupling is the stated simulation assumption.

At each field, the initial state and receiver coil are ¹³C L+. A single-rotor acquisition is run in the lab-frame mode, followed by exponential apodisation with parameter 6 and a Fourier transform. The plotted real spectra are overlaid and labelled 4.7, 9.4, and 14.1 Tesla. This is a field-dependent simulation; the only explicitly identified experimental inputs are the two isotropic shift values.
