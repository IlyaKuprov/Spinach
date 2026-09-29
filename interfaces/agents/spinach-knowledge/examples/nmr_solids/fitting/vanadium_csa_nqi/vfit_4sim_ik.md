# examples/nmr_solids/fitting/vanadium_csa_nqi/vfit_4sim_ik.m

- Signature: vfit_4sim_ik()

## Purpose and data

The example jointly fits four ⁵¹V MAS NMR spectra to chemical-shielding and quadrupolar tensor parameters. It loads v12_29_dec15.spc, v12_31_dec15.spc, v12_33_dec15.spc, and v12_35_dec15.spc, applies a Savitzky–Golay filter of order 3 and window length 51, crops the signal columns to indices 6200–10000, and normalises each retained segment to its own maximum before padding to 4096 points. The source describes calculation time as hours and says it is much faster with a GPU; however, the GPU-enable line in the fitting function is commented out. Neither a timing nor GPU execution was measured here.

## Shared spin model and optimiser inputs

Each trial uses one ⁵¹V spin, a chemical-shielding tensor and a quadrupolar tensor. The code maps the fit vector to isotropic shift, anisotropy, asymmetry and Euler angles; it converts the three angles from degrees to radians, forms the shielding principal values, and adds 456.818 to those values. Quadrupolar coupling is passed through eeqq2nqi with spin 3.5. The code sets sys.magnet to 14.1 without an explicit unit comment and uses an sphten-liouv basis with no approximation and projection +1.

The optimiser starts from [-669.0, 564.0, 0.255, 82.0, 180.0, 19.0, 3.72, 0.62]. These are initial code values, not fitted results. The first two are transformed by the script when constructing the shielding tensor; the seventh is multiplied by 1e6 for the quadrupolar-coupling input. The source does not annotate units for these fit-vector values.

## Four MAS calculations and observable

The same model is simulated at code-set rate values 41000 for the 35 spectrum, 38500 for 33, 36000 for 31, and 34000 for 29; the source does not annotate their units. All four use the rep_2ang_200pts_oct powder grid, maximum rank 30, 4096 points and zero-fill 4096, with the plotted chemical-shift axis in ppm. Each FID is Gaussian-apodised with the code value 13000.0, Fourier-transformed and normalised to its maximum absolute value. The source plots input and simulated real spectra in four panels and minimises the sum of the four squared real-spectrum residual norms using fminsearch.

This is a four-spectrum fitting wrapper around singlerot and acquire, not a CP or HMQC sequence; it specifies no RF/contact-transfer condition. It describes input spectra and a model objective, but contains no saved best-fit result. No DOI is given in the source.

Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/fitting/vanadium_csa_nqi/vfit_4sim_ik.m
