# examples/nmr_solids/mas_powder_suc_gridfree.m

Source: [examples/nmr_solids/mas_powder_suc_gridfree.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/mas_powder_suc_gridfree.m)

## Purpose

Calculates a 13C MAS NMR spectrum for the sucrose system imported from PCM-DFT output. The source describes grid-free Fokker-Planck MAS with a polyadic representation of the evolution generator and cites [Science Advances, DOI 10.1126/sciadv.aaw8962](https://doi.org/10.1126/sciadv.aaw8962). Its estimate of hours on a Tesla V100 GPU and much longer on CPU is a source comment, not a measured timing.

## Model and acquisition

`gparse('../standard_systems/sucrose.log')` and `g2spinach` import the `13C` system using `182.1` as the absolute shielding reference (ppm; `g2spinach` defines reference values as those placed at zero ppm). The source describes the shielding tensors, J-couplings, and coordinates as DFT-derived. The field is `14.1` T. The basis uses `sphten-liouv`, `IK-0`, projection `+1`, and inter-level 3; interaction and proximity cutoffs are 5.0 and 4.0. The source enables `greedy` and `polyadic`; the listed `gpu` option is commented out.

The rotor-axis vector is `1 1 1`, MAS rate is `6000` Hz, maximum rank is 23, sweep is `50000` Hz, and offset is `15000` Hz. The acquisition has 256 points and zero-fills to 1024. It selects `13C`, sets an empty `decouple` list, labels the axis in ppm, and requests axis inversion. The `g2spinach` import selects only carbon atoms, omitting proton spins and C–H couplings; this implements the header’s ideal `1H`-decoupling assumption. The empty `decouple` list adds no RF decoupling sequence to the carbon-only model. It assigns no dipolar tensor value directly; the parsed model also carries molecular coordinates.

## Calculation and display

The executable call is `gridfree(spin_system,@acquire,parameters,'nmr')`, matching the grid-free method named in the source header. The initial state and receiver are both `L+` on `13C`. The code exponentially apodises the calculated FID with parameter 6, Fourier transforms to 1024 points, and plots the real spectrum with `plot_1d`. This is a simulated spectrum, not an experimental measurement.
