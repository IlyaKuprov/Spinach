# examples/nmr_solids/mas_powder_suc_floquet.m

Source: [examples/nmr_solids/mas_powder_suc_floquet.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/mas_powder_suc_floquet.m)

## Purpose

Calculates a 13C MAS NMR spectrum for the sucrose system imported from the PCM-DFT output file. The source labels the method Floquet MAS and estimates “days” of calculation time, or “hours with a Tesla card”; those are source comments, not measured timings.

## Model and acquisition

`gparse('../standard_systems/sucrose.log')` and `g2spinach` import the carbon isotope set `13C`, using `182.1` as the absolute shielding reference (ppm; `g2spinach` defines reference values as those placed at zero ppm). The source describes the shielding tensors, J-couplings, and coordinates as DFT-derived. The field is `14.1` T. The basis uses `sphten-liouv`, `IK-0`, projection `+1`, and inter-level 3; interaction and proximity cutoffs are set to 5.0 and 4.0.

The rotor-axis vector is `1 1 1`, MAS rate is `6000` Hz, maximum rank is 23, sweep is `50000` Hz, and offset is `15000` Hz. The acquisition has 256 points and zero-fills to 1024. It selects `13C`, specifies an empty `decouple` list, and labels the axis in Hz. The header assumes 1H decoupling, but the script selects only 13C, sets no explicit RF field or pulse program, and does not request a decoupling channel. It assigns no dipolar tensor value directly; the parsed model also carries molecular coordinates.

## Calculation and display

The source calls `floquet(spin_system,@acquire,parameters,'nmr')` with grid `leb_2ang_rank_23`. The initial state and receiver are both `L+` on `13C`. It exponentially apodises the calculated FID with parameter 6, Fourier transforms to 1024 points, and plots the real spectrum with `plot_1d`. These are simulated settings and output, not experimental measurements.
