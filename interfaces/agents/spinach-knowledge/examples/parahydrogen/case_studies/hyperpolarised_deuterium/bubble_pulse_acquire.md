# examples/parahydrogen/case_studies/hyperpolarised_deuterium/bubble_pulse_acquire.m

- Signature: bubble_pulse_acquire()
- Source: [examples/parahydrogen/case_studies/hyperpolarised_deuterium/bubble_pulse_acquire.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/parahydrogen/case_studies/hyperpolarised_deuterium/bubble_pulse_acquire.m)

## Model and assumptions

The script predicts a partially negative-line (PNL) deuterium spectrum for ortho-D2 in a parahydrogenation catalyst model. Its four spin-1 deuterons are arranged as two exchanging pairs: spins 1–2 have the D2-like shifts 4.55 and 4.55, while spins 3–4 have shifts −13.5 and −16.5. The scalar couplings are 12.0 Hz within the first pair and 0.24 Hz within the second. The field is set to 7.05 T. These are model inputs. The four-spin system contains only deuterons: no para-H2 protons or explicit SABRE transfer step is represented. The plotted spectrum is simulated, not a report of a measured PNL signal or measured hyperpolarisation.

The two chemical compartments are defined by atom groups [1 2] and [3 4], with concentrations [1 0] and two directed first-order records: pair 1 to pair 2 at 1 inverse second, reverse at 5000 inverse seconds, with spin matching [1 3; 2 4] and its inverse. The catalyst-associated deuterons also have DFT-derived, rotated nuclear-quadrupole-interaction tensors. In the source's stated 10^3 multiplier convention, the tensors are [[106.7, −6.2, 31.9], [−6.2, −55.5, 3.3], [31.9, 3.3, −51.2]] and [[−53.7, 11.8, −19.7], [11.8, −5.3, −73.8], [−19.7, −73.8, 59.0]]. Their units are not written beside these matrices. DFT coordinates for spins 3 and 4 are [−1.962, 0.573, −0.576] and [−0.175, 1.399, −1.630]; coordinates for the D2 pair are absent.

## Relaxation, bubbling, and readout

The `sphten-liouv` model uses no basis approximation and combines the spin Hamiltonian, secular Redfield relaxation, chemical exchange, and a magnetic-pumping term. Relaxation is configured with correlation times 1 ps and 400 ps, R1 values [0.04, 0.04, 0, 0], and R2 values [8, 8, 0, 0]; the equilibrium state is zero and the retained relaxation is secular. The initial state is the unit state. The pumping target is the singlet plus five quintet components, with the identity component removed; the `magpump` rate is 0.1. The source explicitly calls the bubbling-rate guess in need of a proper rate, so this term is a provisional model choice rather than a fitted experimental parameter. The acquisition example evolves this model for 7 seconds before applying a 45-degree deuterium pulse.

The simulated `hp_acquire` acquisition uses a deuterium coil, offset 209.6554 Hz, sweep 60 Hz, 256 points, and 1024-point zero filling. The ppm-axis spectrum is exponentially apodised (parameter 6), Fourier transformed, and normalised to unit maximum magnitude; the displayed region is 4.40–4.70 ppm and intensity is labelled in arbitrary units. The source says a paper link will follow; it supplies no DOI or experimental spectrum.

Zero track elimination is explicitly enabled with `zte` in `sys.enable`.

All local unit coordinates are excluded from the pumping shape and traceless detection coils. Initial units carry [1, 0]; chemical transfer evolves these populations, and magpump scales each substance source by its instantaneous population.

The former shared-unit pumping source remained constant as D2 bound to catalyst. The concentration-weighted source instead follows the free-D2 population (equilibrium 5000/5001 of the total). This is an intentional physical change, not exact parity with that constant-source trajectory.

## Weak-signal numerical sensitivity

A constant-source comparison separates the population-dependent pumping change from floating-point sensitivity. Its prepared traceless state differs from the archived state by 2.93223e-13 relative; the mapped acquisition generator, receiver, and pulse generator agree exactly. The detected signal is much smaller than the prepared state (maximum FID magnitude 1.47252e-8). Using one common dense exponential and arithmetic ordering for both preparations gives a 2.77915e-8 relative FID difference. Reversing only the coordinate order of the same physical model gives a 4.13135e-8 relative FID change. Thus the residual 3.41135e-8 after constant-source restoration is within the measured double-precision sensitivity of this weak observable, not evidence of a different acquisition generator. This numerical effect is separate from the approximately 2e-4 instantaneous-population pumping effect.
