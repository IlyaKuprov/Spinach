# examples/extremes/high_symmetry_1.m

- Signature: `high_symmetry_1()`
- Source: [`examples/extremes/high_symmetry_1.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/extremes/high_symmetry_1.m)

## Purpose and spin system

The source presents a large, highly symmetric `1H` NMR system with two tert-butyl groups, supplied by Eberhard Matern. Its isotope list comprises two `31P` spins and twenty `1H` spins (22 spin-1/2 sites total). The chemical-shift inputs include `-43.844` for both phosphorus sites and proton values of `4.090` and `1.354`; the coupling table includes `301.99` between the phosphorus sites and phosphorus–proton couplings including `-321.62`, `-19.15`, and `15.63`. Spinach interprets these nuclear scalar chemical shifts in ppm and scalar couplings in Hz. The “Magnetic induction” setting `sys.magnet=9.39798` is in tesla.

## Method and acquisition

Despite the molecular symmetry, this file deliberately uses brute-force time-domain propagation in Hilbert space: `bas.formalism='zeeman-hilb'` and `bas.approximation='none'`; no permutation-symmetry groups are declared. The observable is the proton FID/spectrum. Initial state and receiver coil are both proton `L+`, decoupling is empty, and the script invokes `liquid(spin_system,@acquire,parameters,'nmr')` rather than defining an RF-pulse sequence.

The acquisition inputs are `offset=1150`, `sweep=2400`, `npoints=4096`, and `zerofill=32768`; the output axis is explicitly ppm and inverted. Offset and sweep units are not annotated in the source. The FID is exponentially apodised with parameter 5, Fourier transformed, converted to its real part, and plotted with `plot_1d`.

## Resource note and scope

The source warns that the calculation needs 32 or more CPU cores and 128 or more GB of RAM, and estimates a runtime of hours. These are source comments rather than a portable performance guarantee. The page makes no claim about an experimental assignment, computed peak positions, or failure of any Spinach kernel; the script supplies no DOI or external literature link.
