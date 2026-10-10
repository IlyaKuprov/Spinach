# experiments/nmr_liquids/hetcor.m

**Canonical MATLAB source:** [experiments/nmr_liquids/hetcor.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_liquids/hetcor.m) · [Spin Dynamics Wiki: hetcor.m](https://spindynamics.org/wiki/index.php?title=hetcor.m)

## Purpose

Magnitude-mode HETCOR acquisition with F1 indirect evolution and F2 detection. The implementation assembles the supplied Hamiltonian, relaxation, and kinetics matrices as its Liouvillian; it does not construct those operators itself.

## Input contract

Call `fid=hetcor(spin_system,parameters,H,R,K)`. `parameters.sweep=[F1 F2]` gives the two sweep widths in Hz, and `parameters.npoints=[F1 F2]` gives the sample counts in those dimensions. `parameters.spins={F1,F2}` identifies the nuclei in dimension order; the source example is `{'1H','13C'}`. `parameters.decouple` lists nuclei to decouple during detection (source example: `{'1H','15N'}`), and `parameters.J` is the working scalar coupling in Hz. The source requires the `sphten-liouv` formalism and numeric, same-sized matrix inputs `H`, `R`, and `K`.

## Sequence and acquisition

The initial F1 pulse is the source's sum of a +90° F1 x-pulse and an i-weighted +90° F1 y-pulse. F1 evolves for half an indirect increment, an F2 π pulse is applied, and the second half is refocused to form the sampled F1 evolution. The sequence then evolves for `delta_2=abs(1/(2*J))`, adds the two simultaneous F1/F2 x-pulse branches at +90° and −90°, and evolves for `delta_3=abs(1/(3*J))`. These fixed delays use the absolute value of `J`; `J` is in Hz and the delays are in seconds. The listed `parameters.decouple` channels are applied before F2 observation with the F2 `L+` detection state. This describes the implemented pulse and coherence-selection operations without assigning an additional pathway not specified by the source.

## Output

`fid` is a two-dimensional, magnitude-mode free-induction decay: it contains the configured `npoints(1)` F1 evolution samples and `npoints(2)` F2 detection samples. For natural-abundance experiments, the source directs users to isotope dilution via `dilute.m`.

## Reference

- [HETCOR sequence reference](https://doi.org/10.1016/0022-2364(81)90272-9)
