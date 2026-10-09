# experiments/bruker/hsqcetgp.m

## Purpose and sequence

This routine implements Bruker's echo/antiecho, gradient-selected HSQC sequence. Gradient selection is represented analytically by coherence-order filters, not explicit gradient pulses. Starting from F2 longitudinal magnetisation, the function applies a proton 90-degree pulse and performs INEPT transfer with two evolution periods of `delta/2`, where `delta = abs(1/(2*J))`, separated by simultaneous F1/F2 180-degree refocusing. A proton trim pulse and transfer pulses prepare the indirect F1 coherence; the carbon transfer pulse is phase-cycled as the difference of +90- and -90-degree applications.

F1 is sampled through half-step evolution blocks with requested midpoint refocusing pulses. The +1 and -1 F1 coherence filters form the echo and antiecho branches. After selection, the code applies the F1-spin inversion and simultaneous F1/F2 90-degree back-transfer pulse, then two further `delta/2` evolution periods with a simultaneous 180-degree pulse between them. A final filter retains zero F1 order and +1 F2 order. The selected F2 nuclei are decoupled and each branch is detected independently.

## Call contract

`fid=hsqcetgp(spin_system,parameters,H,R,K)`

- `spin_system` must use the `sphten-liouv` formalism.
- `parameters.sweep` is a two-element vector of positive real F1/F2 sweep widths in Hz; the sample time increments are their reciprocals.
- `parameters.npoints` is a two-element vector of positive integer F1/F2 point counts.
- `parameters.spins` is a two-element cell array naming the active F1 and F2 nuclei (the source examples use 13C and 1H).
- `parameters.decouple_f1` is the cell list of nuclei receiving midpoint 180-degree refocusing pulses during F1 evolution.
- `parameters.decouple_f2` is the cell list of nuclei to decouple during F2 acquisition (the source gives 15N and 13C as examples).
- `parameters.J` is the working scalar coupling in Hz; validation requires a nonzero real scalar. The code computes `delta = abs(1/(2*J))` seconds, with each INEPT transfer half lasting `delta/2`.
- `parameters.trim_angle` is the proton trim-pulse angle in radians; validation requires one finite real scalar, without specifying a narrower allowed range.
- `H`, `R`, and `K` are numeric matrices supplied by the context function. The source requires matching dimensions and combines them as `H + 1i*R + 1i*K`.

The function returns a structure with `fid.pos` and `fid.neg`, documented as echo and antiecho signal components. Each is the F2 observable-evolution output for its selected F1 branch. The routine requests `npoints(1)-1` F1 trajectory steps and `npoints(2)-1` F2 acquisition steps; both arrays have `npoints(2)` F2-time rows and `npoints(1)` F1-stack columns; the source does not prescribe a numeric class.

Natural-abundance simulations should use Spinach isotope-dilution functionality, as noted in the source.

## References and source

- Original HSQC sequence: https://doi.org/10.1016/0009-2614(80)80041-8
- HSQC reference: https://doi.org/10.1002/cmr.a.10095
- Spinach documentation: https://spindynamics.org/wiki/index.php?title=hsqcetgp.m
- Source: https://github.com/IlyaKuprov/Spinach/blob/main/experiments/bruker/hsqcetgp.m
