# experiments/bruker/hsqcedetgp.m

## Purpose and sequence

This routine implements Bruker's echo/antiecho, gradient-selected multiplicity-edited HSQC sequence. Gradient selection is represented analytically by coherence-order filters, not explicit gradient pulses. The function begins with F2 longitudinal magnetisation and applies a proton 90-degree pulse. The initial INEPT transfer uses two evolution periods of `delta/2`, where `delta = abs(1/(2*J))`, separated by simultaneous F1/F2 180-degree refocusing. A proton trim pulse and transfer pulses then prepare the F1 coherence.

F1 is sampled through half-step evolution blocks with requested midpoint refocusing pulses. The resulting +1 and -1 F1 coherence branches represent the echo and antiecho pathways. Each branch then undergoes two multiplicity-editing delays of `edit_time`, separated by a simultaneous F1/F2 180-degree pulse. Back-transfer uses two more `delta/2` periods with a simultaneous 180-degree pulse between them. A final coherence filter retains zero F1 order and +1 F2 order; the code decouples the selected F2 nuclei and detects the two branches separately. The source labels the paired `edit_time` intervals as the multiplicity-editing period; it does not give a recommended delay or map output signs to specific multiplicity classes.

## Call contract

`fid=hsqcedetgp(spin_system,parameters,H,R,K)`

- `spin_system` must use the `sphten-liouv` formalism.
- `parameters.sweep` is a two-element vector of positive real F1/F2 sweep widths in Hz; the sample time increments are their reciprocals.
- `parameters.npoints` is a two-element vector of positive integer F1/F2 point counts.
- `parameters.spins` is a two-element cell array naming the active F1 and F2 nuclei (typically a heteronucleus and 1H in the source examples).
- `parameters.decouple_f1` is the cell list of nuclei given midpoint 180-degree refocusing pulses in F1. The source documentation says this list must not contain the active F1 isotope.
- `parameters.decouple_f2` is the cell list of nuclei to decouple during F2 acquisition (examples in the source are 15N and 13C).
- `parameters.J` is the working scalar coupling in Hz; validation requires a nonzero real scalar. The INEPT interval is computed as `abs(1/(2*J))` seconds and split into two equal halves.
- `parameters.trim_angle` is the proton trim-pulse angle in radians; validation requires one finite real scalar, without specifying a narrower allowed range.
- `parameters.edit_time` is the multiplicity-editing interval in seconds; validation requires one positive real scalar. The sequence uses it twice, with the refocusing pulse between the intervals.
- `H`, `R`, and `K` are numeric matrices supplied by the context function. The source requires matching dimensions and combines them as `H + 1i*R + 1i*K`.

The output is a structure containing `fid.pos` and `fid.neg`, documented as echo and antiecho signal components. They are the F2 observable-evolution outputs for the two F1 coherence branches. The routine requests `npoints(1)-1` F1 trajectory steps and `npoints(2)-1` F2 acquisition steps; both arrays have `npoints(2)` F2-time rows and `npoints(1)` F1-stack columns; the source does not prescribe a numeric class or the multiplicity-to-sign assignment.

Natural-abundance simulations should use Spinach isotope-dilution functionality, as noted in the source.

## References and source

- Original HSQC sequence: https://doi.org/10.1016/0009-2614(80)80041-8
- HSQC reference: https://doi.org/10.1002/cmr.a.10095
- Multiplicity-editing reference: https://doi.org/10.1002/mrc.1260310315
- Spinach documentation: https://spindynamics.org/wiki/index.php?title=hsqcedetgp.m
- Source: https://github.com/IlyaKuprov/Spinach/blob/main/experiments/bruker/hsqcedetgp.m
