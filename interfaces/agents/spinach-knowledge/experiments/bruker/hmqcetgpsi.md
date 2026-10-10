# experiments/bruker/hmqcetgpsi.m

## Purpose and sequence

This routine implements the sensitivity-improved echo/antiecho, gradient-selected HMQC sequence based on the Bruker hmqcetgpsi pulse program. The gradient-selection steps are represented by coherence-order filters; the function does not apply explicit gradient pulses. It starts from longitudinal magnetisation on the F2 spin, applies an F2 90-degree pulse, and evolves under the scalar-coupling interval `delta = abs(1/(2*J))`. It then applies the F1-spin inversion and excitation pulses and samples the indirect dimension in half-step evolution blocks around the requested F1 refocusing pulses. The two branches are selected at F1 coherence order +1 and -1, providing the echo and antiecho pathways.

After the first selection, the code applies the documented sensitivity-improvement transfer and pulse sequence, including three further evolution intervals of duration `delta`, phase-specific F1/F2 pulses, and a final F2 echo pulse. Before acquisition it selects zero coherence on the F1 spin and +1 on F2, applies the requested F2 decoupling, and detects each branch separately. These pulse timings and phase operators are the routine's implemented sequence; no gradient waveform or gradient strength is specified here.

## Call contract

`fid=hmqcetgpsi(spin_system,parameters,H,R,K)`

- `spin_system` must use the `sphten-liouv` formalism.
- `parameters.sweep` is a two-element vector of positive real F1/F2 sweep widths in Hz. Time increments are their reciprocals.
- `parameters.npoints` is a two-element vector of positive integer F1/F2 point counts.
- `parameters.spins` is a two-element cell array naming the active F1 and F2 nuclei (the source examples use 13C and 1H).
- `parameters.decouple_f1` is the cell list of nuclei receiving midpoint 180-degree refocusing pulses during F1 evolution. The code applies these pulses between the two half-step evolution blocks.
- `parameters.decouple_f2` is the cell list of nuclei to decouple during F2 acquisition (the source gives 15N and 13C as examples).
- `parameters.J` is the working scalar coupling in Hz; validation requires a nonzero real scalar. The sequence interval is computed as `abs(1/(2*J))` seconds.
- `H`, `R`, and `K` are numeric matrices supplied by the context function; the source requires matching dimensions and combines them as `H + 1i*R + 1i*K`.

The function returns a structure with `fid.pos` and `fid.neg`, the source-described echo and antiecho signal components. Each field is the output of the F2 observable evolution for its selected F1 branch. The implementation uses `npoints(1)-1` F1 trajectory steps and `npoints(2)-1` F2 observable-evolution steps. Both `fid.pos` and `fid.neg` have `npoints(2) × npoints(1)` shape: F2 time samples are rows and the F1 state stack forms columns.

Natural-abundance simulations should use Spinach isotope-dilution functionality, as noted in the source.

## References and source

- Standard HMQC sequence: https://doi.org/10.1016/0022-2364(83)90241-X
- Sensitivity-improved sequence: https://doi.org/10.1016/0022-2364(91)90036-S
- Spinach documentation: https://spindynamics.org/wiki/index.php?title=hmqcetgpsi.m
- Source: https://github.com/IlyaKuprov/Spinach/blob/main/experiments/bruker/hmqcetgpsi.m
