# experiments/bruker/hsqcetgpsi.m

- Signature: `fid=hsqcetgpsi(spin_system,parameters,H,R,K)`

## Purpose

A sensitivity-improved, echo/antiecho gradient-selected HSQC sequence based on the Bruker `hsqcetgpsi` program. It transfers heteronuclear coherence through the specified scalar coupling, samples indirect evolution, and returns separate echo and antiecho FIDs. The source implements gradient selection analytically by selecting F1 coherence orders; it does not simulate physical gradient waveforms.

The source cites the standard HSQC and pulse-program background: [DOI 10.1016/0009-2614(80)80041-8](https://doi.org/10.1016/0009-2614(80)80041-8), [DOI 10.1002/cmr.a.10095](https://doi.org/10.1002/cmr.a.10095), [DOI 10.1016/0022-2364(91)90036-S](https://doi.org/10.1016/0022-2364(91)90036-S), [DOI 10.1021/ja00052a088](https://doi.org/10.1021/ja00052a088), and [DOI 10.1007/BF00175254](https://doi.org/10.1007/BF00175254).

## Inputs

- `spin_system` must use the sphten-liouv formalism. `H`, `R`, and `K` are same-sized Hamiltonian, relaxation, and kinetics matrices supplied by the context function; the routine combines them as `H+1i*R+1i*K`.
- `parameters.spins` is a two-element cell array naming the active indirect (F1) and detected direct (F2) isotopes, in that order. The source example is `{'13C','1H'}`.
- `parameters.sweep` and `parameters.npoints` are two-element F1/F2 vectors: sweep widths in Hz and point counts (positive integers), respectively. The time increments are their reciprocals.
- `parameters.J` is the working scalar coupling in Hz; the sequence sets its coupling-evolution interval from `abs(1/(2*J))`.
- `parameters.decouple_f1` lists nuclei receiving midpoint 180-degree refocusing pulses during F1 evolution; it must not include the active F1 isotope. `parameters.decouple_f2` lists nuclei decoupled during F2 acquisition. Both isotope lists must refer to spins in the system. The source header gives `decouple_f2={'15N','13C'}` as an example.
- `parameters.trim_angle` is the proton trim-pulse angle in radians. `parameters.si_time` is the sensitivity-improvement delay in seconds.

## Sequence and output

The initial operator is longitudinal magnetisation on F2 and the receive state is its raising operator. After the transfer and trim rotations, the routine evolves an F1 trajectory, applies the requested midpoint refocusing pulses, and selects the +1 and -1 F1 coherence pathways into separate states. Subsequent refocusing and sensitivity-improvement delays lead to direct-dimension acquisition with F2 decoupling. The result is a structure with `fid.pos` and `fid.neg`, the echo and antiecho signals from observable-mode evolution. The source recommends using Spinach isotope-dilution functionality for natural-abundance simulations (see `dilute.m`). Both returned arrays have `npoints(2) × npoints(1)` shape: F2 direct-time samples are rows and the F1 trajectory-state stack forms columns.

## Source limits

No gradient amplitudes, durations, or hardware waveform are specified: pathway selection is by coherence filtering. The source supplies no receiver phase-cycle table or calibration procedure beyond the listed parameters.

Source implementation: https://github.com/IlyaKuprov/Spinach/blob/main/experiments/bruker/hsqcetgpsi.m
