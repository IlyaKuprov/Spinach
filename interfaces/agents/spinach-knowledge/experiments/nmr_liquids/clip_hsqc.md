# experiments/nmr_liquids/clip_hsqc.m

- MATLAB source: [experiments/nmr_liquids/clip_hsqc.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_liquids/clip_hsqc.m)
- Spinach Wiki: [clip_hsqc.m](https://spindynamics.org/wiki/index.php?title=clip_hsqc.m)
- Sequence reference: [DOI 10.1016/j.jmr.2008.03.009](https://doi.org/10.1016/j.jmr.2008.03.009)

## Purpose

This is a liquid-state CLIP-HSQC sequence implementation. It builds a two-dimensional, States-quadrature free-induction signal using forward density-state evolution and backward coil-state evolution. Evolution uses the composed Liouvillian `L=H+1i*R+1i*K`; the receiver branch is propagated with `L'` and negative time increments. This is the code's parameterised sequence design; it does not by itself establish an experimental or run-verified result.

## Inputs and parameters

Signature: fid=clip_hsqc(spin_system,parameters,H,R,K)

- parameters.sweep: two positive sweep widths [F1 F2] in Hz; the time increments are their reciprocals, in seconds.
- parameters.npoints: two positive integer point counts [F1 F2].
- parameters.spins: two isotope labels {F1 F2}, e.g. {'13C','1H'}; the code uses element 1 for F1 and element 2 for F2.
- parameters.J: one active scalar coupling in Hz. The code derives delta=abs(1/(2*J)), a time in seconds.
- H, R, and K: same-size Hamiltonian, relaxation, and kinetics matrices supplied by the context function. The implementation requires the sphten-liouv formalism.

## Coherence transfer and detection

The initial state is longitudinal magnetisation on spin F2 and the receiver is its L+ state. The sequence starts with a 90-degree F2 x pulse, evolves for delta/2, applies a simultaneous 180-degree x pulse on F1 and F2, and evolves for another delta/2. Subsequent F2 y and signed F1 x pulses create the two F1 branches. During indirect evolution, the code uses 1/sweep(1) spacing and refocuses with F2 pulses.

The forward states are coherence-selected for F2 order 0 together with F1 order +1 or -1. The backward-evolved coil branches are selected for F1 order 0 and F2 order +1. Inner products combine the paired branches into fid.pos and fid.neg, the two States-quadrature components over the F1 and F2 sampling dimensions.

No MATLAB execution or experimental signal is claimed here.
