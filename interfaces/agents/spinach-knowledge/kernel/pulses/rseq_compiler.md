# kernel/pulses/rseq_compiler.m

- Signature: `[P,T]=rseq_compiler(spin_system,L,Sx,Sy,pulse_phi,pulse_amp,pulse_dur,element_type)`

## Purpose

Precomputes the distinct propagators needed for an R-sequence, so repeated pulse settings share entries rather than requiring the same propagator to be constructed repeatedly.

## Implementation

For `180_pulse`, the function groups the values in `pulse_phi` by unique phase and constructs one propagator per phase using the Liouvillian `L+pulse_amp*(Sx*cos(phi)+Sy*sin(phi))` and the single pulse duration. For `90270_pulse`, it alternates the two supplied segment durations across the phase list, groups unique phase-duration pairs, and constructs one propagator for each pair using the same phase-modulated Liouvillian and that pair's duration. The two supported element types are simple and composite inversion pulses, respectively.

The inputs use radians for `pulse_phi`, rad/s for scalar `pulse_amp`, and seconds for `pulse_dur`. The source checks that `L`, `Sx`, and `Sy` are numeric and that the spin operators are Hermitian; `180_pulse` requires one duration, while `90270_pulse` requires two durations and an even number of phases.

## Outputs

- `P` — cell array of propagator matrices, one for each distinct setting required by the selected element type.
- `T` — index array shaped like `pulse_phi`; each entry selects the corresponding propagator in `P`.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/rseq_compiler.m) · [Spinach wiki page](https://spindynamics.org/wiki/index.php?title=rseq_compiler.m)
