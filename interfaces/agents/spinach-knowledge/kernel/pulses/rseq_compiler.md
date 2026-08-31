# kernel/pulses/rseq_compiler.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/rseq_compiler.m`
- Signature: `[P,T]=rseq_compiler(spin_system,L,Sx,Sy,pulse_phi,...`
- Total lines: 136

## Purpose

R sequence compiler. Uses the fact that R-sequences are very repetitive to pre-compile the minimal number of pulse propa- gators. Syntax: [P,T]=rseq_compiler(spin_system,L,Sx,Sy,pulse_phi,... pulse_amp,pulse_dur,element_type)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 48-49: Check consistency; implemented by `grumble(L,Sx,Sy,pulse_phi,pulse_amp,pulse_dur,element_type)`.
- Lines 51-52: Compute propagators; implemented by `switch element_type`.
- Lines 54-55: Sequence of pi pulses; implemented by `case '180_pulse'`.
- Lines 57-58: Find unique phases; implemented by `[phi,~,T]=unique(pulse_phi)`.
- Lines 60-61: Preallocate propagators; implemented by `P=cell(numel(phi),1)`.
- Lines 63-64: Run matrix exponentials; implemented by `for n=1:numel(phi)`.
- Lines 69-70: Sequence of composite pulses; implemented by `case '90270_pulse'`.
- Lines 72-73: Alternating segment durations; implemented by `seg_durs=repmat([pulse_dur(1); pulse_dur(2)],numel(pulse_phi)/2,1)`.
- Lines 75-76: Find unique phase-duration pairs; implemented by `[phi_dur,~,T]=unique([pulse_phi(:) seg_durs],'rows')`.
- Lines 78-79: Preallocate propagators; implemented by `P=cell(size(phi_dur,1),1)`.
- Lines 81-82: Run matrix exponentials; implemented by `for n=1:size(phi_dur,1)`.
- Lines 89-90: Complain and bomb out; implemented by `error('Unknown pulse element type.')`.

### Control flow inferred from the code

- Line 52: dispatches on `element_type`; cases `'180_pulse'`, `'90270_pulse'`.
- Line 64: `for` loop over `n=1:numel(phi)`.
- Line 82: `for` loop over `n=1:size(phi_dur,1)`.

### Key state/data transformations

- Lines 58: computes `[phi,~,T]` using `[phi,~,T]=unique(pulse_phi)`.
- Lines 61: computes `P` using `P=cell(numel(phi),1)`.
- Lines 65: computes `LP` using `LP=L+pulse_amp*(Sx*cos(phi(n))+Sy*sin(phi(n)))`.
- Lines 66: computes `P{n}` using `P{n}=propagator(spin_system,LP,pulse_dur)`.
- Lines 73: computes `seg_durs` using `seg_durs=repmat([pulse_dur(1); pulse_dur(2)],numel(pulse_phi)/2,1)`.
- Lines 76: computes `[phi_dur,~,T]` using `[phi_dur,~,T]=unique([pulse_phi(:) seg_durs],'rows')`.

### Local helper functions

- Line 97: `grumble()` — `function grumble(L,Sx,Sy,pulse_phi,pulse_amp,pulse_dur,element_type)`.
  - Representative operation: `if ~ischar(element_type)`.
  - Representative operation: `error('element_type must be a character string.')`.

## Parameters / inputs

- L -background Liouvillian
- Sx,Sy -Cartesian spin operators pertaining to
- the spins affected by the pulses
- pulse_phi -the sequence of pulse phases, radians
- pulse_amp -RF nutation frequency in rad/s, a scalar
- because R-sequences are phase-modulated
- pulse_dur -duration of the pulses in the sequence
- element, a vector with the length mat-
- ching the number of pulses in the sequ-
- ence element (seconds)
- element_type -R element needs to be an inversion
- pulse; common ones are:
- '180_pulse' : simple inversion pulse
- '90270_pulse' : composite inversion pulse

## Outputs

- P -unique propagators, a cell array of matrices
- T -an index array of the same dimension as pulse_phi,
- specifying which propagator is to be used at which
- slice of the phase sequence

## Implementation structure

- R sequence compiler. Uses the fact that R-sequences are very
- repetitive to pre-compile the minimal number of pulse propa-
- gators. Syntax:
- [P,T]=rseq_compiler(spin_system,L,Sx,Sy,pulse_phi,...
- pulse_amp,pulse_dur,element_type)
- L -background Liouvillian
- Sx,Sy -Cartesian spin operators pertaining to
- the spins affected by the pulses
- pulse_phi -the sequence of pulse phases, radians
- pulse_amp -RF nutation frequency in rad/s, a scalar
- because R-sequences are phase-modulated
- pulse_dur -duration of the pulses in the sequence

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `phi()`, `propagator()`, `pulse_dur()`, `pulse_phi()`, `phi_dur()`, `ischar()`, `ishermitian()`, `isscalar()`, `any()`, `ismember()`.
