# kernel/pulses/rsequence.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/rsequence.m`
- Signature: `[phases,pulse_amp,pulse_dur]=rsequence(n_rotor_periods,n_blocks_per_period,...`
- Total lines: 197

## Purpose

R-sequences described in Malcolm Levitt's review: Nomenclature is based on the following notation RN_{n}^{\nu}. Syntax: [phases,pulse_amp,pulse_dur]=rsequence(n_rotor_periods,n_blocks_per_period,... phase_factor,n_cycle_repeats,mas_rate,... element_type,supercycle_type)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 73-75: Check consistency; implemented by `grumble(n_rotor_periods,n_blocks_per_period,phase_factor, n_cycle_repeats,mas_rate,element_type,supercycle_type)`.
- Lines 77-78: Get rotor period (seconds); implemented by `rotor_period=1/mas_rate`.
- Lines 80-81: Get R element (either simple or composite) duration (seconds); implemented by `r_element_dur=n_rotor_periods*rotor_period/n_blocks_per_period`.
- Lines 83-84: Decide sequence element; implemented by `switch element_type`.
- Lines 86-87: A train of pi pulses; implemented by `case '180_pulse'`.
- Lines 89-90: Get RF amplitude (rad/s); implemented by `pulse_amp=pi/r_element_dur`.
- Lines 92-93: One duration (seconds); implemented by `pulse_dur=r_element_dur`.
- Lines 95-96: Calculate phase grid; implemented by `phases=zeros(n_blocks_per_period,1)`.
- Lines 101-103: A train of composite pulses with the overall duration of a 360 pulse; implemented by `case '90270_pulse'`.
- Lines 105-106: Get RF amplitude (rad/s); implemented by `pulse_amp=2*pi/r_element_dur`.
- Lines 108-109: Two durations (seconds); implemented by `pulse_dur(1)=(1/4)*r_element_dur`.
- Lines 112-113: Calculate phase grid; implemented by `phases=zeros(n_blocks_per_period*2,1)`.
- Lines 121-122: Complain and bomb out; implemented by `error('unknown sequence element type.')`.
- Lines 126-127: Apply the supercycle; implemented by `switch supercycle_type`.
- Lines 131-132: Replicate with inverted phases; implemented by `phases=[phases; -phases]`.
- Lines 134-135: Replicate with shifted phases; implemented by `phases=[phases; phases+2*pi/3; phases+4*pi/3]`.
- Lines 139-141: Replicate with inverted phases and a two-step pi supercycle; implemented by `phases=[phases; -phases; -phases+pi; phases+pi]`.
- Lines 150-151: Complain and bomb out; implemented by `error('unknown supercycle type.')`.

### Control flow inferred from the code

- Line 84: dispatches on `element_type`; cases `'180_pulse'`, `'90270_pulse'`.
- Line 97: `for` loop over `q=0:(n_blocks_per_period-1)`.
- Line 114: `for` loop over `q=0:(n_blocks_per_period-1)`.
- Line 127: dispatches on `supercycle_type`; cases `'hetero_single_quantum'`, `'homo_double_quantum_nupicycle'`, `'homo_double_quantum_nucycle'`.

### Key state/data transformations

- Lines 78: computes `rotor_period` using `rotor_period=1/mas_rate`.
- Lines 81: computes `r_element_dur` using `r_element_dur=n_rotor_periods*rotor_period/n_blocks_per_period`.
- Lines 90: computes `pulse_amp` using `pulse_amp=pi/r_element_dur`.
- Lines 93: computes `pulse_dur` using `pulse_dur=r_element_dur`.
- Lines 96: computes `phases` using `phases=zeros(n_blocks_per_period,1)`.
- Lines 98: computes `phases(q+1)` using `phases(q+1)=((-1)^q)*pi*phase_factor/n_blocks_per_period`.
- Lines 109: computes `pulse_dur(1)` using `pulse_dur(1)=(1/4)*r_element_dur`.
- Lines 110: computes `pulse_dur(2)` using `pulse_dur(2)=(3/4)*r_element_dur`.
- Lines 115: computes `phases(2*q+1)` using `phases(2*q+1)=((-1)^q)*pi*phase_factor/n_blocks_per_period`.
- Lines 116: computes `phases(2*q+2)` using `phases(2*q+2)=phases(2*q+1)+pi`.

### Local helper functions

- Line 161: `grumble()` — `function grumble(n_rotor_periods,n_blocks_per_period,phase_factor,`.
  - Representative operation: `n_cycle_repeats,mas_rate,element_type,supercycle_type)`.
  - Representative operation: `if (~isnumeric(n_rotor_periods))||(~isreal(n_rotor_periods))|| (~isscalar(n_rotor_periods))||(mod(n_rotor_periods,1)~=0)|| (n_rotor_periods<1)`.

## Parameters / inputs

- n_rotor_periods -"small n" symmetry number, gives number of
- rotor periods required in the R symmetry
- n_blocks_per_period -"capital n" symmetry number, gives the number of
- R elements contained within the R symmetry
- phase_factor -"nu" to calculate the alternating phase in the
- R sequence:
- 180*nu/N = 180*phase_factor/n_blocks_per_period
- n_cycle_repeats -number of times the full R sequence is applied
- mas_rate -rotor spinning rate, Hz
- element_type -R element needs to be an inversion pulse. Common
- R elements are:
- '180_pulse' : simple inversion pulse
- '90270_pulse' : composite inversion pulse
- supercycle_type -The R sequence can be repeated multiple
- times in combination with supercycles, for
- improved performance, removal of undesired
- higher order terms. If the he unmodified R
- is denoted [phase], this can either be in-
- verted, [-phase], or have an overall phase
- added to it, [phase]_addph.
- Common supercycles are:
- 'hetero_single_quantum'
- [phase]_0:[-phase]_0:[phase]_120:[-phase]_120:[phase]_240:[-phase]_240
- 'homo_double_quantum_nucycle'
- [phase]_0:[-phase]_0
- 'homo_double_quantum_nupicycle'
- [phase]_0:[-phase]_0:[-phase]_180:[phase]_180
- Output:
- phases -the sequence of pulse phases, radians
- pulse_amp -RF nutation frequency in rad/s, a scalar because
- R-sequences are phase-modulated
- pulse_dur -duration of the pulses in the sequence element,
- a vector with the length matching the number of
- pulses in the sequence element (seconds)

## Implementation structure

- R-sequences described in Malcolm Levitt's review:
- Nomenclature is based on the following notation RN_{n}^{\nu}. Syntax:
- [phases,pulse_amp,pulse_dur]=rsequence(n_rotor_periods,n_blocks_per_period,...
- phase_factor,n_cycle_repeats,mas_rate,...
- element_type,supercycle_type)
- n_rotor_periods - "small n" symmetry number, gives number of
- rotor periods required in the R symmetry
- n_blocks_per_period - "capital n" symmetry number, gives the number of
- R elements contained within the R symmetry
- phase_factor - "nu" to calculate the alternating phase in the
- R sequence:
- 180*nu/N = 180*phase_factor/n_blocks_per_period

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `phases()`, `pulse_dur()`, `isscalar()`, `ischar()`.
