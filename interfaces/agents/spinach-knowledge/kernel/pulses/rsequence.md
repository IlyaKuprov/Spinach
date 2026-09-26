# kernel/pulses/rsequence.m

- Signature: `[phases,pulse_amp,pulse_dur]=rsequence(n_rotor_periods,n_blocks_per_period,...`

## Purpose

R-sequences described in Malcolm Levitt's review: Nomenclature is based on the following notation RN_{n}^{\nu}. Syntax: [phases,pulse_amp,pulse_dur]=rsequence(n_rotor_periods,n_blocks_per_period,... phase_factor,n_cycle_repeats,mas_rate,... element_type,supercycle_type)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

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
