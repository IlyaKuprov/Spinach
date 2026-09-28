# kernel/pulses/rsequence.m

- Signature: `[phases,pulse_amp,pulse_dur]=rsequence(n_rotor_periods,n_blocks_per_period,phase_factor,n_cycle_repeats,mas_rate,element_type,supercycle_type)`

## Purpose

Generates pulse phases, RF amplitude, and pulse durations for an R-sequence. The source cites Malcolm Levitt's review and uses the notation `RN_{n}^{nu}`.

## Algorithm

The rotor period is `1/mas_rate`; one R element lasts `n_rotor_periods/(mas_rate*n_blocks_per_period)`. For `'180_pulse'`, the function sets `pulse_amp=pi/r_element_dur`, uses that full duration, and alternates phases according to `(-1)^q*pi*phase_factor/n_blocks_per_period`. For `'90270_pulse'`, it sets `pulse_amp=2*pi/r_element_dur`, splits the element into durations `r_element_dur/4` and `3*r_element_dur/4`, and gives the second pulse in each pair a phase offset of `pi`. It then applies the requested supercycle and repeats the resulting phase sequence `n_cycle_repeats` times.

## Parameters / inputs

- `n_rotor_periods` — “small n” symmetry number: rotor periods required in the R symmetry; positive integer.
- `n_blocks_per_period` — “capital N” symmetry number: R elements in the R symmetry; positive integer.
- `phase_factor` — “nu” phase factor; the alternating phase increment is `pi*phase_factor/n_blocks_per_period` radians (equivalently `180*phase_factor/n_blocks_per_period` degrees).
- `n_cycle_repeats` — positive integer number of times to apply the full sequence.
- `mas_rate` — rotor spinning rate, in Hz.
- `element_type` — `'180_pulse'` (simple inversion) or `'90270_pulse'` (composite inversion).
- `supercycle_type` — `'hetero_single_quantum'`, `'homo_double_quantum_nucycle'`, or `'homo_double_quantum_nupicycle'`. Their phase patterns are respectively `[phase]_0:[-phase]_0:[phase]_120:[-phase]_120:[phase]_240:[-phase]_240`, `[phase]_0:[-phase]_0`, and `[phase]_0:[-phase]_0:[-phase]_180:[phase]_180`.

## Outputs

- `phases` — pulse phase sequence, in radians.
- `pulse_amp` — scalar RF nutation frequency, in radians per second.
- `pulse_dur` — pulse duration(s) for the sequence element, in seconds; a scalar for `'180_pulse'`, two durations for `'90270_pulse'`.

For background, see [Malcolm Levitt's review](https://doi.org/10.1002/9780470034590.emrstm0551) and the [Spinach wiki page](https://spindynamics.org/wiki/index.php?title=rsequence.m).
