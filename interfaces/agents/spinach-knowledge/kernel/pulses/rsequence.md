# kernel/pulses/rsequence.m

- Signature: `[phases,pulse_amp,pulse_dur]=rsequence(n_rotor_periods,n_blocks_per_period,phase_factor,n_cycle_repeats,mas_rate,element_type,supercycle_type)`

## Purpose

Generates the pulse phases, RF nutation frequency, and segment duration(s) for an R-sequence, using the `RN_n^nu` notation cited in Malcolm Levitt's review.

## Element construction

The rotor period is `1/mas_rate` seconds, with `mas_rate` in Hz. One R element lasts `r_element_dur=n_rotor_periods/(mas_rate*n_blocks_per_period)` seconds. For `180_pulse`, the amplitude is `pi/r_element_dur`, the duration is the full element duration, and the phase for block index `q` is `(-1)^q*pi*phase_factor/n_blocks_per_period` radians. This is the phase increment written in the cited notation as `180*nu/N = 180*phase_factor/n_blocks_per_period` degrees.

For `90270_pulse`, the amplitude is `2*pi/r_element_dur`; the two durations are `r_element_dur/4` and `3*r_element_dur/4`. Each block contributes a phase and a second phase offset by `pi`, with the first phase following the same alternating formula. The generated phase list is then transformed by the selected supercycle and the resulting list is repeated `n_cycle_repeats` times.

## Supercycle phase patterns

- `hetero_single_quantum` concatenates the phases and their negatives, then the corresponding copies shifted by `2*pi/3` and `4*pi/3`; its notation is `[phase]_0:[-phase]_0:[phase]_120:[-phase]_120:[phase]_240:[-phase]_240`.
- `homo_double_quantum_nucycle` concatenates `[phase]_0:[-phase]_0`.
- `homo_double_quantum_nupicycle` concatenates `[phase]_0:[-phase]_0:[-phase]_180:[phase]_180`, equivalently `[phase, -phase, -phase+pi, phase+pi]`.

These are the implemented phase concatenations; this description does not assert their performance outside the source's stated purpose.

## Outputs and reference

- `phases` — sequence phases in radians.
- `pulse_amp` — scalar RF nutation frequency in rad/s.
- `pulse_dur` — pulse duration(s) in seconds: one value for `180_pulse`, two for `90270_pulse`.

For background, see [Malcolm Levitt's review](https://doi.org/10.1002/9780470034590.emrstm0551), and the [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/rsequence.m) and [Spinach wiki page](https://spindynamics.org/wiki/index.php?title=rsequence.m).
