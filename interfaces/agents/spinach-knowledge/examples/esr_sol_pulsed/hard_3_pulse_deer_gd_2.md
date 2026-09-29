# examples/esr_sol_pulsed/hard_3_pulse_deer_gd_2.m

- Signature: `hard_3_pulse_deer_gd_2()`
- Source: [`hard_3_pulse_deer_gd_2.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/hard_3_pulse_deer_gd_2.m)

## Pulse-selective Gd(III) DEER model

This example calculates three-pulse DEER for two Gd(III) `E8` (spin-7/2) centres using brute-force time propagation and powder averaging. The probe pulse excites an outermost zero-field-splitting transition, while the pump pulse excites the central transition; both are treated as ideal pulses in the source. Unlike the preceding central-transition case, it explicitly constructs transition-selective pulse operators.

The static field is `3.5` T and both isotropic g values are `2.002319`. The zero-field-splitting eigenvalues for each centre are `[1e9, 1e9, -2e9]` Hz; the first and second centres have Euler angles `[0, π/9, 0]` and `[0, 4π/9, 0]` radians, respectively. Their coordinates are `[0,0,0]` and `30*[sind(20),0,cosd(20)]` Å: a 30 Å separation at 20° in the x–z plane. The untruncated basis is `zeeman-hilb`.

## Transition operators and sequence

The shared [`deer_3p_hard_deer` helper](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_dipolar/deer_3p_hard_deer.m) applies a hard `π/2` probe pulse, evolves for the configured interval, applies a hard `π` pump pulse, refocuses the trajectory, applies a hard `π` probe pulse, then records final evolution through probe-channel detection. The callback is passed to `powder` in the `'deer'` context (not `'deer-zz'`), with the `rep_2ang_1600pts_sph` grid and detailed output. No finite pulse widths or independent pulse offsets are configured.

The DEER trace is sampled at 100 intervals of 20 ns (`stepsize=2e-8` s), for 101 points over 0–2 μs. The source sets a `5e10` Hz frequency sweep and 1024 nominal spectrum points, with the FFT evaluated at four times that length. The hard, pump, and probe pulse FIDs receive exponential apodisation parameter 6.

## Observable, plots, and limitation

The figure contains the frequency-swept spectrum, probe and pump excitation profiles, and `-imag(deer.deer_trace)` versus time in seconds; the frequency axes are offset frequency in Hz. It is a plotted result only; no output file is written by the script. The source warns that the Gd spin echo can be very sharp and hard to capture because experimental zero-field-splitting distributions are omitted. The source estimates calculation time in minutes.

## Source

[Spinach example source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/hard_3_pulse_deer_gd_2.m)
