# examples/optimal_control/pattern_pulse_2.m

- Signature: `pattern_pulse_2()`
- Source: [`examples/optimal_control/pattern_pulse_2.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/pattern_pulse_2.m)

## Purpose

Design a phase-modulated, transmitter-offset-selective excitation for a single `13C` spin. The source cites the Glaser-group paper [10.1016/j.jmr.2004.12.005](https://doi.org/10.1016/j.jmr.2004.12.005) and describes the calculation time as minutes.

## Spin model and target pattern

The script sets `sys.magnet=28.18`, places the carbon at zero chemical shift, and uses the full spherical-tensor Liouville-space basis (`sphten-liouv`, approximation `none`). It normalises the `Sx` and `Sz` states and uses them as the desired outcomes across 128 transmitter offsets from −4000 to 4000 Hz. The target arrays assign `Sz` at indices 1–20, 54–73, and 109–128, and `Sx` at the remaining points; these bands encode the requested offset pattern, rather than a reported achieved profile.

## Pulse design

There is one zero drift and two control operators, `Lx` and `Ly`, for `13C`; `Lz` supplies the offset operator. The initial state is `Sz` for every offset, paired with its corresponding target state. The phase profile has 300 slices of 2×10⁻⁵ s each (6 ms total), a constant amplitude profile, and a configured power level of `2*pi*2000`. Starting from a constant phase of `pi/4`, the script configures L-BFGS with a 200-iteration termination limit and calls `fmaxnewton(spin_system,@grape_phase,guess)`. The source also requests phase/control, robustness, and spectrogram plots.

## Evaluation in the example

After optimisation, the script converts the phase profile to Cartesian controls and propagates the initial state separately at each offset with `shaped_pulse_xy` and `expv-pwc`. It uses `parfor` to evaluate offsets, projects each final state onto normalised `Sx` and `Sz`, then plots those values against the target pattern. This describes the scripted diagnostic, not a successful optimisation run or a measured result.
