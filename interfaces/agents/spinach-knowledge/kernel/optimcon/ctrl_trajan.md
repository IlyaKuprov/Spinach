# kernel/optimcon/ctrl_trajan.m

Source: [kernel/optimcon/ctrl_trajan.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/ctrl_trajan.m)
Wiki: [Spinach documentation: ctrl_trajan.m](https://spindynamics.org/wiki/index.php?title=ctrl_trajan.m)

## Purpose

An internal diagnostic plotting routine for optimal-control runs. Plot selection is supplied through <code>spin_system.control.plotting</code>, typically via <code>optimcon.m</code>; this function produces plots and does not return an objective, constraint value, gradient, or adjoint. It accepts trajectory structures from GRAPE, one per ensemble member, and can display waveform diagnostics, trajectory analyses, and fidelity robustness information.

If the plotting selection is empty, the function returns before validating the other inputs. It also removes <code>trajectory</code> from the requested plot list and returns if no plots remain.

## Call and data

<code>ctrl_trajan(spin_system,waveform,traj_data,fidelities)</code>

- <code>waveform</code> is arranged in successive X/Y row pairs, with columns representing time slices. Phase, amplitude, spectrogram, and instantaneous-frequency plots require an even number of rows.
- <code>traj_data</code> is a cell array of trajectory data; trajectory diagnostics trace over spatial degrees of freedom and are not available for every Zeeman formalism.
- <code>fidelities</code> supplies the fidelity data used by the robustness histogram.

After those early returns, the routine requires <code>waveform</code> and <code>fidelities</code> to be real numeric arrays and <code>traj_data</code> to be a cell array. It does not locally validate finiteness or cross-check waveform columns against slice durations. There are no default input values.

The time axis is either slice number or cumulative <code>pulse_dt</code> in seconds. Rectangle controls are drawn with stairs and the final waveform column is appended; trapezium controls are drawn with linear plots. Plotted amplitudes and bounds use <code>mean(spin_system.control.pwr_levels)</code> and are divided by <code>2*pi</code> for Hz labels.

Spectrogram and instantaneous-frequency plots use only the initial run of exactly equal <code>pulse_dt</code> values, and require at least five such slices. The spectrogram is formed from the complex control <code>X-iY</code>; instantaneous frequency is evaluated from <code>X-iY</code> using the first slice duration. The fidelity-robustness display is a probability-density-normalised histogram.

The plotting routine is diagnostic: its displayed amplitudes, bounds, trajectory summaries and fidelity histogram are not a complete specification of optimisation constraints or the objective's gradient.

Trajectory panels trace out spatial degrees of freedom using the full compiled spin dimension `bas.offsets(end)`, then analyse the direct-sum trajectory with `trajan`.
