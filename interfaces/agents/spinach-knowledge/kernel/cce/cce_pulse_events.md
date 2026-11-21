# cce_pulse_events

[MATLAB source](../../../../../kernel/cce/cce_pulse_events.m)

Internal syntax: `events=cce_pulse_events(parameters)`

Compiles finite_drive and refocusing_train physical declarations into an internal
ordered event array. Drive units are rad/s and seconds; row amplitudes and
matching durations describe waveform slices, optionally with explicit starts. A nonempty drive structure row permits independently phased and detuned pulses; a single pulse can use fraction instead of time, with duration still in seconds.
Ideal declarations use integrated generators or Cartesian angle/axis pairs.
XY8 phases and Uhrig fractions are generated here, never inferred from an
experiment's function-name string. Raw public sequence input is refused even when empty. This
service compiles declarations but does not construct a spin Hamiltonian.
