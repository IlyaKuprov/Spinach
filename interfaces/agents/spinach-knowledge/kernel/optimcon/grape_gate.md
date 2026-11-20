# GRAPE gate fidelity

Use control.gate_targets instead of state pairs for terminal rectangle Hilbert-space control. The objective is the squared magnitude of the trace overlap divided by the squared Hilbert dimension. Use rho_drift to match one target to each drift. Cartesian and phase wrappers share the ordinary waveform pullbacks; gate Hessians and density trajectory options are not supported.

Gate targets reject enabled Bloch–Siegert corrections and nonzero dead time. Explicitly disabling the correction or setting zero dead time is equivalent to omitting these settings.

Fidelity-only calls stream the running gate product without retaining slice histories. Gradient calls retain the histories needed by the backward sweep.
