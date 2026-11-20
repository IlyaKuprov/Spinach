# Self-consistent one-axis CRAZED

Use `imaging(spin_system,@crazed_ddf,parameters)` for the complex signal after two ideal X pulses and two instantaneous whole-pitch gradient kicks. The sample is a uniform periodic grid; scalar coupling between like isotopes is supported. Molecular coherence filtering with `crazed` answers a different question.

The sequence calibrates local Cartesian spin projections against the longitudinal projection of `rho_ref` and the specified equilibrium magnetisation density in A/m. Changing the receive phantom changes detection, not feedback. A zero magnetisation density disables the feedback while preserving the prepared state and receiver. Sampling begins immediately after the second kick. Refine both spatial resolution and LG4 substeps to establish convergence.
