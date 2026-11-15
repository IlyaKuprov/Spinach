# cce_observation — actual preparation and detection reference

Syntax: [rho0,coil,reference,meta]=cce_observation(cluster_system,parameters).

Expert rho0 is a positive trace-one Hermitian observed-spin density matrix; coil is an observed detector or a cell row. Otherwise preparation is transition or a structure with states/spins for cce_state; detection is transition or a structure with operators/spins for cce_operator. IDs are physical. A selected transition supplies equal-superposition preparation and alpha*beta' detection unless explicitly overridden. Without transition or preparation/detection, missing physical choices are errors, not guessed polarization.

preparation/detection descriptors and receiver_phase are experiment-service fields; the context must validate them at the experiment boundary rather than passing generated cce_aux back to the stage-1 resolver. receiver_phase is radians and rotates detected signal by exp(-1i*receiver_phase). normalized is logical and defaults false at this service level: experiments own their defaults. reference and meta.initial_signal are the actual complex initial detection; meta.normalization_divisor equals reference if normalized and ones otherwise. Neither rho0 nor coil is divided here. Vanishing normalization references are rejected. Metadata declares time-by-detector output and column-major vectorisation.

The bath density is the retained product state in physical-ID order. Generalised output follows subsystem.formalism. Projected output is the actual prepared alpha/beta coherence times bath density, vectorised, with a trace detector carrying detector phase; only one selected-transition detector is supported. Preparation descriptors retain Spinach's unnormalised deviation-density convention.
