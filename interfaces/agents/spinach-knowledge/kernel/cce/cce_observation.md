# cce_observation — actual preparation and detection reference

## Native role and invocation

Internal preparation/detection compiler used by the registered experiments; its input system is an already constructed retained subsystem.

Public simulations enter through `clexpn(spin_system,pulse_sequence,parameters,assumptions)` with a registered experiment.

Internal/expert signature: `[rho0,coil,reference,meta]=cce_observation(cluster_system,parameters)`. Engine records are not public experiment parameter structures.

Syntax: [rho0,coil,reference,meta]=cce_observation(cluster_system,parameters).

Expert rho0 is a positive trace-one Hermitian observed-spin density matrix; coil is an observed detector or a cell vector. Otherwise preparation is transition or a structure with states/spins for cce_state; detection is transition or a structure with operators/spins for cce_operator. IDs are physical. A selected transition supplies equal-superposition preparation and alpha*beta' detection unless explicitly overridden. Without transition or preparation/detection, missing physical choices are errors, not guessed polarization.

preparation/detection descriptors and receiver_phase are experiment-service fields; the context must validate them at the experiment boundary rather than passing generated cce_aux back to the stage-1 resolver. receiver_phase is radians and rotates detected signal by exp(-1i*receiver_phase). normalized is logical and defaults false at this service level: experiments own their defaults. reference and meta.initial_signal are the actual complex initial detection; meta.normalization_divisor equals reference if normalized and ones otherwise. Neither rho0 nor coil is divided here. Vanishing normalization references are rejected. Metadata declares time-by-detector output and column-major vectorisation.

The bath density is the retained product state in physical-ID order. Generalised output follows subsystem.formalism. Projected output is the actual prepared alpha/beta coherence times bath density, vectorised, with a trace detector carrying detector phase; only one selected-transition detector is supported. Preparation descriptors retain Spinach's unnormalised deviation-density convention.

The named preparation `transverse` is available for exactly one observed spin one-half and constructs `(identity+2*Lx)/2`. It does not require a transition when detection is supplied separately, and is rejected for multiple observed spins or higher spin. Explicit state descriptors, transition preparation and expert densities are also supported.

The default transition reference is evaluated in the bilinear order `exp(-1i*receiver_phase)*(alpha'*rho_observed*beta)`. This preserves the numerical convention of the selected-transition signal. Explicit coils or detection descriptors instead use `trace(coil'*rho_observed)`. The two forms represent the same detection convention but need not have identical floating-point rounding.

Receiver cells may be row or column vectors. Linear cell order defines output-channel order; references and normalisation divisors are always rows. Each normalisation decision uses only that detector's own scale.

Caller-supplied parameters.sequence is rejected, including empty arrays. Use finite_drive, refocusing_train, or the appropriate named echo experiment; internal events are constructed from native physical declarations, never read from that removed field.
