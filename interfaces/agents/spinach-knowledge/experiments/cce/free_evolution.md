# free_evolution

Syntax: [answer,aux]=free_evolution(spin_system,parameters,H,R,K).

Uncontrolled retained-cluster evolution with explicit transition, preparation, or observed density rho0, and transition detection or an expert coil. times contains elapsed seconds. Normalised signals use the actual complex initial detector response; density output is time-by-column-vectorised observed density and is never normalised as a signal.

This experiment rejects sequence, drive, refocusing, and train fields, including empty raw event arrays. Use finite_drive for finite controls, refocusing_train for ideal trains, or hahn_echo for a Hahn experiment. Shared preparation, propagation, detection, and numerical diagnostics are implemented by cce_acquire; that service accepts native declarations only, not a legacy event-field interface.

Conventional generators obey L=H+1i*R+1i*K. Hilbert evolution reuses one eigensystem. Generalised Liouville evolution uses the stock step action or the explicit stock propagator choice, with its established absolute floating-point estimate. Projected evolution carries a coherence block, not a physical density. Selected-transition FIDs retain branch-product ordering, the evaluated bath trace at zero time, and bath-dimension residual bounds. Detection follows the actual complex receiver and normalisation reference.

The optional aux carries absolute error bounds and applicable eigensystem, orthogonality, Hamiltonian, and backend evidence. The cce_acquire knowledge entry retains the detailed numerical conventions for all shared branches. At truncated cluster order, elementwise density reduction need not preserve positivity.

Elapsed times and the absence of control descriptors are checked locally before acquisition.
