# cce_capabilities — registered experiment contracts

`contract=cce_capabilities(pulse_sequence)` is the native registry consulted by `clexpn` before retained generators are built. It compares actual function handles; a matching name or five-argument signature does not register a callback.

| Experiment | Projection | Dissipation | Reduction and outputs |
|---|---|---|---|
| `free_evolution` | Generalized or projected | Yes | Multiplicative signal/observable; generalized observed density |
| `hahn_echo` | Generalized or projected | Yes | Multiplicative signal/observable; generalized observed density |
| `finite_drive` | Generalised only; mixed ideal/finite controls require Hilbert space | Pure finite drives only | Multiplicative signal/observable; observed density |
| `refocusing_train` | Generalised Hilbert only | No | Multiplicative signal/observable; observed density |
| `cpmg_echo` | Generalized or projected | No | Multiplicative signal/observable |
| `bath_pumped_hahn` | Generalized or projected | No | Multiplicative signal/observable |
| `cce_eseem` | Generalized or projected | No | Multiplicative signal/observable |
| `noise_correlation` | Generalized only | Bath-only jumps | Additive correlation |

Registered representations are `zeeman-hilb` and `zeeman-liouv`, subject to each experiment and generator restriction. `cce_generators` requires retained Liouville space for native generalized jumps and separately constructs conditional projected dissipative generators.

The registry returns the experiment handle, evaluator, supported formalisms, projection/dissipation flags, reduction algebra, density-reduction flag, outputs and contract version. The evaluator constructs the retained subsystem and generators and calls the five-argument experiment. It produces internal raw records for `cce_schedule` and `cce_reduce`; these records are not the public first output of `clexpn`.

For multiplicative experiments the raw record carries signal, observed-density and detector data where applicable, absolute error bounds, validity and propagation diagnostics. Hahn records include the actual physical retained IDs and raw response evidence. Noise records retain value, moment, mean_time, mean_zero, valid, and abs_error_bound for additive ensemble reduction.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_capabilities.m)

The noise evaluator is bound to noise_correlation, even when called directly without noise parameters. A foreign schedule or an explicitly different output is rejected with clexpn:noise before subsystem construction. Its contract version is distinct from the general dynamics evaluator.

Free evolution, finite drive, and refocusing trains reject checkpoints predating the native-only control boundary, as well as the preceding master-propagation and controlled-density schemas. The public experiment signatures are unchanged.

Direct registered evaluation checks retained time-history memory before constructing a subsystem or generator, including density components and retained density pages. Hahn registration covers native stock master-equation pulse propagation and generalized selective refocusing.

Hahn registration covers absolute refocusing time and manifold-preserving projected pulse descriptors. Raw records retain the experiment-owned pulse_time and pulse_applied arrays, independently of the observed numeric response.

The registered native evaluator rejects projected propagation when the actual retained Hamiltonian leaks out of either selected observed manifold, unless the explicitly requested second-order treatment applies. The criterion is dimension*eps*norm(H,fro), matching the original propagation check; direct generator inspection still exposes the measured leakage. Hahn registration covers Cartesian exterior-field updates, and bath-pumped Hahn covers the strict frequency window.

Direct registered experiment callbacks reject cce.checkpoint before constructing a subsystem. The scheduler removes its owned checkpoint request before callback invocation; this applies equally to free evolution, finite drive, master-equation propagation, and fluctuation correlations.

Master-equation raw tasks retain the native generator evidence returned by cce_generators, including physical spin and Lindblad operators. The experiment capability versions advance for this raw-record schema, preventing a checkpoint without those records from being resumed as the current schema. Capability versions are implementation-owned schema identifiers, not a user-selected physical approximation.

refocusing_train is an explicitly registered generalized Hilbert experiment with density/observable/signal reduction. Finite-drive and free-evolution checkpoint versions distinguish retained controlled-density evidence. The registered handles, not name heuristics, select evaluation.

Observable output without an explicit coil has one selected-transition detector; memory preflight counts a cell receiver only when supplied. Checkpoint capability versions also distinguish native backend selection and physical projected echo detection, so earlier raw tasks are not reused under the corrected contracts.
