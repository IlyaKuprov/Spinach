# Cluster correlation expansions

## Native context and experiments

Use `clexpn` for observed-spin/bath dynamics represented by a hierarchy of retained bath clusters. Supply a system returned by `create`; do not build a full-bath basis first. Every retained subsystem contains every observed spin, followed by its selected bath spins. Physical spin IDs remain those of the supplied system.

```matlab
[answer,diagnostics]=clexpn(spin_system,@hahn_echo,parameters,assumptions);
```

`parameters.observed` is the ordered row of observed physical IDs; `parameters.times` is a nonempty, strictly increasing row of nonnegative total elapsed times in seconds. `assumptions` is an explicit stock assumption string, not an approximation inferred from isotope names. The context returns the numerical experiment answer first and validity, error, convergence and execution evidence second. Invalid values remain invalid; inspect `diagnostics.valid` and `diagnostics.invalid_reasons` before fitting or plotting an inferred decay.

The registered five-argument experiments in `experiments/cce/` are:

| Handle | Experiment and native restrictions |
|---|---|
| `@free_evolution` | Free evolution; generalized or projected generators; signal, observable and density outputs; registered dissipative support |
| `@hahn_echo` | Ideal midpoint pi pulse; total duration is twice the free interval; signal, observable and density outputs; registered dissipative support |
| `@cpmg_echo` | Ideal CPMG/XY4/XY8/UDD trains via `parameters.decoupling`; signal or observable output |
| `@bath_pumped_hahn` | Hahn echo with a spectrally selected bath pump; `parameters.pump` declares the window and selected sites/species |
| `@cce_eseem` | Two-pulse ESEEM with ideal pi/2 preparation and a midpoint pi pulse |
| `@finite_drive` | Generalized finite-duration rotating-frame control; signal, observable and density outputs; registered dissipative support |
| `@noise_correlation` | Generalized bath correlation with additive cluster cumulants; `output='correlation'`; not dissipative |

`cce_capabilities` checks handle identity and rejects unregistered callbacks before generator construction. A five-argument signature alone does not establish compatibility. Experiments receive `(cluster_system,parameters,H,R,K)`. The context owns retained subsystem construction, state sampling, cluster enumeration and reduction; the experiment owns propagation, preparation/detection and pulse timing.

## Physical approximation and input ownership

Set `parameters.cce.approximation` to `generalized` (default) or `projected`. Generalized evolution retains the observed-spin Hilbert space; projected evolution uses explicit orthogonal alpha/beta states of an observed transition. Use generalized evolution for observed-state mixing, finite drives and density output. Projected coherent generators require the selected manifolds to be invariant within their declared tolerance; projected dissipation has its own conditional cross-density propagation, not a scalar attenuation substitute. A transition selector does not itself prepare a population or polarization.

Use `cce_transition(spin_system,parameters,assumptions)` to resolve an ordered physical observed-spin selector into states and diagnostics. Explicit alpha/beta states are appropriate when a degeneracy prevents an unambiguous automatic choice. Preparation and detection belong at experiment level: `transition`, `rho0`, `coil`, `preparation`, `detection`, `receiver_phase`, `normalized` and `output`. An explicit `rho0` is the observed-spin density matrix, not a full-bath state. `cce_observation` tensors it with the retained bath state and compiles physical-ID operators. Normalized experiments divide once by the actual complex initial detector response; the observation service itself defaults to unnormalized output. State whether a requested amplitude is absolute or normalized.

| Native field | Meaning |
|---|---|
| `cce.bath_spins` | Ordered physical bath IDs, disjoint from observed IDs; default is the complement |
| `cce.order` | Maximum retained order; omission requests the context convergence policy, not an assertion of full-bath convergence |
| `cce.expansion` | `standard`, `pcce` or `hybrid`; partitioned forms require groups or an explicit partition policy |
| `cce.graph` | `complete`, distance, coupling or explicit connectivity, with the corresponding radius/threshold/sparse logical adjacency |
| `cce.bath_state` | Mixed or sampled product bath states; supplied samples/probabilities/seeds define the sampling problem |
| `cce.meanfield`, `cce.meanfield_radius` | Exterior product-state field policy and optional distance cutoff |
| `cce.groups`, `cce.partition` | Physical bath grouping and bounded partition policy |
| `cce.backend` | `auto`, `direct` or `stock` retained-Hamiltonian construction within validated support |
| `cce.numerics` | Absolute/relative accuracy and memory budget, precision and conditioning policy |
| `cce.parallel`, `cce.checkpoint` | Scheduling and checkpoint controls |
| `cce.convergence` | Context convergence controls; inspected through returned diagnostics |
| `cce.basis` | Retained `zeeman-hilb` or `zeeman-liouv` basis declaration |

`cce_options` validates the complete nested vocabulary. Unknown fields are errors. Algorithm controls belong under `parameters.cce`; physical input tensors belong in `sys`/`inter`; pulse, drive, preparation and detector choices belong to the experiment. The internal model and reducer records are not public parameter containers.

## Units and tensor construction

Spinach input frequencies are in Hz; canonical tensors after `create` are already in rad/s and must not receive a second `2*pi`. The native `finite_drive` descriptor uses amplitude and detuning in rad/s, while explicit pulse events declare their units. Coordinates and distance cutoffs are in Angstrom, elapsed times in seconds, Euler angles and pulse angles in radians, and jump rates in inverse seconds. A pump window explicitly declares Hz or rad/s units.

`cce_compile` retains canonical tensors, physical ordering, assumptions and state metadata without a full-bath basis. `cce_subsystem` constructs the retained stock system; `cce_generators` constructs compatible `H,R,K`, including observation metadata and retained jumps. Keep the full tensor physics and explicit assumptions consistent when comparing a direct retained-Hilbert assembly with stock Spinach construction. Do not introduce a second bath Zeeman term or dipolar interaction already present in the canonical input.

## Cluster algebra and advanced state treatment

The family is downward closed: every proper subset needed by the recurrence is present. The empty-bath response is the observed-system baseline, not an invented unit signal. Multiplicative irreducible factors divide each raw response by its proper-subset contributions, with propagated uncertainty and conditioning checks. Small unresolved denominators are not replaced by epsilon. A complete physical full-bath family uses the exact full-cluster response where the reducer contract permits it.

Generalized observables and density matrices require entrywise handling with the declared observed dimension and detector. A non-informative density entry remains invalid rather than being replaced by a population from another cluster. Noise correlation uses additive inclusion-exclusion, never multiplicative coherence division. The requested reduction algebra comes from the registered experiment.

For sampled bath states, one global product assignment is shared by every cluster of a sample. Perform cluster reduction for each assignment, then ensemble-average; averaging raw cluster signals before nonlinear reduction changes the approximation. Static exterior fields must use that same assignment. Interlaced sampling averages explicit neighbor configurations at the appropriate irreducible-cluster stage; it is an advanced engine operation, not a new public experiment handle.

Partitioned CCE counts complete groups, not individual sites; subset closure must not split a group. Hybrid families combine group clusters with cross-group monomer clusters under the selected policy. Partition pilot scores are construction evidence, not a certificate of physical convergence. Refine graph thresholds, bath extent, cluster/group order and state samples independently.

## Pulses, dissipation and noise

Hahn/CPMG timing is expressed against each total elapsed time, not an interpulse delay. CPMG/XY trains declare pulse count and phase cycle; UDD positions follow the declared sequence. The `finite_drive` descriptor declares observed physical `spins`, `axis`, `amplitude`, `detuning` and `duration`; it applies a rectangular drive followed by free evolution. Each requested time must be zero or at least the drive duration. Explicit pulse events separately declare placement and units. Competing drive and pulse specifications are rejected. Bath-pump selection is based on the declared transition-frequency window, not merely an isotope label; a complete ideal observed-spin manifold swap is the restricted projected-control case.

Local dissipators are supplied as `parameters.bath_jumps`, with ordered physical `support`, a matching `operator` matrix and a nonnegative `rate`. Conventional dissipative propagation requires a retained Liouville basis; explicitly tagged projected branches use their conditional representation. Jump supports must lie in the selected observed-plus-bath system. They are included once in the retained generator, separately from stock relaxation; nonzero exchange is rejected by `cce_generators`. For projected dissipative evolution, actual alpha/beta jump blocks propagate the cross-density. A raising operator pumps toward the upper magnetic sublevel; neither its population preparation nor an inferred longitudinal relaxation time should be assumed from the word pumping.

`noise_correlation` owns the noise operator/correlation convention and the context reduces additive cluster cumulants. When requested, `parameters.spectrum` invokes `cce_spectrum` on the completed ensemble result and returns the spectrum in `diagnostics.spectrum`. Specify angular-frequency grid, sign convention, sidedness, time domain, quadrature and window. Positive-time stationary extension is an explicit physical assumption. A correlation spectrum or envelope fit is not automatically a T1 prediction.

## Numerical evidence and resource limits

`cce.numerics.atol`, `rtol` and `max_memory_bytes` define the propagation/reduction budget. Omitted accuracy or resource values are resolved by `clexpn` and recorded in `diagnostics.policy`; lower-level compilation requires resolved values. Graph closure, state storage, retained generators, density propagation and reduction buffers are memory consumers, not just the final answer.

Inspect numerical error, denominator conditioning, density validity and sampling error separately. A converged cluster order does not prove bath-radius, graph-cutoff, partition or sample-count convergence. For controlled comparisons, hold the bath realization, physical spin ordering, transition phase, pulse schedule and normalization fixed; compare complete-bath raw responses with an independently constructed small-system calculation.

`cce_analyse` fits only an explicitly declared observable, envelope, model, interval and weighting. `cce_spectrum`, `cce_import` and `cce_export` are explicit analysis/interchange utilities, not implicit context stages. The standalone `cce_converge` helper has no native `clexpn` invocation contract in the current implementation; perform explicit native runs for independent refinement ladders instead.

## Orientation composition

`clexpn` accepts a fixed `parameters.orientation`. For orientation averaging, obtain a deferred task with `parameters.context.mode='defer'`, then pass that descriptor to `orientation_average` with the declared grid/rank and matching assumptions. Reduction occurs separately at each orientation before the orientation-weighted numerical answers are combined. Never average Hamiltonians or irreducible factors as a replacement for averaging completed experiment results. `powder` can accept the same registered deferred descriptor; its ordinary five-argument callback path is distinct. See [Contexts](contexts.md).

## Native adaptive and sampling controls

Use parameters.cce.partition for bounded adaptive partitions, with separate native pilot_configs and heldout_configs. clexpn supplies the registered experiment to each exact patch calculation and returns the frozen groups in diagnostics.policy.partition. parameters.cce.parallel.max_tasks caps batch concurrency, not the number of physical tasks. Internal-first sampled full-bath uncertainty is computed from individual unaveraged exact responses; truncated internal-first, connected-noise internal-first, and single-draw estimates are NaN. Interlacing rejects a reference density inconsistent with its declared local basis assignment using cce:spin_config and dimension^2*eps, independently of propagation tolerances.

## Selected-transition numerical conventions

Default projected free evolution uses the two conditional bath eigensystems, preserves the branch-product order and treats zero elapsed time exactly. Its residual-derived propagation bounds use the bath dimension. Default generalized Hahn detection evaluates `alpha'*rho_observed*beta`; the preparation reference uses that same bilinear ordering and the declared receiver phase. Explicit coils and detection descriptors retain their Hilbert–Schmidt traces and operator norms. Do not replace one arithmetic ordering with an algebraically equivalent one when bitwise regression identity is required.

Nonprojected free evolution takes the partial trace before applying observed-spin detectors. Liouville signal, observable and density outputs retain the same propagation action and arithmetic-error estimate, with detector and normalization scaling where appropriate. Select `cce.basis.formalism='zeeman-liouv'` explicitly when testing a zero-rate master-equation limit: an empty jump list alone does not select that representation.

Projected second-order graph controls are forwarded from `parameters.cce.second_order`; an observed-mediated pair is retained unless the graph is explicitly approximate. Sampling and projected-pair regression tests use native numeric answers and separate diagnostics, with unchanged registered IDs.

For Liouville free evolution, exactly diagonal generators use sparse stock-step action; non-diagonal generators preserve their established storage and arithmetic. This representation choice uses `isdiag`, not a numerical threshold.
