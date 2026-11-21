# Cluster correlation expansions

## Native context and experiments

Use `clexpn` for observed-spin/bath dynamics represented by a hierarchy of retained bath clusters. Supply a system returned by `create`; do not build a full-bath basis first. Every retained subsystem contains every observed spin, followed by its selected bath spins. Physical spin IDs remain those of the supplied system.

```matlab
[answer,diagnostics]=clexpn(spin_system,@hahn_echo,parameters,assumptions);
```

`parameters.observed` is the ordered row of observed physical IDs; `parameters.times` is a nonempty, strictly increasing row of nonnegative total elapsed times in seconds. `assumptions` is an explicit stock assumption string, not an approximation inferred from isotope names. The context returns the numerical experiment answer first and validity, error, convergence and execution evidence second. Invalid values remain invalid; inspect `diagnostics.reduction.valid` and `diagnostics.reduction.reason` before fitting or plotting an inferred decay.

The registered five-argument experiments in `experiments/cce/` are:

| Handle | Experiment and native restrictions |
|---|---|
| `@free_evolution` | Free evolution; generalised or projected generators; signal, observable, and density outputs; registered dissipative support |
| `@hahn_echo` | Ideal midpoint pi pulse; total duration is twice the free interval; signal, observable, and density outputs; registered dissipative support |
| `@cpmg_echo` | Ideal transition-selective x-pulse CPMG train with positive integer `parameters.nloops`; signal or observable output |
| `@refocusing_train` | Explicit ideal refocusing descriptors or named XY8/UDD trains; generalised Hilbert dynamics; signal, observable, or density output |
| `@bath_pumped_hahn` | Hahn echo with a spectrally selected bath pump; `parameters.pump` declares the window and selected sites/species |
| `@cce_eseem` | Two-pulse ESEEM with ideal pi/2 preparation and a midpoint pi pulse |
| `@finite_drive` | Generalised finite-duration rotating-frame control; signal, observable, and density outputs; registered dissipative support |
| `@noise_correlation` | Generalised bath correlation with additive cluster cumulants; `output='correlation'`; bath-only Lindblad dissipation supported |

`cce_capabilities` checks handle identity and rejects unregistered callbacks before generator construction. A five-argument signature alone does not establish compatibility. Experiments receive `(cluster_system,parameters,H,R,K)`. The context owns retained subsystem construction, state sampling, cluster enumeration, and reduction; the experiment owns propagation, preparation/detection, and pulse timing.

## Physical approximation and input ownership

Set `parameters.cce.approximation` to `generalized` (default) or `projected`. Generalised evolution retains the observed-spin Hilbert space; projected evolution uses explicit orthogonal alpha/beta states of an observed transition. Use generalised evolution for observed-state mixing, finite drives, and density output. Projected coherent generators require the selected manifolds to be invariant within their declared tolerance; projected dissipation has its own conditional cross-density propagation, not a scalar attenuation substitute. A transition selector does not itself prepare a population or polarisation.

Use `cce_transition(spin_system,parameters,assumptions)` to resolve an ordered physical observed-spin selector into states and diagnostics. Explicit alpha/beta states are appropriate when a degeneracy prevents an unambiguous automatic choice. Preparation and detection belong at experiment level: `transition`, `rho0`, `coil`, `preparation`, `normalized`, and `output`. The lower-level `cce_observation` service additionally accepts `detection` and `receiver_phase`; these are not admitted by the public context option validator. For public calls, include the receiver phase in the explicit complex `coil`. An explicit `rho0` is the observed-spin density matrix, not a full-bath state. `cce_observation` tensors it with the retained bath state and compiles physical-ID operators. Normalised experiments divide once by the actual complex initial detector response; the observation service itself defaults to unnormalised output. State whether a requested amplitude is absolute or normalised.

| Native field | Meaning |
|---|---|
| `cce.bath_spins` | Ordered physical bath IDs, disjoint from observed IDs; default is the complement |
| `cce.order` | Maximum retained order; omission uses the included-bath order; use the explicit convergence service for refinement |
| `cce.expansion` | `standard`, `partition`, `disjoint` or `disjoint_inter`; partitioned forms require groups or an explicit partition policy |
| `cce.bath_radius` | Optional observed-to-bath distance cutoff in Angstrom, inclusive at the boundary; distinct from the graph edge radius |
| `cce.graph` | `complete`, distance, coupling, or explicit connectivity, with the corresponding radius/threshold/sparse logical adjacency |
| `cce.bath_state` | Mixed or sampled product bath states; supplied samples/probabilities/seeds define the sampling problem |
| `cce.meanfield`, `cce.meanfield_radius` | Exterior product-state field policy and optional distance cutoff |
| `cce.groups`, `cce.partition` | Physical bath grouping and bounded partition policy |
| `cce.backend` | `auto`, `direct` or `stock` retained-Hamiltonian construction within validated support |
| `cce.numerics` | Absolute/relative accuracy and memory budget, precision, and conditioning policy |
| `cce.parallel`, `cce.checkpoint` | Scheduling and checkpoint controls |
| `cce.convergence` | Independent refinement controls consumed by `cce_converge` |
| `cce.basis` | Retained `zeeman-hilb` or `zeeman-liouv` basis declaration |

Explicit `cce.backend='stock'` forces the independent stock construction; use the reported generator backend to verify a cross-check. Observable output may omit `coil` when the selected transition supplies detection. Isotope-selected preparation descriptors act only on observed spins, even when the bath contains the same isotope. Projected Hahn, bath-pumped Hahn, CPMG, and two-pulse ESEEM all use the fixed physical receiver. After odd observed-spin manifold interchanges, adjoint the bath-coherence block including its initial preparation coefficient; do not conjugate the receiver or the initial normalisation reference. Zero time has no inversion, and manifold-preserving pulses do not change the coherence ordering. Non-spectral Hahn echoes report propagation and numerical pulse-action estimates in detected units; a projected branch swap is exact. The public bath-pumped experiment remains closed-system only.

`cce_options` validates the complete nested vocabulary. Unknown fields are errors. Algorithm controls belong under `parameters.cce`; physical input tensors belong in `sys`/`inter`; pulse, drive, preparation, and detector choices belong to the experiment. The internal model and reducer records are not public parameter containers.

## Units and tensor construction

Spinach input frequencies are in Hz; canonical tensors after `create` are already in rad/s and must not receive a second `2*pi`. The native `finite_drive` descriptor uses amplitude and detuning in rad/s; experiment-compiled event records are internal, not an additional public input form. Coordinates and distance cutoffs are in Angstrom, elapsed times in seconds, Euler angles and pulse angles in radians, and jump rates in inverse seconds. A pump window explicitly declares Hz or rad/s units.

`cce_compile` retains canonical tensors, physical ordering, assumptions, and state metadata without a full-bath basis. `cce_subsystem` constructs the retained stock system; `cce_generators` constructs compatible `H,R,K`, including observation metadata and retained jumps. Keep the full tensor physics and explicit assumptions consistent when comparing a direct retained-Hilbert assembly with stock Spinach construction. Do not introduce a second bath Zeeman term or dipolar interaction already present in the canonical input.

## Cluster algebra and advanced state treatment

The family is downward closed: every proper subset needed by the recurrence is present. The empty-bath response is the observed-system baseline, not an invented unit signal. Multiplicative irreducible factors divide each raw response by its proper-subset contributions, with propagated uncertainty and conditioning checks. Small unresolved denominators are not replaced by epsilon. A cluster spanning all `diagnostics.graph.bath_ids` supplies the exact response of the radius-selected bath, without dividing by baseline zeros. Excluded sites remain in the bath configuration and can still contribute to the exterior mean field; this is not exact evolution of the unpruned bath.

Generalised observables and density matrices require entrywise handling with the declared observed dimension and detector. A non-informative density entry remains invalid rather than being replaced by a population from another cluster. Noise correlation uses additive inclusion-exclusion, never multiplicative coherence division. The requested reduction algebra comes from the registered experiment.

For sampled bath states, one global product assignment is shared by every cluster of a sample. Choose `cce.averaging=after_expansion` to reduce each assignment before averaging, or `internal_then_expansion` to average raw records, including the empty reference, before reduction. These estimators need not agree. Static exterior fields must use that same assignment. Interlaced sampling averages explicit neighbour configurations at the appropriate irreducible-cluster stage; it is an advanced engine operation, not a new public experiment handle.

Partitioned CCE counts complete groups, not individual sites; subset closure must not split a group. Hybrid families combine group clusters with cross-group monomer clusters under the selected policy. Partition pilot scores are construction evidence, not a certificate of physical convergence. Refine graph thresholds, bath extent, cluster/group order, and state samples independently.

## Pulses, dissipation and noise

Hahn/CPMG timing is expressed against each total elapsed time, not an interpulse delay. CPMG uses `nloops` and x-phase transition inversions. `refocusing_train` owns literal XY8 phases, UDD positions, or explicit ideal-pulse descriptors; it is not a CPMG alias. The `finite_drive` descriptor declares observed physical `spins`, `axis`, `amplitude`, `detuning`, and `duration`; it supports rectangular or piecewise-constant drives. Observations may occur before, during, or after a drive. Native ideal descriptors use `time` in seconds or `fraction` of total elapsed time. Raw event arrays are removed from the supported public contract; named experiments compile their own events. Bath-pump selection is based on the declared transition-frequency window, not merely an isotope label; a complete ideal observed-spin manifold swap is the restricted projected-control case.

Local dissipators are supplied as `parameters.bath_jumps`, with ordered physical `support`, a matching `operator` matrix and a nonnegative `rate`. Conventional dissipative propagation requires a retained Liouville basis; explicitly tagged projected branches use their conditional representation. Jump supports must lie in the selected observed-plus-bath system. They are included once in the retained generator, separately from stock relaxation; nonzero exchange is rejected by `cce_generators`. For projected dissipative evolution, actual alpha/beta jump blocks propagate the cross-density. A raising operator pumps towards the upper magnetic sublevel; neither its population preparation nor an inferred longitudinal relaxation time should be assumed from the word pumping.

A dissipative noise estimate includes source and trace arithmetic even at zero elapsed time, and both first moments for connected correlations. It is in squared fluctuation units. `noise_correlation` owns the noise operator/correlation convention and the context reduces additive cluster cumulants. When requested, `parameters.spectrum` invokes `cce_spectrum` on the completed ensemble result and returns the spectrum in `diagnostics.spectrum`. Specify angular-frequency grid, sign convention, sidedness, time domain, quadrature, and window. Positive-time stationary extension is an explicit physical assumption. A correlation spectrum or envelope fit is not automatically a T1 prediction.

## Numerical evidence and resource limits

`cce.numerics.atol`, `rtol` and `max_memory_bytes` define the propagation/reduction budget. Omitted accuracy or resource values are resolved by `clexpn` and recorded in `diagnostics.policy`; lower-level compilation requires resolved values. Graph closure, state storage, retained generators, density propagation, and reduction buffers are memory consumers, not just the final answer.

Inspect numerical error, denominator conditioning, density validity, and sampling error separately. A converged cluster order does not prove bath-radius, graph-cutoff, partition or sample-count convergence. For controlled comparisons, hold the bath realisation, physical spin ordering, transition phase, pulse schedule, and normalisation fixed; compare complete-bath raw responses with an independently constructed small-system calculation.

`cce_analyse` fits only an explicitly declared observable, envelope, model, interval and weighting. `cce_import` and `cce_export` are explicit interchange utilities. `cce_spectrum` may be called explicitly, or requested through `parameters.spectrum` for post-processing the completed context correlation; no per-cluster Fourier transform is performed. Call `[answer,diagnostics]=cce_converge(spin_system,parameters,pulse_sequence,assumptions)` for independent native refinement ladders.

## Orientation grids

For one orientation, set `parameters.orientation` to a radian Euler triple. For a quadrature, supply `parameters.grid` as a shipped grid name or an explicit structure containing `alphas`, `betas`, `gammas` and `weights`. Call `clexpn(spin_system,pulse_sequence,parameters,assumptions)` directly in both cases. Set `parameters.sum_up=false` to return a column cell array of unweighted responses in input order. Zero-weight nodes are skipped; `diagnostics.grid_indices` maps each response to its original grid entry. `parameters.cce.averaging` distinguishes reducing each orientation before averaging from averaging matching raw cluster records before reduction.

Each orientation constructs its own tensors, transition, graph, pump selection, and dynamics. cce.averaging='after_expansion' averages completed expansions and their bounds/masks; 'internal_then_expansion' averages raw records and the observed-only reference across orientations and spin configurations before reduction, requiring a common cluster family. Never average Hamiltonians. Orientations are serial; cluster/configuration scheduling may own the existing pool. Inside an outer worker the nested scheduling is serial. The ordinary five-argument powder callback remains separate; see [Contexts](contexts.md).

## Native adaptive and sampling controls

Use parameters.cce.partition for bounded adaptive partitions, with separate native pilot_configs and heldout_configs. clexpn supplies the registered experiment to each exact patch calculation and returns the frozen groups in diagnostics.policy.partition. parameters.cce.parallel.max_tasks caps batch concurrency, not the number of physical tasks. Internal-first sampled full-bath uncertainty is computed from individual unaveraged exact responses; truncated internal-first, connected-noise internal-first, and single-draw estimates are NaN. Interlacing rejects a reference density inconsistent with its declared local basis assignment using cce:spin_config and dimension^2*eps, independently of propagation tolerances.

## Selected-transition numerical conventions

Default projected free evolution uses the two conditional bath eigensystems, preserves the branch-product order and treats zero elapsed time exactly. Its residual-derived propagation bounds use the bath dimension. Default generalised Hahn detection evaluates `alpha'*rho_observed*beta`; the preparation reference uses that same bilinear ordering and the declared receiver phase. Explicit coils and detection descriptors retain their Hilbert–Schmidt traces and operator norms. Do not replace one arithmetic ordering with an algebraically equivalent one when bitwise regression identity is required.

Nonprojected free evolution takes the partial trace before applying observed-spin detectors. Liouville signal, observable, and density outputs retain the same propagation action and arithmetic-error estimate, with detector and normalisation scaling where appropriate. Select `cce.basis.formalism='zeeman-liouv'` explicitly when testing a zero-rate master-equation limit: an empty jump list alone does not select that representation.

Projected second-order graph controls are forwarded from `parameters.cce.second_order`; an observed-mediated pair is retained unless the graph is explicitly approximate. Sampling and projected-pair regression tests use native numeric answers and separate diagnostics, with unchanged registered IDs.

For Liouville free evolution, coherent generalised generators use sparse stock-step action; active dissipative generators remain dense. This representation choice follows the declared generator contract, not a numerical smallness threshold.

### Native jump and noise validation

Zero-rate bath_jumps records contribute nothing and are skipped before projected manifold and support checks. Active split supports require a common retained group; active exterior jumps cannot drive a frozen mean field. Projected jumps must preserve both selected manifolds and the loss-operator subspaces, with scale-relative leakage checks. Generalised Liouville evolution permits observed-spin jumps.

Noise inputs belong exclusively to noise_correlation, including direct evaluator requests. Explicit noise matrices use the relative Hermiticity criterion 100*eps*norm(B,'fro'), without an absolute unit-scale floor. Both noise moments and first moments remain available for ensemble connected subtraction.

Connected sampled noise retains M_s-mu_s*mean(nu)-mean(mu)*(nu_s-mean(nu)) before reduction. These centred contributions determine after-expansion sampling uncertainty; internal-first connected uncertainty remains unavailable (NaN). Complete-bath values and propagation bounds use the exact raw reference after centring. Free-evolution master dynamics uses stock exponential action with dense generator arithmetic, including projected coherence blocks, and reports its dimension/norm-based floating-point estimate. The registered rates and noise tests use native clexpn and service contracts throughout.

Capability versions distinguish propagation and retained-control-density schemas; regenerate checkpoints from preceding contracts rather than reusing their raw records.

Receiver matrices supplied as a cell vector retain linear channel order for either row or column storage; output references and uncertainty factors are rows. Generalised Hilbert Hahn calculations may explicitly select cce.numerics.propagator='expm' for direct segment evolution and its segment arithmetic estimate. The default spectral Hahn path remains available; propagation-only policy is not a graph or reduction constraint.

For native exact master-equation tests, choose cce.basis.formalism='zeeman-liouv' explicitly for zero jumps. Direct compiled models use their declared model.formalism; observed preparations belong to experiment parameters, not bath-configuration overrides. Hahn refocusing may be a scalar descriptor with physical spins, a Hermitian integrated operator, and fraction; it replaces, rather than adds to, the midpoint observed rotation. Master echoes use stock interval action and preserve the zero-time initial state. Coherent generalised master generators use sparse action even when non-diagonal; active dissipative generators remain dense. Projected selected-transition FIDs use their manifold trace with explicit initial preparations. The direct evaluator checks retained time storage before allocating, and generator jump-storage preflight precedes matrix construction.

For off-centre echo sampling, refocusing accepts either fraction or absolute time (seconds), never both. The pulse is included at a nonzero sample endpoint and omitted before that endpoint or at zero elapsed time. A projected refocusing matrix must preserve or interchange the selected manifolds; an arbitrary projected pulse that mixes other manifolds is unsupported. Native pulse_time and pulse_applied diagnostics expose inclusion without reconstructing an event schedule.

For native physical-group tests, disjoint fixes group order one and disjoint_inter fixes group order two; the latter includes cross-group pair propagation regardless of a generic input order. cce_partition takes resolved native options and a compiled model carrying the registered pulse_sequence. Adaptive graph selection does not reinterpret the experiment-owned propagator preference; explicit stock master-equation propagation reaches every exact patch evaluator unchanged. Native spin configurations carry bath densities and stable state identifiers, while observed preparation remains in experiment parameters.

For explicit physicality auditing, compile the same model with `cce_compile`, then call `cce_diagnostics(model,diagnostics.raw,diagnostics.reduction,diagnostics.parameters)`. Density output uses the registered column-vectorised row ordering; no legacy result wrapper or separate density copy is required. This audit reports negative populations, trace/Hermiticity residuals, coverage, and numerical-budget failures without repairing the result.

For a projected master-equation calculation with no Lindblad operators, explicitly select `parameters.cce.basis.formalism='zeeman-liouv'`. Its zero dissipator is not an instruction to switch to Hilbert-space propagation; the requested representation and numerical error contract remain in force.

Use pump.center, halfwidth, and units for bath_pumped_hahn spectral selection; set boundary=open to preserve the strict original CCE interval, or boundary=closed for the inclusive native default. allow_empty is explicit. The auxiliary selected_sites list is in physical IDs. For a single selective rotation, hahn_echo.refocusing accepts either an integrated Hermitian operator or a Cartesian axis and angle, plus fraction or time. Only Cartesian one-site exterior rotations update product-state static mean fields; arbitrary exterior operators remain unsupported. A finite cce.meanfield_radius requires observed/exterior coordinates at propagation, filters by the minimum observed-spin distance, and does not require included-site geometry. Projected native propagation verifies invariant selected manifolds before evaluating an experiment.

Native convergence uses an explicit registered experiment and assumptions in `[answer,diagnostics]=cce_converge(spin_system,parameters,pulse_sequence,assumptions)`. Its independent native-path ladders preserve the original baseline, work cap, held-out seed requirement, shared-time coverage, and numerical validity. Agreement between later refinements never certifies an inaccurate baseline.

Interlaced native diagnostics distinguish conditional raw responses from irreducible factors: inspect raw.value for the physical average and raw.factor for the quotient, and never reduce a factor twice. The registered evaluator is shared by every denominator. A full-bath average bypasses singular lower ratios; unavailable lower prefixes remain masked, and one truncated reference does not estimate a sampling error. Preparation rho0 is experiment-owned, not part of a bath configuration.

For checkpointed native experiments, place file and optional logical resume under parameters.cce.checkpoint. Resume may change cce.parallel scheduling and batch size, but not physical tensors, tolerances, sampled configurations, time mesh, or experiment. Do not pass checkpoint requests to direct experiment/interlacing evaluators, or combine them with interlaced sampling or independent convergence ladders. Adaptive partition pilot and held-out work run without the production checkpoint.

For native parallel work, choose the pool through sys.parallel before create. cce.parallel.max_tasks limits the concurrent batch, not the number of physical cluster/configuration records. The memory preflight includes unprojected observed-plus-bath construction, simultaneous Hilbert/Liouville workspaces, all detector columns, and the master-generator arrays retained by completed tasks. Increasing worker count can therefore require a larger budget even when every individual cluster fits; no silent task dropping or model truncation occurs.

A native model's source_hash represents physical/numerical inputs, not execution choices. Changing checkpoint/resume, parallel batching, raw-retention preference, or convergence orchestration does not change this fingerprint; changes to tensors, tolerances, bath states, or experiment parameters do. This also keeps retained generator provenance identical across a partial restart.

For a piecewise finite pulse use finite_drive with drive.spins, axis, row amplitude
(rad/s), matching positive duration (s), and scalar detuning (rad/s). Optional
time supplies a scalar start or matching explicit slice starts; otherwise slices
start consecutively at zero. Observations may occur during a finite pulse. Add ideal pulses
with a refocusing cell row: spins, integrated Hermitian operator or axis/angle,
and absolute time or fraction. Do not supply raw event arrays. Pure ideal trains
use refocusing_train with that cell row, or train.name=xy8/udd, count, axis, and
spins. XY8 phases remain literal; neither train is silently treated as CPMG.
These ideal/mixed controls require generalised Hilbert dynamics, retain the full
drift, and reject unsupported exterior dynamics. Single master echoes continue
to use hahn_echo. Native control tests keep experiment-owned preparation and
normalisation, independent dense references, and all physical assertions.

Retained services validate their own consumed records; `cce_options` owns public input validation. `cce_schedule` and `cce_interlace` require the explicit registered evaluator and native parameters, with no method-name fallback. Internal graph, projection, state, control, and scalar-reduction records are numerical-engine inputs, not alternate public context parameters. Density extraction and reference normalisation within the internal control engine are file-local and preserve column-major ordering and complex detector phase.

The legacy entry points `cce`, `cce_bath`, `cce_cluster`, `cce_master`, `cce_noise`, `cce_observe`, `cce_prepare`, `cce_sequence`, and `cce_validate` are removed, with no adapters or alternate parameter forms. Use `clexpn` and the registered experiments; construct baths with `spin_bath`. Lower-level work uses the native compilation, subsystem, generator, and observation services described above.

## Source-pinned physical fixtures

Use [sys,inter,fixture]=cce_fixture('nv_four_c13','6f7282b6') for the fixed four-13C NV input dataset. It supplies full physical tensors, source metadata, coordinates as metadata, and the original 41-point total-time row through 0.5 ms. Apply observed-Sz projection explicitly in the example; select the transition, bath state, order, numerical policy, and execution resources separately. Do not replace this fixed bath with an isotope-abundance sample or shorten its grid when reproducing the fixture.

## Removed input ownership

The old observed_spins spelling is removed; use observed. Old top-level order, graph, bath_state, partition, interlaced, numerical, scheduling, and reduction controls belong under cce, while assumptions is the fourth context argument. method and quantity wrappers are replaced by the experiment capability and output; hamiltonian_backend becomes cce.backend, deer becomes pump with explicit frequency units and boundary, and jumps becomes bath_jumps. Supplied clusters and cluster_cap overrides are removed: the context owns the family, and the numerical policy owns its memory budget. Reserved controls, sampling, run, diagnostics, initial_rho, bath_state.temperature, and numerics.full_bath are not supported public inputs.

Raw sequence/name/event-array forms are rejected at runtime, including empty arrays. Use the registered experiment and its drive, refocusing, train, nloops, or pump controls instead. All finite-drive physics checks use native descriptors: independent-axis overlapping drives, fractional starts, partial observations, coupled drift, and simultaneous GKSL relaxation. Internal events are compiled inside cce_acquire from validated native declarations; no raw parameter field or alternate context API remains.

Use cce_analyse with an explicit record containing times, value, and valid and separate analysis options; that small data record is not a returned legacy simulation wrapper. Use cce_export(filename,model,payload) with an explicitly packed answer/diagnostics payload. Both utilities remain explicit post-processing choices.

finite_drive accepts a nonempty structure row of independent drives. Each has observed spins, axis, amplitude/duration rows, detuning, and optional absolute time; a single finite pulse may instead specify a start fraction in [0,1]. Durations remain seconds, and overlapping contributions add while the full drift and supported dissipation remain active.

Optional transition-derived coils remain optional for refocusing_train and mixed ideal/finite controls. Convergence order ladders use the resolved baseline bath, so bath_spins may be omitted. spin_bath preserves existing electron magnet-frequency defaults when it introduces a Zeeman table, while respecting explicit zero or nonzero Zeeman inputs. JSON/CSV export of a complete answer/diagnostics payload preserves native registered experiment contracts, sparse arrays, character shapes, and missing strings; import rejects a changed capability version.

For native CCE orientation averaging, a checkpoint ledger is supported only for one positive-weight orientation; a larger grid is rejected before any checkpoint write. Repeated quadrature nodes do not create independent spin draws. Single sampled or truncated-interlaced draws retain `NaN` uncertainty, and internal-first complete-bath uncertainty is estimated after averaging each spin draw over orientations. Exact full-bath interlacing retains zero sampling error and shape-consistent per-order reasons. Diagnostic maxima preserve supplied unavailable residuals rather than hiding them.

A resumed calculation rechecks its current memory allowance but does not hash that resource allowance into the physical checkpoint identity. This includes automatically resolved memory budgets on hosts where available memory changes between calls.

Convergence reports preserve a partially unavailable sampling error: a `NaN` at any shared requested time makes the aggregate sampling uncertainty `NaN`, even when other times have finite estimates.

Tensor-free appended electron baths retain standalone Zeeman and dipolar-scaling defaults even when the existing system declares Zeeman or susceptibility data; explicit bath tensors retain their supplied meaning. Finite Liouville drives report accumulated floating-point segment estimates with detector/reference scaling. Incomplete adaptive partitions raise clexpn:partition before production rather than returning an uncertified numeric answer.

Native retained services validate caller-controlled inputs at their own entry boundary; cce_options resolves declared defaults before its local validation, and cce_reduce accepts both indexed records and native ensemble cells, including exact full-bath interlaced records.

The compiler, partition selector, convergence driver, diagnostics, and orientation driver own their consumed input checks. Hahn experiment branches expose their propagation directly, and the native subsystem/generator services are the active state and dissipative-generator paths.

CPMG and ESEEM retain a floating-point propagation estimate per time and receiver, including the actual normalisation reference. A valid full-bath signal can exceed the numerical budget; zero elapsed time incurs no propagation error.

For a paper-scale benchmark, select the material example by figure rather than assuming a common output shape. The diamond NV and SiC clock examples return order-resolved echoes; the VO(TPP) example returns register energies and fields as well as the proton echo and times. The Finland and OX063 examples compute the published powder-averaged pair-approximation curves, not small-cluster orientation demonstrations. The corresponding knowledge entries distinguish archived input realisations from independently generated baths and document numerical discrepancies.
