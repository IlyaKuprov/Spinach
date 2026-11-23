# kernel/cce/cce_options.m

## Native role and invocation

Native parameter validator and default resolver called by clexpn before model or generator construction.

Public simulations enter through `clexpn(spin_system,pulse_sequence,parameters,assumptions)` with a five-argument experiment.

Internal/expert signature: `parameters=cce_options(spin_system,parameters,pulse_sequence,assumptions)`. Engine records are not public experiment parameter structures.

[MATLAB source](../../../../../kernel/cce/cce_options.m)

## Purpose

Validate native clexpn input ownership and resolve nonphysical defaults, with native experiment and algorithm ownership.

## Contract

- All four inputs are required. The canonical system comes from create; its tensors are already in rad/s. No tensor, time, pump frequency, or jump rate is converted by this function.
- Required parameters are observed (distinct physical spin-ID row) and times (nonnegative strictly increasing row in seconds). A transition (ordered spins and ms/levels/explicit alpha-beta selector) is required for projected evolution; generalised experiments with an explicit physical state need not supply one. spins denotes isotope channel names, not site IDs. The explicit fourth argument owns assumptions. Transition resolution occurs during compilation, not in this validator; the compiled model holds the resolved states and frame-dependent metadata.
- Experiment-owned public controls are preparation or the mutually exclusive physical rho0, coil, output, normalized, nloops, drive, refocusing, train, pump, bath_jumps, noise, and spectrum. Raw sequence events are removed from the supported public interface; sequence records consumed by lower-level propagation are experiment-compiled internal data. Detailed experiment validation is still required; a function handle alone does not declare supported generators or experiments.
- Approximation and ensemble controls live under cce. Tolerances live in `sys.tols.cce_atol`, `cce_rtol` and `cce_mem`; the pool is `sys.parallel`; the bath state is given by top-level `bath_state` (cell row of density matrices) or `bath_samples` with `bath_seed` and optional `bath_basis`, `bath_probs`; `cce.mediated` (logical, default true) keeps the cluster edges mediated by the observed spins. Unknown top-level and declared nested-container fields are rejected. Reserved sampling/run/diagnostics fields are not accepted.
- Defaults are generalised evolution; all non-observed physical spins as the bath; clusters joined by retained couplings and by observed-spin mediation; maximally mixed product density; no mean field; after_expansion averaging; powers reducer; auto backend; standard expansion, or partition if groups are supplied; retain_raw=false; second_order.enabled=false; zero_policy=mask. The default zero-angle orientation is reported.
- order remains empty until the context resolves it. Missing numerical tolerances and memory limits are left for the context to resolve. The experiment owns additive versus multiplicative algebra; powers/direct is the reducer implementation choice, not that algebra.
- Explicit approximation is `projected` or `generalized`. Expansion controls preserve standard/partition/disjoint/disjoint_inter; static mean fields may use deterministic product densities or sampled/interlaced states; meanfield_radius is rejected unless the mean field is static. Groups partition selected physical bath IDs. Sampling requires basis, probabilities, samples, and seed; explicit uniform declarations are accepted, otherwise cell rows follow bath-ID order. Densities are Hermitian, positive, and trace one; validation never normalises them.
- `cce.projection='conditional'` opts into the diagonal register approximation for a first-order, closed-system Hilbert Hahn echo with a level-index transition and ideal midpoint refocusing. The default invariant projection is unchanged.
- `periodic_segments` is a truncated, fixed-order projected `free_evolution` signal on a single cycle. It requires power reduction and an explicit static exterior product state without a mean-field radius.
- Bath jumps require explicit support, operator, and nonnegative rate in inverse seconds. Named L+/L-/Lx/Ly/Lz jumps have one physical-spin support; matrix jumps match the ordered support dimension. No jump rate/direction or physical population is invented.
- Pump controls specify center/halfwidth with explicit Hz or rad/s units, fraction, allow_empty, and optional open/closed boundary. Optional sites and isotope-species cell rows intersect the selected bath; absent filters mean that selected bath. No frequency conversion occurs here; the experiment converts the declared input units. Nonzero second-order projection requires projected evolution and explicit gap/overlap criteria.

The context and the grumblers of the experiments jointly enforce this native input contract.

Native `bath_jumps` retains the complete physical GKSL support contract: selected observed spins as well as selected bath spins are allowed. Generalised dynamics may transfer observed populations; projected dynamics still rejects jumps that mix the declared manifolds. IDs outside the selected observed-plus-bath set are rejected.
Transition eigensystem construction belongs to compilation, after native options validation. `cce_options` does not construct the observed Hamiltonian or return resolved transition metadata.

## Array-free fluctuation descriptors

Instead of noise.terms, supply a nonempty structure array noise.descriptors.
Every element has operators (Lx, Ly, Lz, or E), spins (a nonempty row of
unique physical bath IDs), and coefficient (a finite real scalar in the
fluctuation units). It denotes the coefficient times the sum of those
one-spin operators. The experiment constructs only retained-cluster matrices;
selected exterior spins contribute zero in that retained cluster. All physical
IDs are validated before restriction, including for the empty cluster.
Descriptors and explicit matrix terms are mutually exclusive. Unknown fields,
non-Hermitian labels, observed IDs, duplicate IDs, and non-real coefficients
are rejected. Connected ensemble and Fourier conventions remain explicit.

parameters.drive is owned by finite_drive (also admitted during internal cce_compile and cce_generators validation); it is rejected for other experiments. The experiment validates every descriptor field.

## Parallel task ceiling

There are no scheduling options: cluster tasks run on the pool that `create` starts from `sys.parallel`.

The refocusing descriptor belongs to hahn_echo and its direct compiler/generator validation. Its physical matrix, support, and placement are checked by the experiment. Passing it to free_evolution is rejected rather than silently ignored.

The optional pump.boundary selects an open or closed physical frequency interval; only those two values are accepted. The native closed default is unchanged, while explicit open preserves the original strict CCE window and excludes an exact-boundary resonance.

## Checkpoint request ownership

A supplied cce.checkpoint requires a nonempty character filename; resume, if supplied, is a scalar logical. Interlaced bath-state requests and independent convergence controls cannot share a checkpoint ledger and are rejected before physical compilation. Adaptive partition pilot and held-out evaluations remain separate from the production scheduler checkpoint.

drive belongs to finite_drive. A refocusing cell row may accompany finite_drive or define refocusing_train; the single scalar hahn_echo descriptor is unchanged. train belongs to refocusing_train. Detailed waveform, pulse, and train validation remains experiment-owned.

Caller-supplied sequence fields are rejected before compilation, including empty event arrays. Use finite_drive, refocusing_train, hahn_echo, or another experiment with its native physical controls. Native finite-drive descriptors remain valid during generator-side validation of simultaneous GKSL relaxation.

Declared defaults are resolved before a single local validation stage checks ownership, nested vocabulary, and physical values. The seed range remains the exact mt19937ar interval 0 through 2^32-1; no physical population or numerical resource allowance is inferred.
