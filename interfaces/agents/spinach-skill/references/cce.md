# Cluster correlation expansions

## Choose the physical approximation

Use `cce` for central-spin/bath dynamics when a hierarchy of small bath clusters
is more tractable than the complete Hilbert space. Keep every central spin in
every cluster. Conditional `method='cce'` projects the register onto explicit
orthogonal alpha/beta states and returns branch coherence; choose `gcce` when
central-state mixing, register density, arbitrary observables, or generalized
controls must be propagated. A conditional projection is a modelling choice,
not an automatic substitute for generalized dynamics.

`expansion='partition'` selects pCCE: freeze bounded physical groups, solve
within-group interactions exactly, and expand in the number of groups. Use it
for strongly correlated local groups; check grouping and inter-group order
separately. `method='mecce'` and `'megcce'` add explicitly declared local
GKSL jumps to conditional/generalized dynamics. They do not derive rates from
a bath or justify a universal product approximation for population T1.

`method='noise'`, `quantity='autocorrelation'` uses additive connected-cluster
subtraction rather than multiplying coherences. Declare operator ordering and
connected subtraction before computing a spectrum.

## Minimal conditional calculation

This deliberately small example checks the API and exact full-bath limit; it
is not a converged material model. Use a clean MATLAB path containing Spinach,
not EasySpin. `create` converts its public interaction units; do not reconvert
its output.

```matlab
% Define a two-spin lab-frame model
sys.isotopes={'13C','13C'}; sys.magnet=0;
inter.zeeman.scalar={0,0};
inter.coupling.scalar=cell(2); inter.coupling.scalar{1,2}=25;
spin_system=create(sys,inter);

% Declare the register, bath and conditional observable
p.central_spins=1; p.bath_spins=2;
p.method='cce'; p.expansion='standard';
p.quantity='coherence'; p.order=1;
p.central.alpha=[1;0]; p.central.beta=[0;1];

% Declare the Hamiltonian, connectivity and initial bath state
p.assumptions='labframe'; p.orientation=[0 0 0];
p.interaction_cutoff=Inf;
p.graph=struct('mode','complete','score','coupling_norm');
p.bath_state=struct('mode','mixed_product'); p.meanfield='none';

% Supply total sequence times and numerical budgets
p.times=linspace(0,0.02,21); p.sequence=[];
p.numerics=struct('atol',1e-9,'rtol',1e-9,...
                  'max_memory_bytes',2^28);
r=cce(spin_system,p);

% Inspect validity before using values
assert(all(r.valid,'all'));
assert(r.exact_reference);
```

For generalized density dynamics, set `p.method='gcce'`,
`p.quantity='density'` and `p.central.rho` to a normalized register density
matrix. Scalar generalized coherence/observable requests also require explicit
`p.normalized`; `p.operators` supplies one matrix or a cell of matrices for
`quantity='observable'`. Normalization uses `central.rho`, not the first
requested time.

## Parameter structure: units and actual defaults

The driver does **not** fill a broad default parameter structure. Unless a
default is named below, supply the field explicitly and obey its source
`grumble` contract. The numbers in the example are selected values, not defaults.

| Field | Meaning and units | Default/requirement |
|---|---|---|
| `central_spins`, `bath_spins` | Disjoint physical spin IDs in the created system | Required; register tensor order follows central IDs |
| `method`, `expansion`, `quantity`, `order` | Method above; standard/partition; coherence/density/observable/autocorrelation; maximum cluster-group order | Required |
| `central.alpha`, `.beta` | Normalized orthogonal register vectors | Conditional branches and coherence projection |
| `central.rho` | Normalized Hermitian positive register density | Required for generalized dynamics; conditional branches may use the alpha/beta superposition |
| `assumptions`, `orientation` | Stock assumption name or `{name,retention}`; three Euler angles in rad | Required |
| `interaction_cutoff` | Physical Hamiltonian distance cutoff, Angstrom | Required, `Inf` allowed; distinct from graph and mean-field cutoffs |
| `hamiltonian_backend` | auto/direct/stock | auto; unsupported direct models are rejected |
| `graph.mode` | complete/distance/coupling/explicit | Required; radius in Angstrom, threshold in rad/s, or adjacency as appropriate |
| `graph.score` | coupling_norm/resonance/supplied | coupling_norm; resonance additionally requires mismatch matrix and positive mismatch_scale in rad/s |
| `graph.bath_radius` | Bath selection radius in Angstrom | Inf; omitted sites require convergence checks |
| `groups`, `clusters`, `cluster_cap` | Physical group partition; explicit group-ID clusters; per-order caps | Singleton groups; generated connected family; Inf cap |
| `bath_state.mode` | mixed_product/explicit_product/sampled/interlaced | Required; explicit `rho` cells follow bath IDs |
| `bath_state.basis`, `.probabilities`, `.samples`, `.seed` | Product-state sampling declarations | Sampled/interlaced modes require basis/probabilities (cells or `uniform`), sample count and seed |
| `meanfield` | none/static exterior product-state field | Required; same global exterior assignment for every cluster of a state |
| `meanfield_radius` | Exterior-field cutoff in Angstrom | Inf when omitted |
| `times`, `sequence` | Nonnegative sorted total sequence times in seconds; pulse events | Required; `sequence=[]` is FID |
| `numerics.atol`, `.rtol`, `.max_memory_bytes` | Observable error tolerances and memory budget in bytes | Required; budgets can invalidate outputs or refuse allocation |
| `numerics.zero_policy`, `reducer` | mask/error; powers/direct | mask; powers. No epsilon denominator repair |
| `numerics.precision` | Arithmetic precision | Only double accepted when specified |
| `parallel.mode`, `.batch_size`, `.max_tasks` | serial/clusters/states; tasks per batch; intentional stop boundary | serial, one task (parallel: pool size), Inf |
| `checkpoint.file`, `.resume` | MAT checkpoint path and explicit resume flag | No checkpoint when file absent/empty; resume false |
| `analysis.fit` | Enable explicit-window envelope fit | false; fit model/window/weights/confidence must be supplied when enabled |

For second-order conditional projection, enable `central.second_order.enabled`
and supply `gap_tolerance` in rad/s and `overlap_tolerance`; disabled otherwise.
For pCCE supply `partition.mode=user|coupling|spatial|adaptive`,
`max_group_size`, `max_group_dimension`, and mode-specific fields.
`partition.order_semantics` defaults to `max_groups`; `disjoint` and
`disjoint_inter` have restricted order meanings described in `cce_partition`.
Adaptive grouping requires pilot records/scores and production holdout ranges;
it is not an online regrouping of production samples.

Interlacing uses `interlaced.max_states`, `reference_state`, and optional
`condition_sites`. Numerator and denominator subsets must share the same
conditional cylinder set. The helper explicitly enumerates within the budget;
it does not invent an approximate sampler when the space is too large.

## Pulse and open-system conventions

Ideal events declare `type='ideal'`, physical `target` IDs, axis, angle in rad,
`placement='fraction'|'absolute'`, and `at`. Fractional placement is relative
to each **total sequence duration**, not an interpulse delay. Conditional CCE
accepts complete-register ideal pi swaps; arbitrary register or bath pulses
need generalized CCE. Named hahn/cpmg/xy4/xy8/udd schedules are supported by
`cce_controls`; follow its `name,count,axis,target` contract.

Finite events additionally use `duration` in seconds, `amplitude` and
`detuning` in rad/s. They are constant rotating-frame segments; supply separate
events for a piecewise envelope. Simultaneous finite drives, and finite or
entangling exterior pulses, are refused. Static exterior one-site ideal
rotations update the sampled field exactly.

For ME-CCE, `dissipation.jumps` records contain ordered physical `support`, an
`operator` in that tensor order and a nonnegative `rate` in s^-1. A correlated
jump cut across a cluster boundary is refused. Dynamics uses column-major
vectorization; the module does not silently reinterpret a Hamiltonian matrix
as a Liouvillian.

## Convergence and diagnostics

1. First reproduce the full-bath result of a small problem, independently when
   possible. The exact-reference bypass checks the same physical model, not
   the validity of a conditional projection or a chosen bath truncation.
2. Increase order with bath, graph, state sampling and pulses fixed. Repeat
   independently for physical bath size, graph threshold/radius, cluster caps,
   Hamiltonian cutoff and exterior-field radius. For pCCE also vary grouping.
3. Increase product-state samples with a fixed seed, then check a different
   seed. Average complete per-state reductions, not raw factors across global
   states. Inspect reported sample uncertainty.
4. Refine a nested time grid and the finite-pulse discretization independently.
   Check `valid`, `abs_error_bound`, invalid reasons and density trace,
   Hermiticity and minimum eigenvalue. A truncated density is not repaired by
   clipping, renormalization or a positivity projection.
5. Use `cce_converge` with explicit `convergence.ladders` records (`field`,
   cell-row `values`), `atol`, `rtol`, and `max_work`. Supported automated
   ladders are order, bath_spins, graph.threshold, graph.radius,
   bath_state.samples and times. Other checks use `convergence.checks`, a cell
   of complete parameter structures. Sampling ladders require an independent
   seed check; time ladders require exact nested sample points.

An unresolved multiplicative denominator is a convergence/conditioning result,
not permission to divide by epsilon. Full-bath raw values can remain valid at
zeros that make an intermediate product invalid. Preserve the complex phase
and validity mask when plotting or exporting. A plateau at one order is not
evidence that bath size, graph pruning or sampling is converged.

`cce_spectrum` requires explicit convention `exp(+i*omega*t)`, angular
frequency grid in rad/s, one/two-sidedness, full/positive time domain,
rectangle/trapezoid weights and none/Hann window. Positive-time extension
requires `stationary=true`. There is no window-power renormalization.
`cce_analyse` only fits the explicitly declared real/magnitude envelope,
interval, model, weights and confidence choice; it does not label an arbitrary
fit as T1.

Parallel modes require a pool already started by `create` with explicit
`sys.parallel`; never start a nested pool per cluster. Checkpoint resumption
requires matching physical identity and task IDs. CSV export has a JSON
sidecar that is the round-trip authority; retain it.

## Source and regression routes

- [CCE driver and kernel entries](../../spinach-knowledge/kernel/cce/cce.md)
- [Independent conditional reference](../../spinach-knowledge/tests/kernel/test_cce_exact.md)
- [Generalized density reference](../../spinach-knowledge/tests/kernel/test_gcce_density.md)
- [Open-system reference](../../spinach-knowledge/tests/kernel/test_mecce_exact.md)
- [Numerical conditioning](../../spinach-knowledge/tests/kernel/test_cce_conditioning.md)
- [Convergence tests](../../spinach-knowledge/tests/kernel/test_cce_convergence.md)
- [Import/export tests](../../spinach-knowledge/tests/kernel/test_cce_import_export.md)

Read the local `kernel/cce/*.m` headers and `grumble` checks before extending
these structures; the implementation is the authority for supported fields.
