# clexpn

[MATLAB source](../../../../../kernel/contexts/clexpn.m)

Syntax: `[answer,diagnostics]=clexpn(spin_system,pulse_sequence,parameters,assumptions)`

Native context over create output and a five-argument experiment. Physical observed IDs are in `parameters.observed`; algorithm controls are in `parameters.cce`. It compiles subsystem generators, schedules raw records, and returns the natural numeric experiment output separately from diagnostics.

Supports mixed, explicit, and deterministic sampled product states; `after_expansion` and `internal_then_expansion` estimators; explicit, coupling, and spatial partitions; serial or outer parallel scheduling and exact checkpoint resume. Interlaced product states use `cce_evaluate` on every subset and the complete bath, with common parent conditioning. Irreducible factors are consumed once.

Free evolution and Hahn echo support GKSL dissipation; the context selects Liouville representation for jumps unless a representation is explicit. Generalised free evolution supports elementwise observed density reconstruction, which need not remain positive at truncated order. Noise correlation retains first moments and centres after ensemble averaging. No legacy parameter translation occurs. Adaptive partition controls live under parameters.cce.partition; independent convergence uses cce_converge. Omitted order uses the included bath order.

Hahn echo also supports generalised `output='density'`, with column-vectorised observed density rows and the same elementwise reduction/full-bath override as free evolution. Projected coherence blocks do not expose a density.

Sampled assignments use configuration-indexed Threefry substreams, preserving deterministic physical state IDs independently of worker scheduling. Density-only free evolution needs an observed initial density but no receiver operator.
Transition construction happens once, in the compiled model, and is not repeated in options validation. Coherent free-evolution raw records preserve measured spectral residuals and absolute propagation-error estimates. Exact-reference validity and budget satisfaction are distinct.

The polarised-bath, two-observed-centre, and dissipative-bath examples use native hahn_echo/free_evolution calls and return [answer,diagnostics]. Sampled probability tables are rows; the 'uniform' basis selector names the standard Zeeman basis independently of the supplied probabilities. Reduced-density columns use column-major matrix ordering. Independent dense and exhaustive reference calculations are in test_cce_examples.

Native noise descriptors carry physical bath-spin support without building a full-bath matrix. Graph selection uses these supports for unconnected or sampled noise, preserving cross terms. The native context computes optional parameters.spectrum only after ensemble reduction and returns diagnostics.spectrum; the Fourier convention remains explicit.

Use clexpn(spin_system,@finite_drive,parameters,assumptions) for a rectangular observed-spin drive followed by free evolution. parameters.drive supplies spins, axis, amplitude, detuning, and duration; no public pulse-event array is required.

Bath-pumped Hahn evaluation uses `clexpn(spin_system,@bath_pumped_hahn,parameters,assumptions)`. The five-argument experiment optionally returns `[E,aux]`; coherent projected propagation supplies spectral residuals, Hamiltonian round-off bounds, and signal error estimates to the native reducer. The pump selects canonical axial Zeeman frequencies, with explicit Hz or rad/s units, a nonnegative halfwidth and explicit open/closed interval convention, a fraction in [0,1], logical allow_empty, and optional distinct physical bath sites or isotope-label cells. No-field product-density ensembles allow selected spins outside a retained cluster; static exterior mean-field pumping remains explicitly unsupported by bath_pumped_hahn; the Cartesian exterior-rotation support belongs to hahn_echo refocusing and the internal segment service.

## Native adaptive groups and ensemble contracts

parameters.cce.partition.mode=adaptive is executed through the experiment, using explicitly supplied pilot_configs and heldout_configs. The frozen production partition and convergence diagnostics are returned in diagnostics.partition. diagnostics.reduction.sampling_standard_error preserves raw full-bath Monte Carlo uncertainty under internal-first averaging; truncated internal-first and single-draw uncertainty are NaN.

For projected evolution, `parameters.cce.second_order` is included in graph construction. Enabled virtual-level corrections therefore retain observed-mediated bath pairs unless `parameters.cce.mediated` is false.

The experiment owns cce.propagator; graph and cluster-selection calls receive the same resolved parameters as the experiment.

Before propagation, the resolved representation validates the numerical preference: eig for Hilbert dynamics, stock for Liouville dynamics, with expm additionally available for generalised Hahn echo. Unknown or incompatible preferences are rejected rather than discarded.

The disjoint and disjoint_inter expansion families fix the retained physical-group order at one and two, respectively, capped by the number of bath spins, independently of a supplied generic order. This resolution precedes cluster enumeration and reduction: disjoint_inter must retain inter-group pairs even when the caller supplied order one. The resolved order is visible in diagnostics.parameters.cce.order.

A projected request is rejected by the grumbler of an experiment that does not support projected generators. Finite_drive and refocusing_train therefore cannot silently discard driven observed-spin population transfer.

The native input boundary rejects parameters.sequence, including empty arrays. Finite, fractional, overlapping, and partially observed controls belong to finite_drive; ideal trains belong to refocusing_train. No raw-event parameter adapter is provided.

The optional cce.bath_radius is forwarded to graph selection. It retains bath sites within the stated Angstrom distance of any observed spin, including the boundary; diagnostics.graph reports retained and omitted physical IDs. This bath-domain cutoff is distinct from the graph edge radius. A cluster containing every radius-selected site supplies the exact response for that selected model. Omitted sites remain in the physical configuration and can contribute to a static exterior mean field.

An incomplete partition is rejected with `clexpn:partition` before its candidate groups are installed or production tasks are scheduled. Pilot/held-out failure and exhausted probe, capacity, or iteration budgets therefore cannot return an ordinary numeric answer, even for a one-output call. Completed partitions preserve their production grouping and diagnostics.

For a single crystal, supply `parameters.orientation` (Euler angles in radians). Alternatively supply `parameters.grid` for orientation-resolved cluster generators. The same four-argument context call returns the weighted signal or the individual responses (`parameters.sum_up=false`). Cluster reduction before and after orientational averaging are distinct estimators.

`cce.backend='analytical_hahn'` selects the analytical secular Hahn evaluator. It requires a system created with `sys.tols.prox_cutoff=0`: point dipoles are computed from `inter.coordinates`, and user-supplied couplings are added: electron-nuclear tensors through their longitudinal component, and isotropic bath exchange through the flip-flop element. The system must be one electron followed by a located bath of identical spin-1/2 nuclei, without relaxation theories; anisotropic bath tensors are rejected. It supports a complete mixed-state bath of identical spin-½ nuclei, a single projected electron transition, ideal midpoint refocusing, and orders 0–2. It streams irreducible pair factors instead of allocating a graph or raw records; checkpointing and nonstandard cluster policies are not available on this path. The screening bound is exposed in `diagnostics.reduction.abs_error_bound`; powder quadrature combines completed products.
