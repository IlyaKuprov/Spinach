# clexpn

Syntax: `[answer,diagnostics]=clexpn(spin_system,pulse_sequence,parameters,assumptions)`

Native context over create output and a registered five-argument experiment. Physical observed IDs are in `parameters.observed`; algorithm controls are in `parameters.cce`. It compiles subsystem generators, schedules raw records and returns the natural numeric experiment output separately from diagnostics.

Supports mixed, explicit and deterministic sampled product states; `after_expansion` and `internal_then_expansion` estimators; explicit, coupling and spatial partitions; serial or outer parallel scheduling and exact checkpoint resume. Interlaced product states use the supplied evaluator on every subset and the complete bath, with common parent conditioning. Irreducible factors are consumed once.

Free evolution and Hahn echo support GKSL dissipation; the context selects Liouville representation for jumps unless a representation is explicit. Generalized free evolution supports elementwise observed density reconstruction, which need not remain positive at truncated order. Noise correlation retains first moments and centres after ensemble averaging. No legacy parameter translation occurs. Adaptive partition and convergence controllers remain unavailable; omitted order uses the included bath order.

Hahn echo also supports generalized `output='density'`, with column-vectorised observed density rows and the same elementwise reduction/full-bath override as free evolution. Projected coherence blocks do not expose a density.

Sampled assignments use configuration-indexed Threefry substreams, preserving deterministic physical state IDs independently of worker scheduling. Density-only free evolution needs an observed initial density but no receiver operator.
Resolved transition policy metadata is copied from the compiled model; transition construction is not repeated in options validation. Coherent free-evolution raw records preserve measured spectral residuals and absolute propagation-error estimates. Exact-reference validity and budget satisfaction are distinct.

The polarized-bath, two-observed-centre, and dissipative-bath examples use native hahn_echo/free_evolution calls and return [answer,diagnostics]. Sampled probability tables are rows; the 'uniform' basis selector names the standard Zeeman basis independently of the supplied probabilities. Reduced-density columns use column-major matrix ordering. Independent dense and exhaustive reference calculations are in test_cce_examples.
