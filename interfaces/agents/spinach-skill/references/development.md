# Spinach development

## Scope and source of truth

Read and follow the distribution-root `AGENTS.md` as mandatory instructions
for every task, as required by [the skill entrypoint](../SKILL.md#mandatory-repository-instructions).
Consult that file directly for repository policy; this reference does not
maintain an independent version of its instructions. Compliance is a hard
requirement without exception: check the complete diff against `AGENTS.md`
rule by rule before delivering, and correct every departure. Preserve unrelated
changes and compiled MEX files. For source-style auditing, `exorcise('offline')`
skips network checks; `exorcise('online')` also checks documentation URLs
through native MATLAB HTTP requests.
Use a separate branch/worktree when the current checkout belongs to another
task. Do not submit scientific or personal data that were not authorised for
publication.

Resolve the named function in the intended checkout. Read its header, body,
input checks, callers, shipped examples, tests, and matching knowledge entry.
Document the observation before proposing a cause. An attractive derivation
is not proof that a composed numerical algorithm or published result is wrong.

Before alleging such a defect, locate and read the primary source, run the
relevant shipped examples/tests both stock and changed, and cover the actual
application regime including a case outside the simplifying assumption in
the derivation. Without that evidence, report a question, not a defect.

## Make the smallest coherent change

- Keep physical and numerical choices visible. Do not replace a problem by
  an easier one merely to obtain a passing test. Keep sparse/matrix-free
  representations and tensor ordering intact.

## Where information lives

Spinach keeps each kind of information in one place and moves it by one
route. Follow the ownership below when adding or changing code, and do not
open a second route for information that already has an owner. The call
conventions are in [contexts](contexts.md#the-context-call) and the input
fields in [inputs](inputs.md).

- **The spin system owns the physics.** Spin system parameters,
  interactions, and constructed models live in subfields of `spin_system`.
  The user states them in `sys` and `inter`, consumed by `create`, and in
  `bas`, consumed by `basis`. `create` fills `spin_system.comp`, `.inter`,
  `.chem`, `.rlx`, `.tols`, and `.sys`; `basis` fills `spin_system.bas`.
  Code after `create` reads the constructed fields (for example
  `spin_system.inter.magnet` and `spin_system.comp.isotopes`), not the
  original `sys` and `inter`. A function that changes the model takes
  `spin_system` and returns it, as `assume`, `dipolar`, and `kill_spin`
  do; make a temporary model on a local copy.
- **`parameters` owns the experiment.** Channels (`parameters.spins`),
  offsets, sweeps, point counts, pulse and gradient settings, grids, the
  initial and detection states (`rho0`, `coil`), and requests to the
  context (`parameters.needs`) are subfields of `parameters`. The context
  adds defaults and derived fields such as `spc_dim` and `spn_dim` to the
  copy it hands to the sequence. A function that derives settings its
  caller needs returns the updated structure, as `fieldsweep` does.
- **Experiments are functions called by contexts.** A pulse sequence is a
  function, in `experiments/` or local to an example, that a context in
  `kernel/contexts` calls as `sequence(spin_system,parameters,H,R,K)`;
  `imaging` and `meshflow` append `G,F`. The context builds the
  generators and does the averaging; the sequence performs the experiment
  and returns the answer. Keep every argument slot, writing `~` for an
  unused one. Do not write an experiment inline in an example or a driver
  when a context and a sequence can express it; scans over settings,
  processing, and plotting stay in the caller.
- **Operators and matrices are single-letter variables passed as
  arguments.** `H`, `R`, and `K` arrive from the context, the sequence
  assembles `L` from them, and `P` is a propagator. Builders return
  operators as outputs (`hamiltonian`, `relaxation`, `kinetics`), and
  propagation functions take them as inputs. Do not park a generator in
  `spin_system` or `parameters` to carry it between functions.

To scan a physical quantity, change `sys` or `inter` and rebuild with
`create` and `basis`. To scan an experiment setting, keep the spin system
and change the field on a copy of `parameters`.

Other owners:

- Numerical tolerances and physical constants are in `spin_system.tols`.
  Users override them through `sys.tols`; `tolerances` installs them.
  Reuse an existing entry where one fits.
- Feature switches are in `sys.enable` and `sys.disable`. Text output goes
  through `report`.
- An optimal control problem is stated in `control`, which `optimcon`
  validates and returns as `spin_system.control`. Pulse timing, power
  levels, drift and control operators, and target states belong there,
  not in `parameters`. `drifts` obtains drift generators from a context.
- A finite-volume mesh is assigned to `spin_system.mesh`.
- Do not move simulation data through global variables or the base
  workspace; the library does not. `persistent` is reserved for static
  reference data, as in `spin`. Caches (the pool `ValueStore`, scratch
  files under `sys.scratch`) hold derived data that can be recomputed,
  never inputs.

Established exceptions follow. Do not "correct" them, and do not add a new
one without a reason of the same kind.

- **Operators in `parameters`.** The callback signature is fixed, so an
  operator that the caller chooses for a particular sequence travels in a
  documented field: `pulse_op`, `mw_oper`, `ez_oper`, and `homodec_oper`
  are the ones `sim2liouv` converts. `imaging` takes `rlx_op` paired with
  `rlx_ph`, `meshflow` takes `H_op`, `R_op`, and `K_op` paired with
  `H_ph`, `R_ph`, and `K_ph`, and `gridfree` takes `add_terms`. A
  sequence that lists `'zeeman_op'` in `parameters.needs` receives
  `parameters.hzeeman`. The main generators stay positional.
- **Qualified operator names.** When several operators coexist, a
  qualifier keeps them apart: `Lx`, `Ly`, `Hmw`, `K_chem`. Collections
  are cell arrays with descriptive names, such as `control.operators`
  and `control.drifts`.
- **Numerical building blocks.** `evolution`, `krylov`, `step`,
  `propagator`, the pulse functions in `kernel/pulses`, operator algebra,
  grids, and unit conversions take their operands as explicit arguments.
  They are not sequences and need no `parameters` structure.
- **Direct-call experiments.** `fieldsweep`, `rapidscan`, `masdnp`,
  `solid_effect`, and `eqmag` take `(spin_system,parameters)` and build
  their own Hamiltonians; where a header says to call the function
  without a context, do so.
- **Experiment settings that resemble model data.** A working coupling
  that sets a delay (`parameters.J`) is not the coupling in `inter`. A
  swept or evaluation field is not `inter.magnet`. `parameters.tau_c`
  drives rotational diffusion in `gridfree`; `inter.tau_c` belongs to
  Redfield theory.
- **Data that is not a simulation.** Parser records stay in their own
  structures until an importer converts them to `sys` and `inter`.
  Fitting and data-processing functions take their data directly.

Before adding a field or an argument:

1. Search `create`, `basis`, the context, and the sequence header for a
   field that already carries the information.
2. Put a new physical input in `sys` or `inter` and parse it in `create`;
   put a new experiment input in a `parameters` field and check it in the
   `grumble` of the function that reads it.
3. Document the field in the header of the function that reads it.

## Tests, examples, and documentation

For GPU `propagator` work, two-sparse-operand Taylor and squaring products
use custom low-level CUDA CSR arithmetic when `cuda_sparse_by_sparse_mex`
is available. Missing or unloadable platform binaries retain native GPU
multiplication; storage-layout, computation, and validation failures propagate.
Exercise both stages and `clean_up` on real and complex GPU arrays, preserving
its round-to-grid, density, and disable policies. Gateways rebuilt from the current source and the updated Linux binary use
scalar placeholder expansion for nonempty output construction; the shipped
Windows binary requires rebuilding to obtain this optimisation. Memory tests
must include constructor temporaries, not just numerical-kernel workspace.

Select tests from `tests/README.md` and `list_tests()`; add a
focused regression when behaviour changes. Run affected examples when they
exercise a different regime from the unit tests. Capture assertions and
observable comparisons, not only screenshots or log tails. Follow
[execution and validation](validation.md) for runtime evidence.

When replacing a handwritten generator by a model declared through `sys`,
`inter`, and `bas`, confirm the assumption and frame conventions of the
current implementation, and test the full-space generator before trusting
a spectrum from a reduced sector: a projected Hamiltonian can match
exactly while omitted off-sector couplings still give unphysical
trajectories. Compare the commutators with the conserved quantities, and a
propagated trajectory, against an independent small reference. For a
relaxation migration, compare trace preservation, positivity, and a
representative trajectory against an independently assembled small
dissipative model. Keep export and testing machinery out of readable
examples.

Follow the distribution-root `AGENTS.md` for required skill and knowledge-base
updates, content preservation, and documentation policy. It is the sole
source of those instructions.

For a documentation-only change, validate local links, frontmatter, code
examples, and claims against the current source. Do not claim numerical
coverage for unchanged algorithms merely because the documentation checker
passes. Use [skill evaluation](skill-evaluation.md) for skill changes.

For `clebsch_gordan`, preserve the preliminary zero filters and the input ceiling.
Its compensated magnetic recurrence uses native doubles only. Run both unchanged
`examples/fundamentals/clebsch_gordan/cg_test_*` examples
at their shipped strict absolute tolerance; also check nonzero tails relatively
against independent high-precision references, since absolute tolerance alone
can accept a vanished tiny coefficient. The angular-momentum regression suite
contains representative reference, tail, exact-zero, and ceiling cases.

## Review and delivery

Before committing, re-read the checkout `AGENTS.md`, inspect the full diff,
and confirm that only intended files changed. Record the exact checks and
limitations in the pull request. Do not claim tests that were not run.
Address reviewer findings on the actual current head and revalidate affected
behaviour; a review of an earlier commit is not a review of later changes.
Follow the repository and user review/merge policy; opening a PR does not
authorise merging it.

## Cluster-expansion contract checks

Read `cce_capabilities` for supported experiment handles and representations, `cce_options` for public field ownership, and `list_tests()` for registered test IDs. Retained internal engine records are not alternative public parameter forms. Preserve physical IDs, complex receiver convention, total elapsed time, both ensemble/reduction orders, and numerical validity/error channels when changing this boundary. Test a new public-field rejection at `clexpn`, not only at a low-level service. Use the existing source-matched regression notes rather than creating a knowledge copy of the test registry.
