# Spinach development

## Scope and source of truth

Read and follow the distribution-root `AGENTS.md` as mandatory instructions
for every task, as required by [the skill entrypoint](../SKILL.md#mandatory-repository-instructions).
Consult that file directly for repository policy; this reference does not
maintain an independent version of its instructions. Preserve unrelated
changes and compiled MEX files.
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

## Tests, examples, and documentation

For GPU `propagator` work, two-sparse-operand Taylor and squaring products
use custom low-level CUDA CSR arithmetic when `cuda_sparse_by_sparse_mex`
is available. Missing or unloadable platform binaries retain native GPU
multiplication; storage-layout, computation, and validation failures propagate.
Exercise both stages and `clean_up` on real and complex GPU arrays, preserving
its round-to-grid, density, and disable policies.

Select tests from `tests/README.md` and `tests/lib/test_manifest.m`; add a
focused regression when behaviour changes. Run affected examples when they
exercise a different regime from the unit tests. Capture assertions and
observable comparisons, not only screenshots or log tails. Follow
[execution and validation](validation.md) for runtime evidence.

Follow the distribution-root `AGENTS.md` for required skill and knowledge-base
updates, content preservation, and documentation policy. It is the sole
source of those instructions.

For a documentation-only change, validate local links, frontmatter, code
examples, and claims against the current source. Do not claim numerical
coverage for unchanged algorithms merely because the documentation checker
passes. Use [skill evaluation](skill-evaluation.md) for skill changes.

## Review and delivery

Before committing, re-read the checkout `AGENTS.md`, inspect the full diff,
and confirm that only intended files changed. Record the exact checks and
limitations in the pull request. Do not claim tests that were not run.
Address reviewer findings on the actual current head and revalidate affected
behaviour; a review of an earlier commit is not a review of later changes.
Follow the repository and user review/merge policy; opening a PR does not
authorise merging it.
