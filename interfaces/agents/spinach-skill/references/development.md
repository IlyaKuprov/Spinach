# Spinach development

## Scope and source of truth

Read the checkout `AGENTS.md` and applicable local instructions before edits.
This skill routes the work; it does not replace repository coding style or
maintainer instructions. Preserve unrelated changes and compiled MEX files.
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

- Reuse an existing kernel operation or experiment where possible. Do not
  introduce a dispatcher, cache, wrapper framework, or new user API for a
  local calculation. Preserve established defaults and error behaviour unless
  their change is part of the request.
- Follow `AGENTS.md` for MATLAB typography, comments, naming, validation, and
  example structure. Library input checks and example scripts have different
  rules; do not add library boilerplate to examples by habit.
- Keep physical and numerical choices visible. Do not replace a problem by
  an easier one merely to obtain a passing test. Keep sparse/matrix-free
  representations and tensor ordering intact.
- Run `checkcode` on changed MATLAB files and review the diff against the
  house style. A clean analyser report does not prove scientific correctness
  or compliance with every formatting rule.

## Tests, examples, and documentation

Select tests from `tests/README.md` and `tests/lib/test_manifest.m`; add a
focused regression when behaviour changes. Run affected examples when they
exercise a different regime from the unit tests. Capture assertions and
observable comparisons, not only screenshots or log tails. Follow
[execution and validation](validation.md) for runtime evidence.

After a code change, update the relevant skill reference and mirrored
`interfaces/agents/spinach-knowledge` entry when behaviour or operational
guidance changes. Preserve existing valid explanations; replace stale
claims explicitly rather than appending contradictory notes. Knowledge entries
are concise function explanations, not copied source, process logs, or line
number inventories. The repository intentionally has no knowledge entry for
`tests/lib/test_manifest.m` and no section-wide knowledge indexes.

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
