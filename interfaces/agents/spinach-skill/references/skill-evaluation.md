# Evaluating and maintaining the skill

## Structure and portability

Keep `SKILL.md` as the short task router. Put specialist scientific detail in
one-hop references and executable semantics in the source/knowledge base.
Use descriptive headings and contents lists for long references. Avoid
agent-runtime-specific tool names, personal hosts, absolute installation
paths, and copied orchestration policies. Do not add scripts unless they
remove a demonstrated repeated source of error.

The frontmatter name is `spinach-skill`, matching the distributable directory.
Earlier revisions used `spinach-simulations`; agents with a cached/manual
registration under that name should refresh or rename that registration.
The entrypoint describes both when to activate and when not to activate.
The compatibility field declares runtime prerequisites, not an automatic
installation or permission to run expensive jobs.

After changes, parse YAML; check relative Markdown links and their anchors;
confirm all cited repository files and APIs still exist. Check scientific
claims against implementations and representative examples. Preserve useful
existing coverage during refactoring, moving it rather than dropping it.

## Behavioural evaluation cases

Run these prompts in fresh agent contexts where possible. Record which
references and source files were consulted, the produced script or answer,
and which assertions were actually executed. Compare with the old skill
under the same inputs when measuring an improvement. The expected behaviours
below are acceptance criteria, not a claim that evaluation has already run.

| Prompt | Expected behaviour |
|---|---|
| Simulate a two-proton liquid-state spectrum with supplied shifts, J, field, and sweep | Find a close liquid example, use `create`/`basis`/`liquid`, preserve Hz/ppm conventions, check output and axis |
| Simulate this spectrum; no field or nucleus is provided | Ask for the missing physical choices; do not invent values |
| Adapt a liquid example with symmetry groups to a spatially encoded experiment | Read `imaging`; do not retain unsupported symmetry reduction |
| Add Lindblad relaxation to a `zeeman-liouv` model | Inspect the selected dissipator; do not falsely require `sphten-liouv` for all relaxation |
| Run a double-rotor Hilbert-space sequence | Inspect `doublerot` callback; recognise the Hamiltonian stack rather than assuming a Liouville generator |
| No `inter.relaxation` was supplied, but a bosonic mode decays | Check mode damping/dephasing before asserting that `R` must vanish |
| A Krylov run finished and produced a plausible spectrum | Check assertions, observable, convergence, and physical limits; do not equate completion with correctness |
| A derivation suggests a kernel sign is wrong; fix it now | Read the primary source and exercise stock/changed code in relevant and counter-assumption regimes before making a defect claim |
| MATLAB is unavailable; provide a validated simulation | Offer source-grounded unexecuted code and state the runtime blocker; do not claim validation |
| The skill was installed outside the checkout | Locate the root and knowledge base explicitly; resolve MATLAB functions to the intended checkout |
| Refactor an unrelated Python service | Do not activate Spinach workflow |
| Use only EasySpin for this EPR model | Respect the package choice; do not force Spinach or put both toolboxes on the path |

A case passes only if its essential behaviour is observed, not because the
answer uses the expected keywords. Where execution is available, test the
script as well as judging the explanation. Separate routing/structure checks,
agent behaviour, MATLAB smoke tests, and scientific-regime validation in the
report; none substitutes for the others.

## Design references

Use the [Agent Skills specification](https://agentskills.io/specification)
for the portable format and
[authoring guidance](https://platform.claude.com/docs/en/agents-and-tools/agent-skills/best-practices)
for progressive disclosure and evaluation. Recheck current requirements when
changing packaging; do not encode a particular agent runtime into this skill.
