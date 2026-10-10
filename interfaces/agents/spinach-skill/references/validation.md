# Execution and validation

## Before launching

- Record the intended checkout and revision, MATLAB release, input data, and
  requested observable. Resolve `which create -all` and `which basis -all`
  after setup; a successful call from a different checkout is not validation.
- Run `fix_path` from the Spinach root. Its default mode resets the MATLAB
  path; use `fix_path('add')` only when retaining other required paths is
  intentional and conflicts have been checked. Never keep Spinach and
  EasySpin active on the same MATLAB path. Add the selected example folder
  explicitly: `fix_path` does not add the entire example tree.
- Preserve shipped and locally built MEX binaries. Do not delete them to
  clean a worktree or replace them with another platform build. Check
  `mexext` and function resolution when a native dependency matters.
- Set `sys.parallel={'processes',N}` before `create`, with `N` chosen from
  the allocated resources. `create` owns pool setup and can
  replace a foreign pool; do not pre-open one or assume a requested worker
  count was accepted. Inspect its report. Follow the actual execution host
  allocation, not a hard-coded laptop or cluster worker count.
- Inspect large arrays, orientation grids, time steps, ensembles, and output
  size before production. Start with a physically meaningful reduced case.
  Keep diagnostic plots out of tight objective/gradient loops.

## Run and collect evidence

Use MATLAB, not Octave. For unattended execution prefer `matlab -batch` with
a bounded runtime and a captured log. Confirm that MATLAB actually started;
collect the eventual exit status, error/warning text, and expected outputs.
A process identifier or scheduler submission is not evidence of completion.

Place a distinctive success marker only after scientific assertions and
required file writes. Check the marker, expected artefacts, and process exit
together. A nonzero exit is not a clean run: inspect whether computation
failed or a later shutdown failed. Independently verified saved results may
remain useful after a shutdown failure, but report that distinction and keep
the failing log; never silently waive the exit status or catch an error
merely to print success. Record which checks ran; static inspection cannot
substitute for execution.

For a reproducible script, retain its physical parameters, data sources,
initial conditions, numerical tolerances, random seed when relevant, and
processing conventions. Save numerical arrays as well as plotted images
when the deliverable will be compared or reused.

## Choose checks that can falsify the result

Use the checks relevant to the observable and model; this is not a demand
to run every check on every job. Define acceptance tolerances from the
requested accuracy, numerical precision, or experimental resolution.

| Calculation | Useful independent check |
|---|---|
| Closed-system propagation | Norm/trace conservation in the applicable representation; a single-spin precession limit |
| Dissipative propagation | Expected decay and stationary state; do not demand norm conservation |
| Coupled-spin spectrum | Zero-coupling limit or a small full-basis calculation |
| Restricted basis | Increase correlation order/connectivity range and compare the requested observable |
| Powder or rotor calculation | Independently refine orientation grid, rotor ranks, and time discretisation as applicable |
| Spectral processing | Check frequency sign, carrier, sweep, FFT dimension, and the acquisition observable |
| Optimal control | Independently propagate the final waveform over the stated ensemble; objective improvement alone is insufficient |
| Stochastic calculation | Repeat with controlled seeds and estimate sampling uncertainty |
| Matrix-free calculation | Compare actions against a small explicit representation only where materialisation is supported |

Change one numerical approximation at a time. Do not compensate a wrong
axis sign by conjugating data until the plot looks familiar. Zero filling
interpolates a sampled spectrum; it does not establish physical resolution.
Distinguish model error, basis truncation, sampling error, and processing.

## Tests and examples serve different purposes

From the Spinach root, after path setup:

```matlab
addpath('tests');
list_tests();
result=run_test('kernel/pauli_spin_half_algebra');
assert(strcmp(result.status,'PASS'));
```

Read `tests/README.md` and `tests/run_test.m` in the checkout for the current
interface. Select relevant registered tests using `list_tests` or inspect
`tests/lib/test_manifest.m`; do not assume the registry has a knowledge-base
entry. Broader selection uses `run_tests`, for example
`run_tests('pattern','relaxation')`. The runner throws on failed tests and
returns `PASS` statuses on success; retain the selected test IDs and results.

Examples demonstrate intended usage and may be expensive or produce only
figures. Tests provide explicit assertions but do not necessarily cover the
user regime. For an algorithm change, retain stock and changed results and
include a case that breaks the simplifying assumption used in the argument.

## When execution is blocked

Say exactly what was checked and what could not run: missing MATLAB/license,
unsupported binary, absent input data, or unavailable allocated resources.
Deliver a source-grounded script when useful, labelled unexecuted. Do not
invent measured accuracy, expected output files, or a successful test result.

For concentration-weighted spherical-tensor states, `test_cwdm_states` checks
exact block-wise weighting and independent thermal references, while
`test_cwdm_thermalisation` checks 20-second population conservation and IME
recovery, including zero-population and spin-free blocks. Its shared-identity
comparison isolates the unit-column source difference algebraically; an actual
old-kernel comparison additionally requires a separately captured stock result
and an explicit descriptor embedding. Do not treat the algebraic test as a
stock execution.
