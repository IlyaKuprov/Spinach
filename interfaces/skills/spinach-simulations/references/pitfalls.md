# Failure modes and diagnostics

Almost every kernel function has a `grumble()` input checker that stops with a
specific message, and the no-defaults policy means a missing physical input is
an error, not a guess. Crashes are therefore usually self-explanatory once the
message is matched to its cause. The dangerous failures are the silent ones -
a simulation that runs to completion and produces a plausible but wrong
spectrum. This file covers both, ordered by how often a newcomer hits them,
then scaling and memory, then the recommended validation procedure.

## Crash catalogue

### Calling things in the wrong order

```
basis set information is missing, run basis() before calling this function.
```

The most common beginner error. `state`, `operator`, `superop`, and
`hamiltonian` all return objects expressed in the current basis and refuse to
run before `basis()` has been called. The fix is always the same ordering:
`create` → `basis` → everything else. The related `assumption information is
missing, run assume() before calling this function.` means `assume()` was
skipped; contexts call it internally, so it only appears when bypassing the
context layer.

```
Never run Matlab in script mode. All your variables persist, and
there will be no end to confusion. Add a [function ... end] block.
```

`create` refuses to run from the interactive prompt because stale workspace
variables (an old `inter` with leftover fields) are a classic source of
irreproducible behaviour. Put the simulation in a `function ... end` block,
or run it as a script file under `matlab -batch`.

```
paths have not been correctly set - please follow installation
instructions (make sure you had included the subdirectories).
```

Run `fix_path('add')` from the Spinach root (or add the kernel subdirectories
per the installation instructions) before calling `create`.

### Cell array dimensions

```
both dimensions of inter.coupling.scalar array must match the number of spins.
```

MATLAB grows cell arrays to the smallest size that fits the assigned index, so
`inter.coupling.scalar{1,2}=7.0` alone creates a 1×2 cell, not N×N. Either
preallocate with `inter.coupling.scalar=cell(N,N);` or assign the corner
element (`inter.coupling.scalar{N,N}=0;`) so the array reaches full size. The
same N×N requirement applies to `inter.coupling.matrix` and
`inter.coupling.eigs`; the 1×N requirement applies to `inter.zeeman.scalar`,
`inter.zeeman.matrix`, `inter.coordinates`, and the per-spin relaxation rate
arrays, each with its own exactly-worded message.

### Missing required fields

`create` and the contexts demand every physical input explicitly:

| Message | Fix |
|---|---|
| `sys.isotopes subfield must be present.` | Supply the isotope list; it defines N for everything else. |
| `magnet induction must be specified in sys.magnet varaible.` | Set `sys.magnet` in tesla (the typo is in the source; quote when grepping). |
| `working spins must be specified in parameters.spins field.` | Set `parameters.spins={'1H'}` (or the relevant channel list). |
| `sweep width should be specified in parameters.sweep variable.` | Set `parameters.sweep` in Hz. |
| `number of points should be specified in parameters.npoints variable.` | Set `parameters.npoints`. |
| `initial state must be specified in parameters.rho0 variable.` | Build with `state()` after `basis()`. |
| `detection state must be specified in parameters.coil variable.` | Build with `state()` after `basis()`. |
| `spherical averaging grid must be specified in parameters.grid variable.` | Powder context only; e.g. `parameters.grid='rep_2ang_6400pts_sph'`. |
| `list of decoupled spins (or an empty cell array) should be supplied in parameters.decouple variable.` | Set `parameters.decouple={}` when calling `acquire` directly; contexts default it with a logged warning. |
| `transition moment tolerance must be specified in parameters.tm_tol variable.` | `fieldsweep` also needs `rspt_order` and `int_tol`; copy values from the nearest example. |

When a required-field error surfaces from deep in the call stack, the fastest
fix is to diff your `parameters` structure against the nearest working
example rather than reading the whole stack.

### Basis specification

```
basis specification in bas.formalism is required.
unrecognized formalism - see the basis preparation section of the manual.
unrecognized approximation - see the basis preparation section of the manual.
```

Both `bas.formalism` and `bas.approximation` are mandatory strings. The IK
approximations then have their own required fields and cross-rules:

| Message | Meaning |
|---|---|
| `connectivity tracing depth must be specified in bas.level variable.` | IK-0, IK-1, and IK-DNP need `bas.level` (IK-DNP as a 1x3 vector); IK-2 does not use it. |
| `connectivity type must be specified in bas.connectivity variable.` | IK-1/2 need `'scalar_couplings'` or `'full_tensors'`. |
| `proximity tracing depth must be specified in bas.space_level variable.` | IK-1/2 need `bas.space_level`. |
| `bas.level is only applicable to IK-0,1,2,DNP basis sets.` | Do not set IK fields with `approximation='none'`; the same rule exists for `bas.space_level` and `bas.connectivity` (IK-1/2 only). |
| `bas.approximation should be set to 'none' in zeeman-hilb formalism.` | Hilbert space is always complete. |
| `IK-DNP approximation requires both electrons and nuclei.` | IK-DNP is for electron-nuclear systems only. |
| `multiplicities above 16 are not supported by sphten-liouv formalism.` | Very high spins need a different formalism. |

### Formalism and approximation mismatches

Several features exist only in Liouville space, most only in `sphten-liouv`;
picking `zeeman-hilb` or `zeeman-wavef` for convenience and then asking for
relaxation is a frequent dead end:

```
Redfield relaxation theory is only available for sphten-liouv formalism.
extended T1,T2 relaxation theory is only available for sphten-liouv formalism.
analytical decoupling is only available for sphten-liouv, zeeman-liouv, and zeeman-hilb formalisms.
chemical reaction modelling is only available for sphten-liouv formalism.
bas.projections option is only available for sphten-liouv formalism.
```

(Lindblad, Nottingham, Weizmann, SRFK, SRSK relaxation and IME thermalisation
are refused outside Liouville space with parallel wording; the `kite` and
`secular` values of `inter.rlx_keep` are also `sphten-liouv` only.) The
fix is to move to `sphten-liouv`, which is the correct default anyway.

### Relaxation setup

Requesting a relaxation theory pulls in its own required fields:

| Message | Fix |
|---|---|
| `correlation time(s) must be specified with Redfield theory.` | Set `inter.tau_c` in seconds. |
| `relaxation superoperator term retention policy must be specified in inter.rlx_keep field.` | `'diagonal'`, `'kite'`, `'secular'`, or `'labframe'`. |
| `relaxation destination must be specified in inter.equilibrium variable.` | `'zero'`, `'IME'`, or `'dibari'`. |
| `inter.temperature is required when relaxation has a target.` | Set the temperature in kelvin. |
| `R1 and R2 rates must be specified with extended T1,T2 relaxation theory.` | Set `inter.r1_rates`, `inter.r2_rates` in Hz. |

Physics-guard errors in the same area mean the model is invalid, not that a
field is missing: `T1,2>>tau_c validity condition violation in Redfield
theory` (outside the Redfield validity regime); `Lindblad relaxation theory
forbids R2 < 0.5*R1 on thermodynamic grounds.` (hard bound on input rates);
`numerical accuracy issues, temperature too low - switch to Hilbert space.`
(Liouville-space thermal state underflows at very low temperature).

### State and spin specification

```
the requested state is not present in the basis.
```

The state you asked `state()` for was excluded by the basis restriction -
typically a spin-correlation order outside `bas.level`, or a projection
excluded by `bas.projections`. Raise the restriction or fix the state
specification. Related messages: `spins and operators cell arrays should have
the same number of elements.`, `parameters.spins refers to a spin that is not
present in the system.`.

Spin indices always refer to positions in `sys.isotopes`. After importing a
system (from `g2spinach` or similar), print `sys.isotopes` and confirm which
index is which before writing index-based parameters - a spin list pointing
at the wrong index runs cleanly and computes the wrong experiment.

### Chemical kinetics

In exchange simulations each chemical species is a separate spin system:
`couplings detected between spins in different chemical species.` means the
coupling arrays connect spins across species, and `concentrations
(inter.chem.concs) must be provided.` means the required concentrations are
missing.

## Silent wrong-result traps

These produce a finished, plausible-looking spectrum. Check for them by
construction, not by inspection of the output alone.

**Unit mistakes.** Chemical shifts are ppm, all couplings are Hz, Euler angles
are radians, times are seconds (full table in the main skill file). Classic
slips: a shift entered in Hz into `inter.zeeman.scalar`; a hyperfine coupling
in gauss passed through `gauss2mhz` without converting its MHz output to Hz, or
a millitesla value not converted with `mt2hz`; degrees instead of radians in
Euler angles, which distorts powder patterns smoothly - no crash, wrong
lineshape.

**Doubled couplings.** `create` adds `inter.coupling.scalar{i,j}` and
`inter.coupling.scalar{j,i}` independently into the coupling matrix.
Specifying the same coupling in both triangles doubles it. Enter each coupling
once; import scripts that fill both directions must halve the values.

**Wrong assumptions string.** The `assumptions` argument decides which
Hamiltonian terms survive the rotating-frame treatment. `'nmr'` applied to an
electron-nuclear problem silently discards the electron physics; `'labframe'`
omitted where needed silently keeps the rotating-frame approximation. The log
prints the full assumption set applied (e.g. "generic high-field NMR
assumption set" followed by the retained terms) - read it and check it matches
the physics intended.

**Sign of the Liouvillian.** The convention is `L=H+1i*R+1i*K` propagated as
`exp(-1i*L*t)`. Flipping either sign turns dissipation into exponential
growth. Symptom: a FID whose magnitude grows with time, or relaxation that
heats instead of cools. Contexts get this right; hand-assembled Liouvillians
in custom sequences are where it goes wrong.

**Aliasing.** Signals outside the window `parameters.offset ± sweep/2` fold
back into the spectrum at wrong positions - peaks at plausible but incorrect
shifts, often mirrored about a window edge. The log prints the applied offset
per channel and each Zeeman term in Hz; every printed frequency must lie
inside ±sweep/2 after the offset is applied.

**Warnings that are really defaults.** Spinach logs `WARNING` lines whenever
it assumes something: `WARNING - spin temperature not specified, assuming 298
Kelvin`, `WARNING - no relaxation theory specified`, `WARNING - no coordinates
given, point dipolar interactions assumed to be zero.` Each is a physical
assumption - a solid-state simulation without coordinates has no dipolar
couplings at all and happily produces a liquid-like spectrum. Grep the log
for `WARNING` and justify every line.

**Axis conventions in plotting.** `plot_1d` defaults `parameters.axis_units`
to `'ppm'` and `parameters.invert_axis` to the NMR tradition (axis reversed)
when unset, and logs that it did. A mirrored spectrum, or peaks at the
negatives of the expected shifts, is usually an axis convention mismatch, not
a physics error - fix the plot parameters before touching the simulation.

**Unconverged results.** An insufficient `bas.level`, a too-coarse powder
grid, or a too-long timestep does not crash; it produces a smooth, wrong
spectrum (missing multiplet components, distorted powder lineshapes, phase
errors). Convergence is established by refinement, never by appearance - see
the validation procedure below.

## Scaling and memory

The Liouville-space dimension grows as 4^N for N spin-1/2 particles (2^N in
Hilbert space). A complete `sphten-liouv` basis is practical to roughly ten
spins; beyond that, the basis restriction is the tool, not a bigger machine:

- `IK-1`/`IK-2` handles large organic molecules routinely. For `IK-1`, converge
  by incrementing `bas.level` and, when relevant, `bas.space_level`; for
  `IK-2`, `bas.level` does not control basis construction, so increment
  `bas.space_level` instead until the spectrum stops changing.
- `bas.longitudinals` and `bas.projections` are physical approximations, not
  free reductions: the former deletes transverse states on selected spins,
  while the latter deletes total-coherence blocks that are not retained. Use
  them only when the initial state, pulse sequence, Hamiltonian, relaxation,
  and observable cannot enter the discarded blocks, and validate against an
  unfiltered basis when practical.
- `bas.sym_group`/`bas.sym_spins` exploit permutation symmetry - large
  savings for methyl groups and symmetric aromatics.

At run time Spinach reduces dimension further on its own: zero-track
elimination and path tracing routinely cut the active space by an order of
magnitude, visible in the log as "state space dimension reduced from 64 to
15". These reductions can be switched off through
`sys.disable={'zte','pt','symmetry',...}` for debugging, at a large cost.

Matrices switch to sparse algebra automatically, and above a state-space
dimension of 10000 (a tunable tolerance) Spinach uses Krylov propagation
rather than forming the propagator. When writing custom sequences for large
systems, use `krylov` or `step` - both apply the action of the matrix
exponential without materialising `expm(-1i*L*dt)`, which is what exhausts
memory first. GPU arithmetic is enabled with `sys.enable={'gpu'}`.
Powder-grid cost is linear in the number of orientations and parallelises
over the pool that `create` starts.

## The headless validation harness

Run every unattended simulation through this pattern:

1. Write a self-contained script that sets the path explicitly, runs the
   simulation, prints machine-checkable invariants, and ends with a unique
   success marker that appears nowhere else:

   ```matlab
   cd('/path/to/Spinach'); fix_path('add');
   % ... simulation ...
   fprintf('FID length: %d\n',numel(fid));
   fprintf('FID max abs: %.6e\n',max(abs(fid)));
   fprintf('FID finite: %d\n',all(isfinite(fid)));
   disp('SPINACH_SMOKE_OK');
   ```

2. Run it non-interactively and capture everything:

   ```bash
   matlab -batch "run('/path/to/probe.m')" > probe.log 2>&1
   ```

3. Judge success on two conditions together: the unique marker is present in
   the log, and the log contains no error text. Exit status alone is not
   evidence - MATLAB can exit zero after a partial run.

4. Read the log as physics, not as noise. Spinach narrates the entire
   construction: every Zeeman and coupling tensor with its trace and
   anisotropy norms, the assumption set actually applied, every Hamiltonian
   term with its coefficient in Hz, the basis dimension before and after
   reduction, and every `WARNING` default. Most silent wrong-result traps
   above are visible in this narration before the spectrum is ever plotted.

If the script must plot, headless rendering needs
`set(0,'DefaultFigureVisible','off')`, `opengl('save','software')`, and
`QT_QPA_PLATFORM=offscreen`; a crash inside figure creation under `-batch` is
an environment problem, not a Spinach one. For non-plotting runs, `-nojvm`
is faster.

## What counts as evidence that a simulation is right

The four checks from the main skill file, as concrete procedures:

**1. Positions match theory.** Compute expected peak positions by hand before
looking at the spectrum: shift in ppm times the Larmor frequency in MHz gives
the offset in Hz; multiplet splittings must equal the input J values exactly;
for EPR, the centre field follows from the g-value and microwave frequency.
Then read positions off the axis vector and compare numerically, not
visually. The log's per-spin Zeeman terms in Hz are the same numbers from the
other end and must agree.

**2. Invariants are intact.** `all(isfinite(fid))` must be true. A detected FID
is not generally monotonic: coherent transfer, J modulation, echoes, and an
observable initially orthogonal to the state can make its magnitude grow even
with dissipative relaxation. Require bounded or monotonic decay only when the
specific dynamics guarantee it; thermal states must have positive populations,
and in Hilbert-space formalisms the density matrix trace is conserved. Print
these checks in the script so they are run every time, not once.

**3. The result is converged.** Rerun with the active basis restriction
parameter incremented (`bas.level` for `IK-1`, `bas.space_level` for `IK-2`),
the powder grid one rank finer, and the time step halved, one factor at a time,
and compare spectra numerically - `norm(s1-s2)/norm(s1)` below a stated
tolerance, not "looks the same". A converged result stops moving under
refinement; an unconverged one can look entirely reasonable.

**4. A limiting case reproduces a known answer.** A single spin must give one
line at its offset; all J set to zero must collapse multiplets to singlets;
weak coupling must give first-order patterns with Pascal's-triangle
intensities. Best practice: run the nearest `examples/` file unchanged first
and confirm its documented output, then introduce changes one at a time -
when something breaks, the cause is in the last change.

A simulation is trustworthy when all four hold and every `WARNING` line in the
log is accounted for. Anything less is a run, not a result.
