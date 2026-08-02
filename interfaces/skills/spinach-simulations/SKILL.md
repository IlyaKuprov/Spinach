---
name: spinach-simulations
description: Expert use of the Spinach library for magnetic resonance simulation - NMR, EPR, DNP, MRI, and optimal control. Covers spin system setup, basis selection, simulation contexts, pulse sequences, relaxation, and validation. Use whenever a spin dynamics or magnetic resonance simulation is requested.
---

# Spinach simulations

Spinach is a MATLAB library for spin dynamics simulation, covering liquid- and
solid-state NMR, EPR, DNP, MRI, and optimal control. This skill is about *using*
Spinach to answer physical questions. Editing the Spinach source is a different
task with different rules: see `AGENTS.md` in the repository root.

## Loading Spinach

`fix_path.m` in the repository root puts Spinach on the MATLAB path; run it
from the root directory. `fix_path()` resets the MATLAB path to default before
adding Spinach - the right call on a machine that may carry conflicting
toolboxes; `fix_path('add')` adds Spinach to the existing path;
`fix_path('remove')` takes it off again. It finishes with existential checks
and prints `Spinach is ready to run.` The idiom for unattended work:

```matlab
cd('<spinach-root>'); fix_path('add');
```

## The two rules that matter most

**1. Never invent physical parameters.** Spinach has a deliberate policy of never
guessing: there are no default values, and a missing input is an error rather
than an assumption. The same discipline applies to you. Chemical shifts,
couplings, g-tensors, and correlation times come from the literature, from
experiment, or from a quantum chemistry calculation - never from plausibility.
If a parameter is unavailable, say so instead of inventing a number.

**2. Start from the nearest example.** The `examples/` tree holds over 700
working simulations across 37 problem areas, and almost every request resembles
one of them. Locating the closest example and adapting it is faster and far more
reliable than writing from scratch, because the example already encodes the
correct basis, context, assumptions, and processing chain for that physics.

```bash
ls examples/                                  # 37 problem areas
grep -rl "hsqc" examples/ | head              # find by experiment
grep -rl "sys.isotopes={'E'" examples/ | head # find by system type
```

## The canonical simulation

Every Spinach simulation has the same seven-part shape. This is a real, working
liquid-state NMR simulation; the section comments are the house convention and
should be kept.

```matlab
% Magnetic induction, Tesla
sys.magnet=5.9;

% Isotopes and interactions
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={1.0 1.5};
inter.coupling.scalar=cell(2,2);
inter.coupling.scalar{1,2}=7.0;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Sequence parameters
parameters.spins={'1H'};
parameters.rho0=state(spin_system,'L+','1H');
parameters.coil=state(spin_system,'L+','1H');
parameters.decouple={};
parameters.offset=300;
parameters.sweep=300;
parameters.npoints=1024;
parameters.zerofill=4096;
parameters.axis_units='Hz';
parameters.invert_axis=1;

% Simulation
fid=liquid(spin_system,@acquire,parameters,'nmr');

% Apodisation and Fourier transform
fid=apodisation(spin_system,fid,{{'exp',10}});
spectrum=fftshift(fft(fid,parameters.zerofill));

% Plotting
kfigure(); plot_1d(spin_system,real(spectrum),parameters);
```

The order is not negotiable. `create` absorbs and validates the physics,
`basis` chooses the state space, and only then do `state` and `operator`
mean anything, because they return objects expressed in that basis. Calling
`state` before `basis` is the single most common beginner error and produces
`basis set information is missing, run basis() before calling this function`.

On output, the frequency axis is `ft_axis(parameters.offset,parameters.sweep,
parameters.zerofill)` and runs in absolute frequency: after `fftshift(fft(...))`
of an `L+`-detected signal, a spin with chemical shift delta lands at
delta times the base frequency, in a window centred at `parameters.offset`.
Query the base frequency from Spinach itself rather than external constants:
`[gamma,~]=spin('1H')` returns the magnetogyric ratio in rad/(s*T), so
`gamma*sys.magnet/(2e6*pi)` is the Hz-per-ppm conversion factor. In the
apodisation call, `{'exp',k}` multiplies the FID by `exp(-k*x)` with `x`
running from 0 to 1 across the FID, i.e. k e-folds of decay by the last
point. `kfigure` is an interactive figure helper; for unattended runs, plot
to an invisible figure and `print`, or skip plotting - see the headless
harness in `references/pitfalls.md`.

## Units on input

Getting these wrong produces a plausible spectrum that is quantitatively wrong,
which is far more dangerous than a crash. Spinach converts everything to rad/s
internally, but the *input* units are fixed by convention:

| Quantity | Field | Unit on input |
|---|---|---|
| Magnetic induction | `sys.magnet` | tesla |
| Nuclear chemical shift | `inter.zeeman.scalar/eigs/matrix` | ppm |
| Electron g-tensor | `inter.zeeman.*` for `'E'` spins | dimensionless g-value |
| All couplings (J, dipolar, hyperfine, quadrupolar, ZFS) | `inter.coupling.*` | hertz |
| Coordinates | `inter.coordinates` | angstrom |
| Correlation time | `inter.tau_c` | seconds |
| Temperature | `inter.temperature` | kelvin |
| Relaxation rates | `inter.r1_rates`, `inter.r2_rates` | hertz |
| Chemical exchange rates | `inter.chem.rates` | hertz |
| Offsets, sweeps, J in sequences | `parameters.offset/sweep/J` | hertz |
| Euler angles | `inter.*.euler`, `parameters.orientation` | radians |
| Times | `parameters.tau`, `tmix`, `timestep` | seconds |
| Control power levels | `control.pwr_levels` | rad/s |

Two conversions catch people out: hyperfine couplings quoted in gauss need
`gauss2mhz(...)` followed by multiplication by `1e6` to convert its MHz output
to hertz, while couplings quoted in millitesla can use `mt2hz(...)` directly.
A coupling tensor entered as a matrix is still in hertz even when its
magnitude looks like a frequency in rad/s.

## Choosing the state space

`bas.formalism` chooses the mathematical space; `bas.approximation` chooses how
much of it to keep. Together they decide whether a simulation is feasible.

| `bas.formalism` | Use when |
|---|---|
| `sphten-liouv` | Default for almost everything. Liouville space in irreducible spherical tensors; required for relaxation, chemical kinetics, and state-space restriction. |
| `zeeman-hilb` | Small systems needing exact Hilbert-space treatment; field-swept EPR; propagator inspection. |
| `zeeman-liouv` | Liouville space in the Zeeman basis; occasional diagnostic and dissipative work. |
| `zeeman-wavef` | Wavefunction (state-vector) dynamics; no relaxation. |

| `bas.approximation` | Meaning |
|---|---|
| `none` | Complete basis. Correct by construction, but the state space grows as 4^N for spin-1/2 in Liouville space; practical to roughly ten spins. |
| `IK-0` | Keeps all states up to a given spin correlation order (`bas.level`), irrespective of distance. |
| `IK-1` | Correlation order `bas.level` restricted to spins within `bas.space_level` bonds; needs `bas.connectivity`. The workhorse for large molecules. |
| `IK-2` | Uses direct coupling connectivity with proximity subgraphs controlled by `bas.space_level`; standard for strychnine-class organic molecules. |
| `IK-DNP` | Tailored to electron-nuclear DNP systems. |

`bas.connectivity` is `'scalar_couplings'` or `'full_tensors'`. The filters
`bas.longitudinals` and `bas.projections` are physical approximations:
`bas.longitudinals` deletes transverse states on selected spins, while
`bas.projections` deletes total-coherence blocks that are not retained. Use
them only when the initial state, pulse sequence, Hamiltonian, relaxation, and
observable cannot enter the discarded blocks. Permutation symmetry via
`bas.sym_group` and `bas.sym_spins` factorises the problem into irreducible
representations and can be a large saving for methyl groups and symmetric
aromatics.

## Choosing the context

The context supplies the physical setting: it builds the Hamiltonian,
relaxation and kinetics superoperators, performs any orientational or spatial
averaging, and calls your pulse sequence. All contexts share the signature
`answer=context(spin_system,@pulse_sequence,parameters,assumptions)`, except
`imaging` and `meshflow`, which take no assumptions argument.

| Context | Physical situation |
|---|---|
| `liquid` | Isotropic tumbling; no orientational averaging. |
| `crystal` | Single orientation, given by `parameters.orientation`. |
| `powder` | Static powder; averages over `parameters.grid`. |
| `singlerot` | Magic angle and other single-axis spinning, Fokker-Planck. |
| `doublerot` | Double rotation with two rotors. |
| `floquet` | Spinning treated by Floquet theory. |
| `gridfree` | Spinning and the stochastic Liouville equation without a grid. |
| `imaging` | Spatial dynamics with gradients, diffusion, and flow. |
| `meshflow` | Flow and reaction on an imported finite-element mesh. |

The `assumptions` string selects which terms survive the rotating-frame
treatment inside `hamiltonian`: `'nmr'`, `'esr'`, `'deer'`, `'deer-zz'`,
`'labframe'`, `'qnmr'`, and the DNP variants `'se_dnp_h+'`, `'se_dnp_h-'`,
`'se_dnp_h0'`. Using `'nmr'` for an EPR problem silently discards the physics
you wanted.

## Propagation

Spinach propagates with `exp(-1i*L*t)`, and the Liouvillian is assembled as

```matlab
L=H+1i*R+1i*K;
```

with `R` from `relaxation` and `K` from `kinetics`. The factors of `1i` are
not decoration: they make relaxation and chemical exchange dissipative rather
than oscillatory, and getting the sign wrong produces exponential growth.
Contexts do this for you; you only assemble `L` by hand when writing a pulse
sequence or bypassing the context system.

`evolution` is the general time-propagation routine, with output modes
`'final'`, `'trajectory'`, `'total'`, `'refocus'`, `'observable'`, and
`'multichannel'`. For very large systems where the propagator will not fit in
memory but the Liouvillian will, use `krylov`; `step` applies a matrix
exponential without forming it.

## Reference material

Consult these for detail; they are dense and specific.

- `references/inputs.md` - complete `sys`, `inter`, and `bas` field reference, importing from Gaussian, ORCA, CASTEP, PDB/BMRB, and GISSMO.
- `references/contexts.md` - context and pulse-sequence API, `parameters` fields per context, `state` and `operator` grammar, 2D acquisition and processing.
- `references/recipes.md` - starting points by physical problem, across NMR, EPR, DNP, MRI, and optimal control.
- `references/relaxation.md` - relaxation theories, thermal equilibrium, chemical kinetics, powder grids, and the convention set.
- `references/pitfalls.md` - failure modes, diagnostics, and what counts as evidence that a simulation is right.
- `references/literature.md` - the citation for Spinach, the methodology papers behind its algorithms, and published applications.

## Validating a simulation

A simulation that runs is not a simulation that is right. Before reporting a
result, check what physics says it must satisfy.

- Does the spectrum appear where theory puts it? Chemical shifts should land at
  their known positions, multiplet splittings should equal the input couplings
  in hertz, and EPR lines should sit at the field the g-value implies.
- Are the invariants intact? Populations should be positive, the trace of a
  density matrix preserved where it must be, and the signal finite everywhere;
  a detected FID need not be monotonic or bounded by its initial magnitude
  when coherent transfer contributes.
- Does it converge? Increase the active basis restriction (`bas.level` and,
  when relevant, `bas.space_level` for `IK-1`; `bas.space_level` for `IK-2`),
  refine the grid, and halve the time step; a converged result stops moving.
  An unconverged simulation can look entirely reasonable.
- Does a limiting case reproduce a known answer? Weak coupling should give
  first-order multiplets, a single spin should give one line, and a
  well-studied system should reproduce its published spectrum.

Run MATLAB non-interactively for unattended work, and require a unique success
marker at the end of the script rather than trusting exit status alone.
