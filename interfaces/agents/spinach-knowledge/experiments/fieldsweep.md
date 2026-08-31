# experiments/fieldsweep.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/fieldsweep.m`
- Signature: `[spec,parameters]=fieldsweep(spin_system,parameters)`
- Total lines: 216

## Purpose

Field swept powder EPR spectra. A rough implementation with ex- pensive eigenfields algorithm, an explicit spherical grid, and a hard-coded Lorentzian line shape. Syntax: [b_axis,spec]=fieldsweep(spin_system,parameters)

## Physical / mathematical content

- This file belongs to the `experiments` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 58-59: Check consistency; implemented by `grumble(spin_system,parameters)`.
- Lines 61-62: Peak position tolerance is quarter a grid pixel; implemented by `window_size=max(parameters.window)-min(parameters.window)`.
- Lines 65-66: Get the Hamiltonians and tidy up their isotropic parts; implemented by `[Ic,Qc]=hamiltonian(assume(spin_system,'labframe','couplings'))`.
- Lines 70-71: Get the microwave operator (state in Liouville space); implemented by `Hmw=state(spin_system,'Lx',parameters.spins{1})`.
- Lines 73-75: Get the initial grid and compute its convex hull; implemented by `load([spin_system.sys.root_dir filesep 'kernel' filesep 'grids' filesep parameters.grid],'betas','gammas')`.
- Lines 78-81: Make the magnetic field axis; implemented by `parameters.b_axis=linspace(parameters.window(1), parameters.window(2), parameters.npoints)`.
- Lines 83-84: Preallocate asynchronous work outputs; implemented by `eigensets(1:grid_size)=parallel.Future`.
- Lines 86-87: Schedule work; implemented by `for n=1:grid_size`.
- Lines 89-90: Create a local copy and specify system orientation; implemented by `localpar=parameters; localpar.orientation=[0 betas(n) gammas(n)]`.
- Lines 92-93: Assemble Zeeman and coupling Hamiltonians; implemented by `Hz=Iz+orientation(Qz,localpar.orientation)`.
- Lines 96-97: Submit asynchronous eigenfields calculation job for this vertex; implemented by `eigensets(n)=parfeval(@eigenfields,1,spin_system,localpar,Hz,Hc,Hmw)`.
- Lines 101-102: Retrieve vertex eigensets; implemented by `eigensets=fetchOutputs(eigensets)`.
- Lines 104-105: Add coordinates; implemented by `for n=1:grid_size`.
- Lines 111-112: Preallocate asynchronous work outputs; implemented by `spec(1:size(hull,1))=parallel.Future`.
- Lines 114-115: Schedule work; implemented by `for n=1:size(hull,1)`.
- Lines 117-118: Extract triangle vertex indices; implemented by `a=hull(n,1); b=hull(n,2); c=hull(n,3)`.
- Lines 120-121: Extract vertex eigensets; implemented by `triangle=eigensets([a b c])`.
- Lines 123-125: Call recursive Voitlander integrator; implemented by `spec(n)=parfeval(@voitlander,1,spin_system,parameters, triangle,Ic,Iz,Qc,Qz,Hmw)`.

### Control flow inferred from the code

- Line 87: `for` loop over `n=1:grid_size`.
- Line 105: `for` loop over `n=1:grid_size`.
- Line 115: `for` loop over `n=1:size(hull,1)`.

### Key state/data transformations

- Lines 62: computes `window_size` using `window_size=max(parameters.window)-min(parameters.window)`.
- Lines 63: computes `parameters.pp_tol` using `parameters.pp_tol=0.25*window_size/(parameters.npoints-1)`.
- Lines 66: computes `[Ic,Qc]` using `[Ic,Qc]=hamiltonian(assume(spin_system,'labframe','couplings'))`.
- Lines 67: computes `[Iz,Qz]` using `[Iz,Qz]=hamiltonian(assume(spin_system,'labframe','zeeman'))`.
- Lines 68: computes `Ic` using `Ic=(Ic+Ic')/2; Iz=(Iz+Iz')/2`.
- Lines 71: computes `Hmw` using `Hmw=state(spin_system,'Lx',parameters.spins{1})`.
- Lines 76: computes `hull` using `hull=get_hull(betas,gammas); grid_size=numel(betas)`.
- Lines 79-81: computes `parameters.b_axis` using `parameters.b_axis=linspace(parameters.window(1), parameters.window(2), parameters.npoints)`.
- Lines 84: computes `eigensets(1:grid_size)` using `eigensets(1:grid_size)=parallel.Future`.
- Lines 90: computes `localpar` using `localpar=parameters; localpar.orientation=[0 betas(n) gammas(n)]`.
- Lines 93: computes `Hz` using `Hz=Iz+orientation(Qz,localpar.orientation)`.
- Lines 94: computes `Hc` using `Hc=Ic+orientation(Qc,localpar.orientation)`.
- Lines 97: computes `eigensets(n)` using `eigensets(n)=parfeval(@eigenfields,1,spin_system,localpar,Hz,Hc,Hmw)`.
- Lines 102: computes `eigensets` using `eigensets=fetchOutputs(eigensets)`.
- Lines 106: computes `eigensets(n).xyz` using `eigensets(n).xyz=[sin(betas(n))*cos(gammas(n))`.
- Lines 112: computes `spec(1:size(hull,1))` using `spec(1:size(hull,1))=parallel.Future`.
- Lines 118: computes `a` using `a=hull(n,1); b=hull(n,2); c=hull(n,3)`.
- Lines 121: computes `triangle` using `triangle=eigensets([a b c])`.

### Local helper functions

- Line 135: `grumble()` — `function grumble(spin_system,parameters)`.
  - Representative operation: `if spin_system.inter.magnet~=1`.
  - Representative operation: `error('field swept experiment, set sys.magnet=1')`.

## Parameters / inputs

- parameters.grid -initial spherical grid, ideally
- a non-symmetric one to avoid
- transition degeneracies, a good
- start is
- 'rep_2ang_100pts_sph'
- parameters.spins -a one-element cell array speci-
- fying the spin that is coupled
- couple to the microwave field,
- e.g. {'E'}
- parameters.mw_freq -microwave frequency, Hz
- parameters.fwhm -Lorentzian line FWHM, Tesla
- parameters.window -field sweep window in Tesla,
- as a vector [Bmin Bmax]
- parameters.npoints -number of points in the sweep
- parameters.tm_tol -relative transition moment to-
- lerance, 0.01 is a good start
- parameters.rspt_order -perturbation theory order for
- eigenfields calculation, 2 is
- a good start; specify Inf for
- exact diagonalisation
- parameters.int_tol -powder integration tolerance,
- a balance between speed and
- integration accuracy

## Outputs

- b_axis -magnetic field axis for plotting
- spec -field-swept EPR spectrum
- Note: irrespective of the actual sweep extents, the magnetic field
- in sys.magnet should be set to 1 Tesla.
- Note: this experiment should be called directly without a context.

## Implementation structure

- Field swept powder EPR spectra. A rough implementation with ex-
- pensive eigenfields algorithm, an explicit spherical grid, and
- a hard-coded Lorentzian line shape. Syntax:
- [b_axis,spec]=fieldsweep(spin_system,parameters)
- parameters.grid - initial spherical grid, ideally
- a non-symmetric one to avoid
- transition degeneracies, a good
- start is
- 'rep_2ang_100pts_sph'
- parameters.spins - a one-element cell array speci-
- fying the spin that is coupled
- couple to the microwave field,

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `hamiltonian()`, `assume()`, `state()`, `load()`, `get_hull()`, `eigensets()`, `betas()`, `gammas()`, `orientation()`, `parfeval()`, `fetchOutputs()`, `spec()`, `hull()`, `isfield()`, `ischar()`.
