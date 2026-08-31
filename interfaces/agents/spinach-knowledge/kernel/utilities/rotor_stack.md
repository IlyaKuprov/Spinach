# kernel/utilities/rotor_stack.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/rotor_stack.m`
- Signature: `[L,rotor_phases]=rotor_stack(spin_system,parameters,assumptions)`
- Total lines: 246

## Purpose

Returns a rotor stack of Liouvillians or Hamiltonians. The stack is needed for the traditional style calculation of MAS dynamics. Syntax: L=rotor_stack(spin_system,parameters,assumptions)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 65-66: Check consistency; implemented by `grumble(spin_system,parameters,assumptions)`.
- Lines 68-69: Get the Hamiltonian; implemented by `[H,Q]=hamiltonian(assume(spin_system,assumptions))`.
- Lines 71-72: Apply offsets; implemented by `H=frqoffset(spin_system,H,parameters)`.
- Lines 74-77: Get rotor axis orientation; implemented by `[rotor_phi,rotor_theta,~]=cart2sph(parameters.axis(1), parameters.axis(2), parameters.axis(3))`.
- Lines 80-81: Compute rotor angles; implemented by `rotor_phases=fourdif(2*parameters.max_rank+1,1)`.
- Lines 83-84: Get carrier operators; implemented by `C=cell(size(parameters.rframes))`.
- Lines 89-90: Preallocate Liouvillian blocks; implemented by `L=cell(2*parameters.max_rank+1,1)`.
- Lines 92-93: Build Liouvillian blocks; implemented by `parfor n=1:(2*parameters.max_rank+1)`.
- Lines 95-96: Get the block started; implemented by `L{n}=H`.
- Lines 98-99: Loop over spherical ranks; implemented by `for r=1:numel(Q)`.
- Lines 101-102: Compute rotor axis tilt; implemented by `D_rot2lab=wigner(r,+rotor_phi,+rotor_theta,0)`.
- Lines 105-108: Compute the initial crystallite orientation; implemented by `D_initial=wigner(r,parameters.orientation(1), parameters.orientation(2), parameters.orientation(3))`.
- Lines 110-111: Compute rotor rotation; implemented by `D_rotor=wigner(r,0,0,rotor_phases(n))`.
- Lines 113-114: Decide the reference frame; implemented by `if strcmp(parameters.masframe,'magnet')`.
- Lines 116-117: Compose rotations; implemented by `D=D_rot2lab*D_rotor*D_lab2rot*D_initial`.
- Lines 121-122: Compose rotations; implemented by `D=D_rot2lab*D_rotor*D_initial`.
- Lines 126-127: Complain and bomb out; implemented by `D=0; error('unknown MAS frame.')`.
- Lines 131-132: Build the block; implemented by `for k=1:(2*r+1)`.

### Control flow inferred from the code

- Line 85: `for` loop over `n=1:numel(parameters.rframes)`.
- Line 93: `parfor` loop over `n=1:(2*parameters.max_rank+1)`.
- Line 99: `for` loop over `r=1:numel(Q)`.
- Line 114: conditional branch on `strcmp(parameters.masframe,'magnet')`.
- Line 132: `for` loop over `k=1:(2*r+1)`.
- Line 133: `for` loop over `m=1:(2*r+1)`.
- Line 141: `for` loop over `k=1:numel(parameters.rframes)`.

### Key state/data transformations

- Lines 69: computes `[H,Q]` using `[H,Q]=hamiltonian(assume(spin_system,assumptions))`.
- Lines 72: computes `H` using `H=frqoffset(spin_system,H,parameters)`.
- Lines 75-77: computes `[rotor_phi,rotor_theta,~]` using `[rotor_phi,rotor_theta,~]=cart2sph(parameters.axis(1), parameters.axis(2), parameters.axis(3))`.
- Lines 78: computes `rotor_theta` using `rotor_theta=pi/2-rotor_theta`.
- Lines 81: computes `rotor_phases` using `rotor_phases=fourdif(2*parameters.max_rank+1,1)`.
- Lines 84: computes `C` using `C=cell(size(parameters.rframes))`.
- Lines 86: computes `C{n}` using `C{n}=carrier(spin_system,parameters.rframes{n}{1})`.
- Lines 90: computes `L` using `L=cell(2*parameters.max_rank+1,1)`.
- Lines 96: computes `L{n}` using `L{n}=H`.
- Lines 102: computes `D_rot2lab` using `D_rot2lab=wigner(r,+rotor_phi,+rotor_theta,0)`.
- Lines 103: computes `D_lab2rot` using `D_lab2rot=wigner(r,0,-rotor_theta,-rotor_phi)`.
- Lines 106-108: computes `D_initial` using `D_initial=wigner(r,parameters.orientation(1), parameters.orientation(2), parameters.orientation(3))`.
- Lines 111: computes `D_rotor` using `D_rotor=wigner(r,0,0,rotor_phases(n))`.
- Lines 117: computes `D` using `D=D_rot2lab*D_rotor*D_lab2rot*D_initial`.

### Local helper functions

- Line 154: `grumble()` — `function grumble(spin_system,parameters,assumptions)`. Spinning axis
  - Representative operation: `if ~isfield(parameters,'axis')`.
  - Representative operation: `error('parameters.axis subfield must be present.')`.

## Parameters / inputs

- parameters.axis -spinning axis, given as a normalized
- 3-element vector
- parameters.offset -a cell array giving transmitter off-
- sets in Hz on each of the spins listed
- in parameters.spins array
- parameters.spins -a cell array giving the spins that
- the offsets refer to, e.g. {'1H','13C'}
- parameters.max_rank -maximum harmonic rank to retain in
- the solution (increase till conver-
- gence is achieved, approximately
- equal to the number of spinning si-
- debands in the spectrum)
- parameters.rframes -rotating frame specification, e.g.
- {{'13C',2},{'14N,3}} requests second
- order rotating frame transformation
- with respect to carbon-13 and third
- order rotating frame transformation
- with respect to nitrogen-14. When
- this option is used, the assumptions
- on the respective spins should be
- laboratory frame.
- parameters.orientation -the orientation of the spin system
- at rotor phase zero, a vector of
- three Euler angles in radians.
- parameters.masframe -the frame in which the rotations
- are applied. The possibilities are:
- 'magnet' -the initial orientation in the lab frame
- (three-angle powder grids will be required)
- 'rotor' -the initial orientation in the rotor frame
- (two-angle powder grids will be required)
- assumptions -assumption set to be used in generating the
- Hamiltonian, see assume.m

## Outputs

- L -a cell array of Hamiltonian or Liouvillian matrices,
- one for each tick of the rotor.
- rotor_phases -rotor phases at each tick, radians
- Note: relaxation and chemical kinetics are not included.

## Implementation structure

- Returns a rotor stack of Liouvillians or Hamiltonians. The stack is
- needed for the traditional style calculation of MAS dynamics. Syntax:
- L=rotor_stack(spin_system,parameters,assumptions)
- parameters.axis -spinning axis, given as a normalized
- 3-element vector
- parameters.offset -a cell array giving transmitter off-
- sets in Hz on each of the spins listed
- in parameters.spins array
- parameters.spins -a cell array giving the spins that
- the offsets refer to, e.g. {'1H','13C'}
- parameters.max_rank -maximum harmonic rank to retain in
- the solution (increase till conver-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `hamiltonian()`, `assume()`, `frqoffset()`, `cart2sph()`, `fourdif()`, `carrier()`, `wigner()`, `rotor_phases()`, `strcmp()`, `rotframe()`, `clean_up()`, `isfield()`, `elseif()`, `isrow()`, `iscell()`.
