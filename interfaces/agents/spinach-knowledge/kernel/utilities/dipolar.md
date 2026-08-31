# kernel/utilities/dipolar.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/dipolar.m`
- Signature: `spin_system=dipolar(spin_system)`
- Total lines: 204

## Purpose

Computes dipolar couplings in the presence or absence of periodic boundary conditions. This is an auxiliary function of Spinach ker- nel, direct calls are discouraged. Use xyz2dd and xyz2hfc to con- vert Cartesian coordinates into dipolar and hyperfine couplings respectively. Syntax: spin_system=dipolar(spin_system)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 28-29: Check consistency; implemented by `grumble(spin_system)`.
- Lines 31-33: Report to the user; implemented by `report(spin_system,['dipolar interaction distance threshold: ' num2str(spin_system.tols.prox_cutoff) ' Angstrom.'])`.
- Lines 36-37: Preallocate distance vector array; implemented by `distvects=cell(spin_system.comp.nspins,spin_system.comp.nspins)`.
- Lines 39-40: Loop over chemical species; implemented by `for m=1:numel(spin_system.chem.parts)`.
- Lines 42-43: Extract the spin list; implemented by `spin_list=spin_system.chem.parts{m}`.
- Lines 45-46: Loop over atom pairs; implemented by `for n=spin_list`.
- Lines 49-51: Only proceed if coordinates are specified; implemented by `if (~isempty(spin_system.inter.coordinates{n}))&& (~isempty(spin_system.inter.coordinates{k}))&&(n~=k)`.
- Lines 53-54: Determine possible distance vectors from n to k; implemented by `if isempty(spin_system.inter.pbc)`.
- Lines 56-58: Compute the distance vector for standalone system; implemented by `dv=spin_system.inter.coordinates{k}- spin_system.inter.coordinates{n}`.
- Lines 62-63: Preallocate distance vector array; implemented by `nvecs=2*spin_system.tols.dd_ncells+1; dv=zeros(nvecs,3)`.
- Lines 65-66: Compute distance vectors with 1D periodic boundary; implemented by `linear_index=1`.
- Lines 76-77: Preallocate distance vector array; implemented by `nvecs=(2*spin_system.tols.dd_ncells+1)^2; dv=zeros(nvecs,3)`.
- Lines 79-80: Compute distance vectors with 2D periodic boundary; implemented by `linear_index=1`.
- Lines 93-94: Preallocate distance vector array; implemented by `nvecs=(2*spin_system.tols.dd_ncells+1)^3; dv=zeros(nvecs,3)`.
- Lines 96-97: Compute distance vectors with 3D periodic boundary; implemented by `linear_index=1`.
- Lines 113-114: Complain and bomb out; implemented by `error('PBC translation vector array has invalid dimensions.')`.
- Lines 118-119: Ignore distance vectors that are longer than the threshold; implemented by `dv=dv(sqrt(sum(dv.^2,2))<spin_system.tols.prox_cutoff,:)`.
- Lines 121-122: Detect atomic collisions; implemented by `if any(sqrt(sum(dv.^2,2))<0.5)`.

### Control flow inferred from the code

- Line 40: `for` loop over `m=1:numel(spin_system.chem.parts)`.
- Line 46: `for` loop over `n=spin_list`.
- Line 47: `for` loop over `k=spin_list`.
- Line 50: conditional branch on `(~isempty(spin_system.inter.coordinates{n}))&&`.
- Line 54: conditional branch on `isempty(spin_system.inter.pbc)`.
- Line 67: `for` loop over `p=-spin_system.tols.dd_ncells:spin_system.tols.dd_ncells`.
- Line 81: `for` loop over `p=-spin_system.tols.dd_ncells:spin_system.tols.dd_ncells`.
- Line 82: `for` loop over `q=-spin_system.tols.dd_ncells:spin_system.tols.dd_ncells`.
- Line 98: `for` loop over `p=-spin_system.tols.dd_ncells:spin_system.tols.dd_ncells`.
- Line 99: `for` loop over `q=-spin_system.tols.dd_ncells:spin_system.tols.dd_ncells`.
- Line 100: `for` loop over `r=-spin_system.tols.dd_ncells:spin_system.tols.dd_ncells`.
- Line 122: conditional branch on `any(sqrt(sum(dv.^2,2))<0.5)`.
- Line 131: conditional branch on `~isempty(dv), distvects{n,k}=dv; end`.
- Line 150: `for` loop over `n=1:numel(rows)`.

### Key state/data transformations

- Lines 37: computes `distvects` using `distvects=cell(spin_system.comp.nspins,spin_system.comp.nspins)`.
- Lines 43: computes `spin_list` using `spin_list=spin_system.chem.parts{m}`.
- Lines 57-58: computes `dv` using `dv=spin_system.inter.coordinates{k}- spin_system.inter.coordinates{n}`.
- Lines 63: computes `nvecs` using `nvecs=2*spin_system.tols.dd_ncells+1; dv=zeros(nvecs,3)`.
- Lines 66: computes `linear_index` using `linear_index=1`.
- Lines 68-70: computes `dv(linear_index,:)` using `dv(linear_index,:)=spin_system.inter.coordinates{k}+ p*spin_system.inter.pbc{1}- spin_system.inter.coordinates{n}`.
- Lines 140: computes `spin_system.inter.proxmatrix` using `spin_system.inter.proxmatrix=sparse(~cellfun(@isempty,distvects))`.
- Lines 143: computes `[rows,cols,~]` using `[rows,cols,~]=find(spin_system.inter.proxmatrix)`.
- Lines 156: computes `distance` using `distance=norm(distvects{rows(n),cols(n)}(k,:),2)`.
- Lines 159: computes `ort` using `ort=distvects{rows(n),cols(n)}(k,:)/distance`.
- Lines 162-163: computes `A` using `A=0.5*spin_system.inter.gammas(rows(n))*spin_system.inter.gammas(cols(n))* spin_system.tols.hbar*spin_system.tols.mu0/(4*pi*(distance*1e-10)^3)`.
- Lines 166: computes `D` using `D=A*[1-3*ort(1)*ort(1) -3*ort(1)*ort(2) -3*ort(1)*ort(3)`.
- Lines 181: computes `spin_system.inter.coupling.matrix{rows(n),cols(n)}` using `spin_system.inter.coupling.matrix{rows(n),cols(n)}=D`.

### Local helper functions

- Line 194: `grumble()` — `function grumble(spin_system)`. What worries me about religion is that it teaches people to be satisfied with not understanding.
  - Representative operation: `if ~all(isfield(spin_system,{'comp','inter','chem','tols'}))`.
  - Representative operation: `error('spin_system object is missing essential information.')`.

## Parameters / inputs

- spin_system -Spinach data object containing infor-
- mation about chemical subsystems, ato-
- mic coordinates, and periodic bounda-
- ry conditions

## Outputs

- spin_system -Spinach data object with the interac-
- tion arrays updated with dipolar and
- hyperfine coupling information

## Implementation structure

- Computes dipolar couplings in the presence or absence of periodic
- boundary conditions. This is an auxiliary function of Spinach ker-
- nel, direct calls are discouraged. Use xyz2dd and xyz2hfc to con-
- vert Cartesian coordinates into dipolar and hyperfine couplings
- respectively. Syntax:
- spin_system=dipolar(spin_system)
- spin_system -Spinach data object containing infor-
- mation about chemical subsystems, ato-
- mic coordinates, and periodic bounda-
- ry conditions
- spin_system -Spinach data object with the interac-
- tion arrays updated with dipolar and

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `report()`, `num2str()`, `isscalar()`, `any()`, `cellfun()`, `pair()`, `rows()`, `cols()`, `ort()`, `ismember()`, `all()`, `isfield()`.
