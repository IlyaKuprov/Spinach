# kernel/cache/sle_operators.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/cache/sle_operators.m`
- Signature: `[Lx,Ly,Lz,D,space_basis]=sle_operators(max_rank,int_ranks)`
- Total lines: 347

## Purpose

Wigner D function basis set and rotation generators required by the SLE module. Syntax: [Lx,Ly,Lz,D,space_basis]=sle_operators(max_rank,int_ranks)

## Physical / mathematical content

- Cache-management utilities. These files maintain Spinach temporary or persistent cache state used to avoid repeated expensive construction of large operators or metadata.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `clebsch_gordan_bypass()`, `clebsch_gordan_general()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 48-49: Check consistency; implemented by `grumble(max_rank,int_ranks)`.
- Lines 51-52: Generate cache record name; implemented by `own_path=mfilename('fullpath')`.
- Lines 58-59: Check the cache; implemented by `if exist(cache_file,'file')`.
- Lines 61-62: Lift data from the cache if the file is already available; implemented by `load(cache_file,'space_basis','Lx','Ly','Lz','D')`.
- Lines 66-67: Determine the spatial problem dimension; implemented by `basis_dim=(1+max_rank)*(1+2*max_rank)*(3+2*max_rank)/3`.
- Lines 69-70: Preallocate the spatial basis descriptor; implemented by `space_basis=zeros(basis_dim,3)`.
- Lines 72-73: Populate the spatial basis descriptor; implemented by `for rank=0:max_rank`.
- Lines 75-76: Determine the index extents; implemented by `extents_upper=rank*(2*rank-1)*(2*rank+1)/3+1`.
- Lines 79-80: Assign Wigner function ranks; implemented by `space_basis(extents_upper:extents_lower,1)=rank`.
- Lines 82-83: Assign Wigner function left projection indices; implemented by `space_basis(extents_upper:extents_lower,2)=kron((rank:-1:-rank)',ones(2*rank+1,1))`.
- Lines 85-86: Assign Wigner function right projection indices; implemented by `space_basis(extents_upper:extents_lower,3)=kron(ones(2*rank+1,1),(rank:-1:-rank)')`.
- Lines 90-91: Build spatial L+ generator: pick the states that can be raised; implemented by `source_states=space_basis(space_basis(:,2)<space_basis(:,1),:); new_dim=size(source_states,1)`.
- Lines 93-94: Build spatial L+ generator: find out what they are raised into; implemented by `destin_states=source_states+[zeros(new_dim,1) ones(new_dim,1) zeros(new_dim,1)]`.
- Lines 96-98: Build spatial L+ generator: compute the corresponding matrix elements; implemented by `matrix_elements=sqrt(source_states(:,1).*(source_states(:,1)+1)- source_states(:,2).*(source_states(:,2)+1))`.
- Lines 100-101: Build spatial L+ generator: determine linear indices of source states; implemented by `L=source_states(:,1); M=source_states(:,2); N=source_states(:,3)`.
- Lines 104-105: Build spatial L+ generator: determine linear indices of destination states; implemented by `L=destin_states(:,1); M=destin_states(:,2); N=destin_states(:,3)`.
- Lines 108-109: Build spatial L+ generator: form the L+ matrix; implemented by `Lp=sparse(destin_state_indices,source_state_indices,matrix_elements,basis_dim,basis_dim)`.
- Lines 111-112: Build Lx and Ly rotation generators from L+; implemented by `Lx=(Lp+Lp')/2; Ly=(Lp-Lp')/2i`.

### Control flow inferred from the code

- Line 59: conditional branch on `exist(cache_file,'file')`.
- Line 73: `for` loop over `rank=0:max_rank`.
- Line 125: `for` loop over `int_rank=int_ranks`.
- Line 134: `for` loop over `m=1:n_prj`.
- Line 137: `for` loop over `n=1:n_prj`.
- Line 146: `for` loop over `k=1:basis_dim`.
- Line 174: conditional branch on `int_rank==2`.

### Key state/data transformations

- Lines 52: computes `own_path` using `own_path=mfilename('fullpath')`.
- Lines 54-56: computes `cache_file` using `cache_file=[own_path 'sle_operators_rank_' num2str(max_rank) '_int' sprintf('_%d',int_ranks) '.mat']`.
- Lines 67: computes `basis_dim` using `basis_dim=(1+max_rank)*(1+2*max_rank)*(3+2*max_rank)/3`.
- Lines 70: computes `space_basis` using `space_basis=zeros(basis_dim,3)`.
- Lines 76: computes `extents_upper` using `extents_upper=rank*(2*rank-1)*(2*rank+1)/3+1`.
- Lines 77: computes `extents_lower` using `extents_lower=(1+rank)*(1+2*rank)*(3+2*rank)/3`.
- Lines 80: computes `space_basis(extents_upper:extents_lower,1)` using `space_basis(extents_upper:extents_lower,1)=rank`.
- Lines 83: computes `space_basis(extents_upper:extents_lower,2)` using `space_basis(extents_upper:extents_lower,2)=kron((rank:-1:-rank)',ones(2*rank+1,1))`.
- Lines 86: computes `space_basis(extents_upper:extents_lower,3)` using `space_basis(extents_upper:extents_lower,3)=kron(ones(2*rank+1,1),(rank:-1:-rank)')`.
- Lines 91: computes `source_states` using `source_states=space_basis(space_basis(:,2)<space_basis(:,1),:); new_dim=size(source_states,1)`.
- Lines 94: computes `destin_states` using `destin_states=source_states+[zeros(new_dim,1) ones(new_dim,1) zeros(new_dim,1)]`.
- Lines 97-98: computes `matrix_elements` using `matrix_elements=sqrt(source_states(:,1).*(source_states(:,1)+1)- source_states(:,2).*(source_states(:,2)+1))`.
- Lines 101: computes `L` using `L=source_states(:,1); M=source_states(:,2); N=source_states(:,3)`.
- Lines 102: computes `source_state_indices` using `source_state_indices=L.*(4*L.^2+6*(L-M)+5)/3-M-N+1`.
- Lines 106: computes `destin_state_indices` using `destin_state_indices=L.*(4*L.^2+6*(L-M)+5)/3-M-N+1`.
- Lines 109: computes `Lp` using `Lp=sparse(destin_state_indices,source_state_indices,matrix_elements,basis_dim,basis_dim)`.
- Lines 112: computes `Lx` using `Lx=(Lp+Lp')/2; Ly=(Lp-Lp')/2i`.
- Lines 115: computes `Lz` using `Lz=spdiags(space_basis(:,2),0,basis_dim,basis_dim)`.

### Local helper functions

- Line 208: `clebsch_gordan_bypass()` — `function cg=clebsch_gordan_bypass(L_array,~,~,M1_array,L2_array,M2_array)`. Preallocate the answer
  - Representative operation: `cg=zeros(size(L_array))`.
  - Representative operation: `parfor n=1:numel(L_array)`.
- Line 317: `clebsch_gordan_general()` — `function cg=clebsch_gordan_general(L_array,M_array,L1,M1_array,L2_array,M2_array)`. Preallocate the answer
  - Representative operation: `cg=zeros(size(L_array))`.
  - Representative operation: `parfor n=1:numel(L_array)`.
- Line 331: `grumble()` — `function grumble(max_rank,int_ranks)`.
  - Representative operation: `if (numel(max_rank)~=1)||(~isnumeric(max_rank))||(~isreal(max_rank))|| (max_rank<1)||(mod(max_rank,1)~=0)`.
  - Representative operation: `(max_rank<1)||(mod(max_rank,1)~=0)`.

## Parameters / inputs

- max_rank -maximum L rank for Wigner D functions
- int_ranks -row vector of interaction ranks for which
- product superoperators are required; may
- be empty if only rotation generators are
- needed

## Outputs

- space_basis -lab space basis set descriptor, in
- [L M N] format, giving indices of
- each Wigner function in the basis.
- Lx,Ly,Lz -representations of lab space rotation
- generators in the Wigner function basis,
- to be used in the building of the lab
- space diffusion operator.
- D -a cell array with one element per interac-
- tion rank r in int_ranks, each a cell array
- of Wigner function product superoperators,
- corresponding to multiplication by D[r,M,N]
- of the basis Wigner functions, to be used
- in the building of the spin Hamiltonian
- operator; D{r} has dimensions (2r+1)x(2r+1)
- Automatic caching is implemented -the function would not re-
- compute operator sets that it can find on disk.
- Note: building product superoperators for interaction ranks other
- than 2 calls clebsch_gordan.m, which requires the Java virtual
- machine; cached operator sets load without it.

## Implementation structure

- Wigner D function basis set and rotation generators required by
- the SLE module. Syntax:
- [Lx,Ly,Lz,D,space_basis]=sle_operators(max_rank,int_ranks)
- max_rank -maximum L rank for Wigner D functions
- int_ranks -row vector of interaction ranks for which
- product superoperators are required; may
- be empty if only rotation generators are
- needed
- space_basis - lab space basis set descriptor, in
- [L M N] format, giving indices of
- each Wigner function in the basis.
- Lx,Ly,Lz -representations of lab space rotation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `mfilename()`, `own_path()`, `num2str()`, `exist()`, `load()`, `space_basis()`, `source_states()`, `destin_states()`, `spdiags()`, `clear()`, `destinations()`, `sources()`, `clebsch_gordan_bypass()`, `clebsch_gordan_general()`, `save()`.
