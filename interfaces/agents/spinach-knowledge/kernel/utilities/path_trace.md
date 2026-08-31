# kernel/utilities/path_trace.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/path_trace.m`
- Signature: `projectors=path_trace(spin_system,L,rho)`
- Total lines: 201

## Purpose

Liouvillian path tracing. Treats the user-supplied Liouvillian as the adjacency matrix of a graph, computes the weakly connect- ed subgraphs of that graph and returns a cell array of project- ors into independently evolving populated subspaces. Syntax: projectors=path_trace(spin_system,L,rho)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 38-39: Check the input; implemented by `grumble(L,rho)`.
- Lines 41-42: Check run conditions; implemented by `if ismember('pt',spin_system.sys.disable)`.
- Lines 44-45: Return a unit projector if path tracing is disabled; implemented by `report(spin_system,'WARNING - path tracing disabled by the user.')`.
- Lines 50-51: Return a unit projector if the space is small anyway; implemented by `report(spin_system,'small state space - path tracing skipped.')`.
- Lines 56-57: Report to the user; implemented by `report(spin_system,['analyzing ' num2str(size(L,1)) '-dimensional state space.'])`.
- Lines 63-64: Get the connectivity matrix; implemented by `G=(abs(L)>spin_system.tols.liouv_zero)`.
- Lines 66-67: Make sure isolated states do not get lost; implemented by `G=or(G,transpose(G)); G=or(G,speye(size(G)))`.
- Lines 69-70: Get the weakly connected subgraphs; implemented by `member_states=scomponents(G)`.
- Lines 72-73: Determine the number of subspaces; implemented by `n_subspaces=max(member_states)`.
- Lines 77-78: All subspaces matter by default; implemented by `subspace_important=true(n_subspaces,1)`.
- Lines 80-81: Screen by population; implemented by `if ~isempty(rho)`.
- Lines 83-84: Extract switches for parfor loop; implemented by `tolerance=spin_system.tols.subs_drop`.
- Lines 87-88: Look at actual populations; implemented by `parfor n=1:n_subspaces`.
- Lines 90-91: Determine the importance; implemented by `if ismember(formalism,{'sphten-liouv','zeeman-liouv'})`.
- Lines 93-94: Liouville space has state vectors; implemented by `subspace_important(n)=(norm(rho.*(member_states==n),1)>tolerance)`.
- Lines 98-100: Hilbert space has state matrices; implemented by `subspace_important(n)=(norm(rho(member_states==n,:),1)>tolerance)| (norm(rho(:,member_states==n),1)>tolerance)`.
- Lines 104-105: Complain and bomb out; implemented by `error('unexpected formalism specification.')`.
- Lines 113-114: Ignore unpopulated subspaces; implemented by `significant_subspaces=find(subspace_important)`.

### Control flow inferred from the code

- Line 42: conditional branch on `ismember('pt',spin_system.sys.disable)`.
- Line 81: conditional branch on `~isempty(rho)`.
- Line 88: `parfor` loop over `n=1:n_subspaces`.
- Line 91: conditional branch on `ismember(formalism,{'sphten-liouv','zeeman-liouv'})`.
- Line 121: `parfor` loop over `n=1:n_subspaces`.
- Line 139: `for` loop over `n=1:numel(unique_dims)`.
- Line 149: conditional branch on `ismember('merge',spin_system.sys.disable)`.
- Line 165: `for` loop over `n=1:numel(bins)`.
- Line 175: `for` loop over `n=1:numel(dims)`.

### Key state/data transformations

- Lines 46: computes `projectors` using `projectors={1}; return`.
- Lines 64: computes `G` using `G=(abs(L)>spin_system.tols.liouv_zero)`.
- Lines 70: computes `member_states` using `member_states=scomponents(G)`.
- Lines 73: computes `n_subspaces` using `n_subspaces=max(member_states)`.
- Lines 78: computes `subspace_important` using `subspace_important=true(n_subspaces,1)`.
- Lines 84: computes `tolerance` using `tolerance=spin_system.tols.subs_drop`.
- Lines 85: computes `formalism` using `formalism=spin_system.bas.formalism`.
- Lines 114: computes `significant_subspaces` using `significant_subspaces=find(subspace_important)`.
- Lines 127: computes `subspace_dim` using `subspace_dim=numel(state_index)`.
- Lines 130: computes `projectors{n}` using `projectors{n}=sparse(state_index,1:subspace_dim,ones(1,subspace_dim),size(L,1),subspace_dim)`.
- Lines 135: computes `subspace_dims` using `subspace_dims=cellfun(@(x)size(x,2),projectors)`.
- Lines 136: computes `unique_dims` using `unique_dims=unique(subspace_dims)`.
- Lines 161: computes `bins` using `bins=binpack(subspace_dims,spin_system.tols.merge_dim)`.
- Lines 164: computes `new_projectors` using `new_projectors=cell(numel(bins),1)`.
- Lines 166: computes `new_projectors{n}` using `new_projectors{n}=[projectors{bins{n}}]`.
- Lines 171: computes `dims` using `dims=cellfun(@(x)size(x,2),projectors)`.

### Local helper functions

- Line 185: `grumble()` — `function grumble(L,rho)`.
  - Representative operation: `if (~isnumeric(L))||(~isnumeric(rho))`.
  - Representative operation: `error('both inputs must be numeric.')`.

## Parameters / inputs

- L -Hamiltonian or Liouvillian matrix
- rho -the initial state (source state screening)
- or the detection state (destination state
- screening); pass [] to disable screening

## Outputs

- projectors -a cell array of projectors into independently
- evolving populated subspaces. The projectors
- are to be used as follows:
- L_reduced=P'*L*P; (for matrices)
- rho_reduced=P'*rho; (for state vectors)
- Note: further information on how this function works is availa-
- ble in our papers on this subject

## Implementation structure

- Liouvillian path tracing. Treats the user-supplied Liouvillian
- as the adjacency matrix of a graph, computes the weakly connect-
- ed subgraphs of that graph and returns a cell array of project-
- ors into independently evolving populated subspaces. Syntax:
- projectors=path_trace(spin_system,L,rho)
- L - Hamiltonian or Liouvillian matrix
- rho - the initial state (source state screening)
- or the detection state (destination state
- screening); pass [] to disable screening
- projectors -a cell array of projectors into independently
- evolving populated subspaces. The projectors
- are to be used as follows:

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ismember()`, `report()`, `num2str()`, `transpose()`, `speye()`, `scomponents()`, `true()`, `subspace_important()`, `rho()`, `significant_subspaces()`, `cellfun()`, `nnz()`, `unique_dims()`, `binpack()`, `int2str()`.
