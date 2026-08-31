# kernel/trajsimil.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/trajsimil.m`
- Signature: `score=trajsimil(spin_system,trajectory_1,trajectory_2,scorefcn)`
- Total lines: 183

## Purpose

Computes trajectory similarity scores. Returns a function representing "similarity" of the two state space trajectories at different points in time. See http://dx.doi.org/10.1016/j.jmr.2013.02.012 for further infor- mation. Syntax: score=trajsimil(spin_system,trajectory_1,trajectory_2,scorefcn)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 59-60: Check consistency; implemented by `grumble(spin_system,trajectory_1,trajectory_2,scorefcn)`.
- Lines 62-63: Run state grouping; implemented by `if strncmp(scorefcn,'SG-',3)||strncmp(scorefcn,'BSG',3)`.
- Lines 65-66: Grab the description of the current basis; implemented by `state_list=spin_system.bas.basis`.
- Lines 68-69: Tell the user we're started; implemented by `report(spin_system,'collapsing equivalent subspaces ')`.
- Lines 71-72: Decide how to proceed; implemented by `if strncmp(scorefcn,'SG-',3)`.
- Lines 74-75: Rename all T(l,-m) states into T(l,m) states; implemented by `[L,M]=lin2lm(state_list); state_list=lm2lin(L,abs(M))`.
- Lines 77-78: Update the method variable; implemented by `scorefcn=scorefcn(4:end)`.
- Lines 82-83: Rename all non-identity states into Lz; implemented by `state_list(state_list~=0)=2`.
- Lines 85-86: Update the method variable; implemented by `scorefcn=scorefcn(5:end)`.
- Lines 90-91: Index all unique and repeated states on the modified state list; implemented by `[grouped_state_list,index_forward,index_backward]=unique(state_list,'rows')`.
- Lines 93-94: Preallocate state-grouped trajectories; implemented by `grouped_trajectory_1=zeros(length(index_forward),size(trajectory_1,2))`.
- Lines 97-99: Group trajectory tracks corresponding to the states that are flagged as identical in the indices (root-sum-square); implemented by `for n=1:length(index_forward)`.
- Lines 104-105: Update the variables; implemented by `trajectory_1=grouped_trajectory_1`.
- Lines 108-109: Tell the user we're done; implemented by `report(spin_system,[num2str(size(state_list,1)) ' states collected into ' num2str(size(grouped_state_list,1)) ' groups.'])`.
- Lines 113-114: Preallocate the result array; implemented by `trajectory_length=size(trajectory_1,2)`.
- Lines 117-118: Compute the score function; implemented by `switch scorefcn`.
- Lines 122-123: Compute running scalar product scores; implemented by `for n=1:trajectory_length`.
- Lines 129-130: Compute running difference norm scores; implemented by `for n=1:trajectory_length`.

### Control flow inferred from the code

- Line 63: conditional branch on `strncmp(scorefcn,'SG-',3)||strncmp(scorefcn,'BSG',3)`.
- Line 72: conditional branch on `strncmp(scorefcn,'SG-',3)`.
- Line 99: `for` loop over `n=1:length(index_forward)`.
- Line 118: dispatches on `scorefcn`; cases `'RSP'`, `'RDN'`.
- Line 123: `for` loop over `n=1:trajectory_length`.
- Line 130: `for` loop over `n=1:trajectory_length`.

### Key state/data transformations

- Lines 66: computes `state_list` using `state_list=spin_system.bas.basis`.
- Lines 75: computes `[L,M]` using `[L,M]=lin2lm(state_list); state_list=lm2lin(L,abs(M))`.
- Lines 78: computes `scorefcn` using `scorefcn=scorefcn(4:end)`.
- Lines 91: computes `[grouped_state_list,index_forward,index_backward]` using `[grouped_state_list,index_forward,index_backward]=unique(state_list,'rows')`.
- Lines 94: computes `grouped_trajectory_1` using `grouped_trajectory_1=zeros(length(index_forward),size(trajectory_1,2))`.
- Lines 95: computes `grouped_trajectory_2` using `grouped_trajectory_2=zeros(length(index_forward),size(trajectory_2,2))`.
- Lines 105: computes `trajectory_1` using `trajectory_1=grouped_trajectory_1`.
- Lines 106: computes `trajectory_2` using `trajectory_2=grouped_trajectory_2`.
- Lines 114: computes `trajectory_length` using `trajectory_length=size(trajectory_1,2)`.
- Lines 115: computes `score` using `score=zeros(1,trajectory_length)`.
- Lines 124: computes `score(n)` using `score(n)=dot(trajectory_1(:,n),trajectory_2(:,n))`.

### Local helper functions

- Line 144: `grumble()` — `function grumble(spin_system,trajectory_1,trajectory_2,scorefcn)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
  - Representative operation: `error('this function is only available for Liouville space formalisms.')`.

## Parameters / inputs

- trajectory_1,2 -spin system trajectories, supplied as nstates
- x nsteps matrices.
- scorefcn -similarity scoring method; possibilities are:
- 'RSP' -running scalar product. Computes
- scalar products between the cor-
- responding vectors of the trajec-
- tories.
- 'RDN' -running difference norm. The two
- trajectories are subtracted and
- difference 2-norms returned.
- 'SG-' -prefix that turns on state grou-
- ping. T(l,m) and T(l,-m) states
- of each spin (standalone or in
- direct products with other ope-
- rators)will be considered equva-
- lent.
- 'BSG-' -prefix that turns on broad state
- grouping. All states of a given
- spin (standalone or in direct
- products with other operators)
- will be considered equivalent.
- The possible combinations are: 'RSP','RDN',
- 'SG-RSP','SG-RDN','BSG-RSP','BSG-RDN'.
- State grouping consists in summing the absolute squares of the coeffi-
- cients to be grouped and taking the square root. The trajectories would
- usually come out of the evolution.m or krylov.m run from a given star-
- ting point under a given Liouvillian.
- Output:
- score -the similarity score vector, one element
- per time slice
- Note: SG and BSG options require sphten-liouv formalism.

## Implementation structure

- Computes trajectory similarity scores. Returns a function representing
- "similarity" of the two state space trajectories at different points in
- time. See http://dx.doi.org/10.1016/j.jmr.2013.02.012 for further infor-
- mation. Syntax:
- score=trajsimil(spin_system,trajectory_1,trajectory_2,scorefcn)
- trajectory_1,2 -spin system trajectories, supplied as nstates
- x nsteps matrices.
- scorefcn -similarity scoring method; possibilities are:
- 'RSP' -running scalar product. Computes
- scalar products between the cor-
- responding vectors of the trajec-
- tories.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `strncmp()`, `report()`, `lin2lm()`, `lm2lin()`, `scorefcn()`, `state_list()`, `grouped_trajectory_1()`, `trajectory_1()`, `grouped_trajectory_2()`, `trajectory_2()`, `num2str()`, `score()`, `dot()`, `ismember()`, `ischar()`.
