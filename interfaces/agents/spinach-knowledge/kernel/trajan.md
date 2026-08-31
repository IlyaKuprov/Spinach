# kernel/trajan.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/trajan.m`
- Signature: `trajan(spin_system,traj,property,time_axis)`
- Total lines: 347

## Purpose

Trajectory analysis function. Plots the time dependence of the densi- ty matrix norm, partitioned into user-specified property classes. See for further information. Syntax: trajan(spin_system,traj,property,time_axis)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 68-69: Check the input; implemented by `if ~exist('time_axis','var'), time_axis=[]; end`.
- Lines 72-73: Project out unit state; implemented by `if ~strcmp(property,'level_populations')`.
- Lines 78-79: Determine how to proceed; implemented by `switch property`.
- Lines 83-84: Determine the correlation order of each state; implemented by `correlation_orders=sum(logical(spin_system.bas.basis),2)`.
- Lines 86-87: Find out which correlation orders are present; implemented by `unique_correlation_orders=unique(correlation_orders)`.
- Lines 89-90: Eliminate zero-spin order; implemented by `unique_correlation_orders=setdiff(unique_correlation_orders,0)`.
- Lines 92-93: Preallocate the norm trajectory array; implemented by `result=zeros(numel(unique_correlation_orders),size(traj,2))`.
- Lines 95-96: Loop over the unique correlation orders that are present; implemented by `for n=1:numel(unique_correlation_orders)`.
- Lines 98-99: Find the subspace of the given correlation order; implemented by `subspace_mask=(correlation_orders==unique_correlation_orders(n))`.
- Lines 101-102: Get the part of the trajectory belonging to the subspace; implemented by `subspace_trajectory=traj(subspace_mask,:)`.
- Lines 104-105: Get the norm of the trajectory; implemented by `result(n,:)=sqrt(sum(subspace_trajectory.*conj(subspace_trajectory),1))`.
- Lines 109-110: Create a legend; implemented by `legend_text=cell(numel(unique_correlation_orders),1)`.
- Lines 115-116: Create labels; implemented by `label_text='correlation order amplitude'`.
- Lines 119-120: Compute Y axis extents; implemented by `max_val=max(result,[],'all')`.
- Lines 125-126: Determine projection quantum numbers of the basis; implemented by `[~,M]=lin2lm(spin_system.bas.basis)`.
- Lines 128-129: Determine the coherence order of each state; implemented by `coherence_orders=sum(M,2)`.
- Lines 131-132: Find out which coherence orders are present; implemented by `unique_coherence_orders=unique(coherence_orders)`.
- Lines 134-135: Preallocate the norm trajectory array; implemented by `result=zeros(numel(unique_coherence_orders),size(traj,2))`.

### Control flow inferred from the code

- Line 69: conditional branch on `~exist('time_axis','var'), time_axis=[]; end`.
- Line 73: conditional branch on `~strcmp(property,'level_populations')`.
- Line 79: dispatches on `property`; cases `'correlation_order'`, `'coherence_order'`, `'total_each_spin'`, `'local_each_spin'`, `'level_populations'`.
- Line 96: `for` loop over `n=1:numel(unique_correlation_orders)`.
- Line 111: `for` loop over `n=1:numel(unique_correlation_orders)`.
- Line 138: `for` loop over `n=1:numel(unique_coherence_orders)`.
- Line 153: `for` loop over `n=1:numel(unique_coherence_orders)`.
- Line 171: `for` loop over `n=1:spin_system.comp.nspins`.
- Line 186: `for` loop over `n=1:spin_system.comp.nspins`.
- Line 187: conditional branch on `~isempty(spin_system.comp.labels{n})`.
- Line 214: `for` loop over `n=1:spin_system.comp.nspins`.
- Line 230: `for` loop over `n=1:spin_system.comp.nspins`.
- Line 231: conditional branch on `~isempty(spin_system.comp.labels{n})`.
- Line 265: `for` loop over `n=1:size(traj,2)`.

### Key state/data transformations

- Lines 74: computes `unit` using `unit=unit_state(spin_system)`.
- Lines 75: computes `traj` using `traj=traj-(unit*unit')*traj`.
- Lines 84: computes `correlation_orders` using `correlation_orders=sum(logical(spin_system.bas.basis),2)`.
- Lines 87: computes `unique_correlation_orders` using `unique_correlation_orders=unique(correlation_orders)`.
- Lines 93: computes `result` using `result=zeros(numel(unique_correlation_orders),size(traj,2))`.
- Lines 102: computes `subspace_trajectory` using `subspace_trajectory=traj(subspace_mask,:)`.
- Lines 105: computes `result(n,:)` using `result(n,:)=sqrt(sum(subspace_trajectory.*conj(subspace_trajectory),1))`.
- Lines 110: computes `legend_text` using `legend_text=cell(numel(unique_correlation_orders),1)`.
- Lines 112: computes `legend_text{n}` using `legend_text{n}=[num2str(unique_correlation_orders(n)) '-spin']`.
- Lines 116: computes `label_text` using `label_text='correlation order amplitude'`.
- Lines 117: computes `title_text` using `title_text='correlation orders'`.
- Lines 120: computes `max_val` using `max_val=max(result,[],'all')`.
- Lines 121: computes `y_axis_extents` using `y_axis_extents=[-0.05*max_val 1.05*max_val]`.
- Lines 126: computes `[~,M]` using `[~,M]=lin2lm(spin_system.bas.basis)`.
- Lines 129: computes `coherence_orders` using `coherence_orders=sum(M,2)`.
- Lines 132: computes `unique_coherence_orders` using `unique_coherence_orders=unique(coherence_orders)`.
- Lines 259: computes `nlevels` using `nlevels=sqrt(size(traj,1))`.
- Lines 266: computes `result(:,n)` using `result(:,n)=real(diag(reshape(traj(:,n),[nlevels nlevels])))`.

### Local helper functions

- Line 319: `grumble()` — `function grumble(spin_system,trajectory,property,time_axis)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('trajectory analysis is only available for sphten-liouv formalism.')`.

## Parameters / inputs

- traj -a stack of state vectors of any length. The
- number of rows in the trajectory array must
- match the number of states in the basis.
- property -if set to 'correlation_order', returns the
- time dependence of the total populations of
- one-spin, two-spin, three-spin, etc. corre-
- lations in the system.
- if set to 'coherence_order', returns the ti-
- me dependence of different orders of coheren-
- ce in the system, where a coherence order is
- defined as the sum of projection quantum num-
- bers in the spherical tensor representation
- of each state.
- if set to 'total_each_spin', returns the ti-
- me dependence of total state space populati-
- on that involves each individual spin in the
- system in any way (all local populations and
- coherences of the spin as well as all of its
- correlations to all third party spins).
- if set to 'local_each_spin', returns the ti-
- me dependence of the population of the sub-
- space of states that are local to each indi-
- vidual spin and do not involve any correla-
- tions to other spins in the system.
- if set to 'level_populations', returns the
- populations of the Zeeman energy levels.
- time_axis -(optional) user specified time axis, a row
- vector of time positions of each state vec-
- tor inthe trajectory array.
- The trajectory would usually come out of the evolution.m run from a
- given starting point under a given Liouvillian.
- Output:
- this function writes into the current figure
- Note: this function is only applicable to the trajectories recorded
- in sphten-liouv formalism.
- Note: unit state population is ignored.

## Implementation structure

- Trajectory analysis function. Plots the time dependence of the densi-
- ty matrix norm, partitioned into user-specified property classes. See
- for further information. Syntax:
- trajan(spin_system,traj,property,time_axis)
- traj -a stack of state vectors of any length. The
- number of rows in the trajectory array must
- match the number of states in the basis.
- property -if set to 'correlation_order', returns the
- time dependence of the total populations of
- one-spin, two-spin, three-spin, etc. corre-
- lations in the system.
- if set to 'coherence_order', returns the ti-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `strcmp()`, `unit_state()`, `logical()`, `setdiff()`, `unique_correlation_orders()`, `traj()`, `result()`, `conj()`, `num2str()`, `lin2lm()`, `unique_coherence_orders()`, `sphten2zeeman()`, `kxlabel()`, `hsv2rgb()`.
