# kernel/optimcon/wrappers/grape_coop.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/wrappers/grape_coop.m`
- Signature: `[traj_data,fidelity,gradient]=grape_coop(phi_profile,spin_system)`
- Total lines: 131

## Purpose

Pairs of cooperative pulses that may be used as components of a phase cycle. The pulses are designed to produce as much of the destination state as they can, and to have imputities of opposite sign. Adding the outcomes of the two experiments then destroys the impurities. Syntax: [traj_data,fidelity,gradient]=grape_coop(phi_profile,spin_system)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 29-30: Check consistency; implemented by `grumble(spin_system)`.
- Lines 32-33: Extract phase profiles; implemented by `n_channels=spin_system.control.ncontrols/2`.
- Lines 37-38: Get target and impurity projectors; implemented by `rho_targ=spin_system.control.rho_targ{1}`.
- Lines 50-51: Make sure final states are available; implemented by `spin_system.control.return_traj=true()`.
- Lines 53-54: Run both experiments; implemented by `[traj_data_a,fidelity_a,gradient_a]=grape_phase(profile_a,spin_system)`.
- Lines 57-58: Project out the impurities; implemented by `dirt_a=cell(numel(traj_data_a),1)`.
- Lines 84-85: Replicate the initial state; implemented by `rho=spin_system.control.rho_init{1}`.
- Lines 89-90: Impurity cancellation gradients; implemented by `spin_system.control.ens_corrs={'rho_ens'}`.
- Lines 96-97: Average fidelity of the two pulses; implemented by `fidelity=(fidelity_a+fidelity_b)/2`.
- Lines 99-100: Penalty on the squared norm of the dirt; implemented by `switch spin_system.bas.formalism`.
- Lines 107-108: Assemble the gradient, impurity term only enters the fidelity slice; implemented by `gradient=cat(1,gradient_a,gradient_b)/2`.
- Lines 111-112: Return both trajectories; implemented by `traj_data={traj_data_a,traj_data_b}`.

### Control flow inferred from the code

- Line 39: dispatches on `spin_system.bas.formalism`; cases `'zeeman-hilb'`.
- Line 42: conditional branch on `abs(targ_norm)==0`.
- Line 59: `for` loop over `n=1:numel(traj_data_a)`.
- Line 60: dispatches on `spin_system.bas.formalism`; cases `'zeeman-hilb'`.
- Line 69: `for` loop over `n=1:numel(traj_data_b)`.
- Line 70: dispatches on `spin_system.bas.formalism`; cases `'zeeman-hilb'`.
- Line 79: `for` loop over `n=1:numel(dirt_sum)`.
- Line 100: dispatches on `spin_system.bas.formalism`; cases `'zeeman-hilb'`.

### Key state/data transformations

- Lines 33: computes `n_channels` using `n_channels=spin_system.control.ncontrols/2`.
- Lines 34: computes `profile_a` using `profile_a=phi_profile(1:n_channels,:)`.
- Lines 35: computes `profile_b` using `profile_b=phi_profile((n_channels+1):end,:)`.
- Lines 38: computes `rho_targ` using `rho_targ=spin_system.control.rho_targ{1}`.
- Lines 41: computes `targ_norm` using `targ_norm=hdot(rho_targ,rho_targ)`.
- Lines 46: computes `P_targ` using `P_targ=rho_targ*rho_targ'`.
- Lines 47: computes `P_dirt` using `P_dirt=eye(size(P_targ))-P_targ`.
- Lines 51: computes `spin_system.control.return_traj` using `spin_system.control.return_traj=true()`.
- Lines 54: computes `[traj_data_a,fidelity_a,gradient_a]` using `[traj_data_a,fidelity_a,gradient_a]=grape_phase(profile_a,spin_system)`.
- Lines 55: computes `[traj_data_b,fidelity_b,gradient_b]` using `[traj_data_b,fidelity_b,gradient_b]=grape_phase(profile_b,spin_system)`.
- Lines 58: computes `dirt_a` using `dirt_a=cell(numel(traj_data_a),1)`.
- Lines 62: computes `rho_a` using `rho_a=traj_data_a{n}.forward{end}`.
- Lines 63: computes `dirt_a{n}` using `dirt_a{n}=rho_a-rho_targ*hdot(rho_targ,rho_a)/targ_norm`.
- Lines 68: computes `dirt_b` using `dirt_b=cell(numel(traj_data_b),1)`.
- Lines 72: computes `rho_b` using `rho_b=traj_data_b{n}.forward{end}`.
- Lines 73: computes `dirt_b{n}` using `dirt_b{n}=rho_b-rho_targ*hdot(rho_targ,rho_b)/targ_norm`.
- Lines 78: computes `dirt_sum` using `dirt_sum=cell(size(dirt_a))`.
- Lines 80: computes `dirt_sum{n}` using `dirt_sum{n}=dirt_a{n}+dirt_b{n}`.

### Local helper functions

- Line 117: `grumble()` — `function grumble(spin_system)`. Morally authoritarian movements are attractive to ugly, miserable, talentless people.
  - Representative operation: `if (numel(spin_system.control.rho_targ)~=1)|| (numel(spin_system.control.rho_init)~=1)`.
  - Representative operation: `(numel(spin_system.control.rho_init)~=1)`.

## Parameters / inputs

- phi_profile -phase profiles of the two pulses,
- concatenated horizontally

## Outputs

- traj_data -trajectory information structure
- fidelity -cooperative fidelity measure
- gradient -cooperative fidelity gradient
- Note: only phase-modulated point-to-point transformations are supported.

## Implementation structure

- Pairs of cooperative pulses that may be used as components of a phase
- cycle. The pulses are designed to produce as much of the destination
- state as they can, and to have imputities of opposite sign. Adding the
- outcomes of the two experiments then destroys the impurities. Syntax:
- [traj_data,fidelity,gradient]=grape_coop(phi_profile,spin_system)
- phi_profile - phase profiles of the two pulses,
- concatenated horizontally
- traj_data - trajectory information structure
- fidelity - cooperative fidelity measure
- gradient - cooperative fidelity gradient
- Note: only phase-modulated point-to-point transformations are supported.
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `phi_profile()`, `hdot()`, `true()`, `grape_phase()`, `ens_catalog()`, `fidelity()`, `cellfun()`, `cat()`, `gradient()`, `gradient_c()`, `gradient_d()`.
