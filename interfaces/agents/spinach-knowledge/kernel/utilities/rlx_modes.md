# kernel/utilities/rlx_modes.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/rlx_modes.m`
- Signature: `R=rlx_modes(spin_system)`
- Total lines: 136

## Purpose

Bosonic mode dissipation superoperator. Builds thermalised GKSL dissipators for the amplitude damping and the pure dephasing of the bosonic modes declared in inter.modes, using the amplitude damping rates and the pure dephasing rates ingested by create.m and the Bose-Einstein thermal occupation numbers computed from the physical mode frequencies, meaning the sum of the declared carrier and the declared frequency wher

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 35-36: Check consistency; implemented by `grumble(spin_system)`.
- Lines 38-39: Preallocate the answer; implemented by `R=mprealloc(spin_system,1)`.
- Lines 41-42: Locate bosonic modes; implemented by `mode_list=find(ismember(spin_system.comp.types,{'C','V','T'}))`.
- Lines 44-45: Loop over bosonic modes; implemented by `for k=mode_list`.
- Lines 47-48: Get the dissipation rates; implemented by `kappa=spin_system.inter.modes.damp(k)`.
- Lines 51-52: Move on if the mode is not dissipative; implemented by `if (kappa==0)&&(gphi==0), continue; end`.
- Lines 54-55: Get the thermal occupation number, only damping needs it; implemented by `if (kappa>0)&&(spin_system.rlx.temperature>0)`.
- Lines 57-58: Physical frequency, laboratory carrier included where declared; implemented by `if spin_system.inter.modes.carriers(k)>0`.
- Lines 65-66: Catch modes whose thermal occupation is undefined; implemented by `if phys_frq<2*pi*spin_system.tols.inter_cutoff`.
- Lines 71-73: Bose-Einstein statistics at the physical frequency; implemented by `beta_factor=spin_system.tols.hbar*phys_frq/ (spin_system.tols.kbol*spin_system.rlx.temperature)`.
- Lines 76-78: Inform the user which frequency has been used; implemented by `report(spin_system,['bosonic mode ' num2str(k) ': thermal occupation from a ' 'physical frequency of ' num2str(phys_frq/(2*pi)) ' Hz, nbar = ' num2str(nbar)])`.
- Lines 82-83: Zero occupation without damping or at zero temperature; implemented by `nbar=0`.
- Lines 87-89: Inform the user; implemented by `report(spin_system,['bosonic mode ' num2str(k) ': kappa = ' num2str(kappa) ' s^-1, gamma_phi = ' num2str(gphi) ' s^-1, nbar = ' num2str(nbar)])`.
- Lines 91-92: Get the ladder superoperators; implemented by `a_left=operator(spin_system,{'A'},{k},'left')`.
- Lines 99-100: Add the cooling dissipator; implemented by `R=R+kappa*(1+nbar)*(a_left*c_right-0.5*(n_left+n_right))`.
- Lines 102-103: Add the heating dissipator; implemented by `if nbar>0`.
- Lines 109-110: Add the pure dephasing dissipator; implemented by `if gphi>0`.

### Control flow inferred from the code

- Line 45: `for` loop over `k=mode_list`.
- Line 52: conditional branch on `(kappa==0)&&(gphi==0), continue; end`.
- Line 55: conditional branch on `(kappa>0)&&(spin_system.rlx.temperature>0)`.
- Line 58: conditional branch on `spin_system.inter.modes.carriers(k)>0`.
- Line 66: conditional branch on `phys_frq<2*pi*spin_system.tols.inter_cutoff`.
- Line 103: conditional branch on `nbar>0`.
- Line 110: conditional branch on `gphi>0`.

### Key state/data transformations

- Lines 39: computes `R` using `R=mprealloc(spin_system,1)`.
- Lines 42: computes `mode_list` using `mode_list=find(ismember(spin_system.comp.types,{'C','V','T'}))`.
- Lines 48: computes `kappa` using `kappa=spin_system.inter.modes.damp(k)`.
- Lines 49: computes `gphi` using `gphi=spin_system.inter.modes.dephase(k)`.
- Lines 59-60: computes `phys_frq` using `phys_frq=spin_system.inter.modes.carriers(k)+ spin_system.inter.modes.frqs(k)`.
- Lines 72-73: computes `beta_factor` using `beta_factor=spin_system.tols.hbar*phys_frq/ (spin_system.tols.kbol*spin_system.rlx.temperature)`.
- Lines 74: computes `nbar` using `nbar=1/(exp(beta_factor)-1)`.
- Lines 78: computes `'physical frequency of ' num2str(phys_frq/(2*pi)) ' Hz, nbar` using `'physical frequency of ' num2str(phys_frq/(2*pi)) ' Hz, nbar = ' num2str(nbar)])`.
- Lines 88-89: computes `report(spin_system,['bosonic mode ' num2str(k) ': kappa` using `report(spin_system,['bosonic mode ' num2str(k) ': kappa = ' num2str(kappa) ' s^-1, gamma_phi = ' num2str(gphi) ' s^-1, nbar = ' num2str(nbar)])`.
- Lines 89: computes `' s^-1, gamma_phi` using `' s^-1, gamma_phi = ' num2str(gphi) ' s^-1, nbar = ' num2str(nbar)])`.
- Lines 92: computes `a_left` using `a_left=operator(spin_system,{'A'},{k},'left')`.
- Lines 93: computes `a_right` using `a_right=operator(spin_system,{'A'},{k},'right')`.
- Lines 94: computes `c_left` using `c_left=operator(spin_system,{'C'},{k},'left')`.
- Lines 95: computes `c_right` using `c_right=operator(spin_system,{'C'},{k},'right')`.
- Lines 96: computes `n_left` using `n_left=operator(spin_system,{'N'},{k},'left')`.
- Lines 97: computes `n_right` using `n_right=operator(spin_system,{'N'},{k},'right')`.
- Lines 104: computes `ac_left` using `ac_left=operator(spin_system,{'AC'},{k},'left')`.
- Lines 105: computes `ac_right` using `ac_right=operator(spin_system,{'AC'},{k},'right')`.

### Local helper functions

- Line 119: `grumble()` — `function grumble(spin_system)`.
  - Representative operation: `if ~isfield(spin_system.inter,'modes')`.
  - Representative operation: `error('bosonic mode information is missing from the spin_system structure.')`.

## Parameters / inputs

- spin_system -Spinach spin system description object
- with bosonic mode information present

## Outputs

- R -bosonic mode dissipation superoperator
- Note: the dissipators are kappa*(1+nbar)*D[a], kappa*nbar*D[c],
- and 2*gamma_phi*D[n], where D[x] is the GKSL dissipator of
- the operator x, built from ladder operators truncated to
- the level count of each mode. The Spinach convention of
- zero temperature meaning the high-temperature limit is not
- applicable to bosonic modes: zero temperature here produces
- zero thermal occupation numbers.

## Implementation structure

- Bosonic mode dissipation superoperator. Builds thermalised GKSL
- dissipators for the amplitude damping and the pure dephasing of
- the bosonic modes declared in inter.modes, using the amplitude
- damping rates and the pure dephasing rates ingested by create.m
- and the Bose-Einstein thermal occupation numbers computed from
- the physical mode frequencies, meaning the sum of the declared
- carrier and the declared frequency where inter.modes.carriers
- is present, and the system temperature. Syntax:
- R=rlx_modes(spin_system)
- spin_system -Spinach spin system description object
- with bosonic mode information present
- R -bosonic mode dissipation superoperator

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `mprealloc()`, `ismember()`, `num2str()`, `report()`, `operator()`, `isfield()`.
