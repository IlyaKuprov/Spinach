# kernel/coherent.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/coherent.m`
- Signature: `rho=coherent(spin_system,mode,alpha)`
- Total lines: 95

## Purpose

Coherent state of a bosonic mode. Builds the normalised trunca- tion of the coherent state with the specified amplitude on the specified bosonic mode, with unit operators on all other parti- cles of the system. Syntax: rho=coherent(spin_system,mode,alpha)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 29-30: Check consistency; implemented by `grumble(spin_system,mode,alpha)`.
- Lines 32-33: Fock state amplitudes up to the truncation level; implemented by `nlev=spin_system.comp.mults(mode)`.
- Lines 36-37: Report the Poisson distribution weight lost to the truncation; implemented by `lost_weight=1-exp(-abs(alpha)^2)*sum(abs(amps).^2)`.
- Lines 41-42: Normalise the truncated state; implemented by `amps=amps/norm(amps,2)`.
- Lines 44-45: Build the mode density matrix; implemented by `rho_mode=sparse(amps.'*conj(amps))`.
- Lines 47-48: Kron up the unit operators of the other particles; implemented by `rho=1`.
- Lines 57-58: Convert into the current formalism; implemented by `switch spin_system.bas.formalism`.
- Lines 62-63: Density matrix is already built; implemented by `rho=full(rho)`.
- Lines 67-68: Stretch the density matrix; implemented by `rho=rho(:)`.
- Lines 72-73: Complain and bomb out; implemented by `error('coherent states are only available in Zeeman formalisms.')`.

### Control flow inferred from the code

- Line 49: `for` loop over `n=1:spin_system.comp.nspins`.
- Line 50: conditional branch on `n==mode`.
- Line 58: dispatches on `spin_system.bas.formalism`; cases `'zeeman-hilb'`, `'zeeman-liouv'`.

### Key state/data transformations

- Lines 33: computes `nlev` using `nlev=spin_system.comp.mults(mode)`.
- Lines 34: computes `amps` using `amps=alpha.^(0:(nlev-1))./sqrt(factorial(0:(nlev-1)))`.
- Lines 37: computes `lost_weight` using `lost_weight=1-exp(-abs(alpha)^2)*sum(abs(amps).^2)`.
- Lines 45: computes `rho_mode` using `rho_mode=sparse(amps.'*conj(amps))`.
- Lines 48: computes `rho` using `rho=1`.

### Local helper functions

- Line 80: `grumble()` — `function grumble(spin_system,mode,alpha)`.
  - Representative operation: `if ~isfield(spin_system,'bas')`.
  - Representative operation: `error('basis set information is missing, run basis() before calling this function.')`.

## Parameters / inputs

- mode -index of a bosonic mode in sys.isotopes
- alpha -coherent state amplitude, a complex scalar

## Outputs

- rho -coherent state density matrix (zeeman-hilb)
- or its vectorisation (zeeman-liouv)
- Note: the Fock space truncation of the mode chops the tail of
- the Poisson distribution; the state is renormalised after
- the truncation and the lost weight is reported.

## Implementation structure

- Coherent state of a bosonic mode. Builds the normalised trunca-
- tion of the coherent state with the specified amplitude on the
- specified bosonic mode, with unit operators on all other parti-
- cles of the system. Syntax:
- rho=coherent(spin_system,mode,alpha)
- mode -index of a bosonic mode in sys.isotopes
- alpha -coherent state amplitude, a complex scalar
- rho -coherent state density matrix (zeeman-hilb)
- or its vectorisation (zeeman-liouv)
- Note: the Fock space truncation of the mode chops the tail of
- the Poisson distribution; the state is renormalised after
- the truncation and the lost weight is reported.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `factorial()`, `report()`, `num2str()`, `conj()`, `speye()`, `rho()`, `isfield()`, `basis()`, `isscalar()`, `ismember()`.
