# experiments/relaxan.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/relaxan.m`
- Signature: `[r1,r2,t1,t2,R]=relaxan(spin_system,euler_angles)`
- Total lines: 98

## Purpose

Automated relaxation theory analysis. Prints longitudinal and transverse relaxation rates and times for all spins in the system. Syntax: [r1,r2,t1,t2,R]=relaxan(spin_system,euler_angles)

## Physical / mathematical content

- This file belongs to the `experiments` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 37-38: Move into adjoint representation if needed; implemented by `spin_system=sim2liouv(spin_system,struct(),[],[],[])`.
- Lines 40-41: Check consistency; implemented by `if nargin==2`.
- Lines 47-48: Compute the relaxation superoperator; implemented by `if nargin==2`.
- Lines 54-55: Preallocate outputs; implemented by `r1=zeros(spin_system.comp.nspins,1)`.
- Lines 60-61: Fill the arrays; implemented by `parfor n=1:spin_system.comp.nspins`.
- Lines 68-69: Do the printing; implemented by `report(spin_system,'===============================================================')`.

### Control flow inferred from the code

- Line 41: conditional branch on `nargin==2`.
- Line 48: conditional branch on `nargin==2`.
- Line 61: `parfor` loop over `n=1:spin_system.comp.nspins`.
- Line 72: `for` loop over `n=1:spin_system.comp.nspins`.

### Key state/data transformations

- Lines 38: computes `spin_system` using `spin_system=sim2liouv(spin_system,struct(),[],[],[])`.
- Lines 49: computes `R` using `R=relaxation(spin_system,euler_angles)`.
- Lines 55: computes `r1` using `r1=zeros(spin_system.comp.nspins,1)`.
- Lines 56: computes `r2` using `r2=zeros(spin_system.comp.nspins,1)`.
- Lines 57: computes `t1` using `t1=zeros(spin_system.comp.nspins,1)`.
- Lines 58: computes `t2` using `t2=zeros(spin_system.comp.nspins,1)`.
- Lines 62: computes `Lz` using `Lz=state(spin_system,{'Lz'},{n},'cheap')`.
- Lines 63: computes `r1(n)` using `r1(n)=-real((Lz'*R*Lz)/(Lz'*Lz)); t1(n)=1/r1(n)`.
- Lines 64: computes `Lp` using `Lp=state(spin_system,{'L+'},{n},'cheap')`.
- Lines 65: computes `r2(n)` using `r2(n)=-real((Lp'*R*Lp)/(Lp'*Lp)); t2(n)=1/r2(n)`.

### Local helper functions

- Line 83: `grumble()` — `function grumble(spin_system,euler_angles)`. By the grace of reality and the nature of life, man -every man -is an end in himself, he exists for his own sake, and the achievement of his
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
  - Representative operation: `error('this function is only available in Liouville space.')`.

## Parameters / inputs

- euler_angles -optional euler angles for situations
- when relaxation properties are orien-
- tation-dependent

## Outputs

- r1 -a vector of longitudinal relaxation rates
- for each spin
- r2 -a vector of transverse relaxation rates
- for each spin
- t1 -a vector of longitudinal relaxation times
- for each spin
- t2 -a vector of transverse relaxation times
- for each spin
- R -complete relaxation superoperator
- Note: dynamic frequency shifts are dropped.

## Implementation structure

- Automated relaxation theory analysis. Prints longitudinal
- and transverse relaxation rates and times for all spins in
- the system. Syntax:
- [r1,r2,t1,t2,R]=relaxan(spin_system,euler_angles)
- euler_angles -optional euler angles for situations
- when relaxation properties are orien-
- tation-dependent
- r1 -a vector of longitudinal relaxation rates
- for each spin
- r2 -a vector of transverse relaxation rates
- t1 -a vector of longitudinal relaxation times
- t2 -a vector of transverse relaxation times

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `sim2liouv()`, `grumble()`, `relaxation()`, `state()`, `report()`, `pad()`, `num2str()`, `ismember()`.
