# kernel/rotframe.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/rotframe.m`
- Signature: `Hr=rotframe(spin_system,H0,H,isotope,order)`
- Total lines: 98

## Purpose

Rotating frame transformation with respect to specified spins to specified order in perturbation theory, using the formalism described in https://doi.org/10.1063/1.4928978 Syntax: Hr=rotframe(spin_system,H0,H,isotope,order)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 35-36: Check consistency; implemented by `grumble(spin_system,H0,H,isotope,order)`.
- Lines 38-39: Compute the period; implemented by `switch spin_system.bas.formalism`.
- Lines 43-44: Liouville space period for H0; implemented by `T=-2*pi/(spin(isotope)*spin_system.inter.magnet)`.
- Lines 48-49: Hilbert space period for H0; implemented by `T=-4*pi/(spin(isotope)*spin_system.inter.magnet)`.
- Lines 53-54: Complain and bomb out; implemented by `error('unknown formalism specification.')`.
- Lines 58-59: Run the interaction representation transformation; implemented by `Hr=intrep(spin_system,H0,H,T,order)`.

### Control flow inferred from the code

- Line 39: dispatches on `spin_system.bas.formalism`; cases `{'zeeman-liouv','sphten-liouv'}`, `{'zeeman-hilb','zeeman-wavef'}`.

### Key state/data transformations

- Lines 44: computes `T` using `T=-2*pi/(spin(isotope)*spin_system.inter.magnet)`.
- Lines 59: computes `Hr` using `Hr=intrep(spin_system,H0,H,T,order)`.

### Local helper functions

- Line 64: `grumble()` — `function grumble(spin_system,H0,H,isotope,order)`.
  - Representative operation: `if ~ischar(isotope)`.
  - Representative operation: `error('isotope must be a character string')`.

## Parameters / inputs

- H0 -carrier Hamiltonian with respect to which the
- rotating frame transformation is to be done
- H -laboratory frame Hamiltonian H0+H1 that is to
- be transformed into the rotating frame
- isotope -string, such as '1H', specifying the spins
- with respect to which the transformation is
- being computed
- order -perturbation theory order in the rotating
- frame transformation, this may be inf

## Outputs

- Hr -rotating frame Hamiltonian
- Notes: the auxiliary matrix method is massively faster than
- either commutator series or diagonalisation.

## Implementation structure

- Rotating frame transformation with respect to specified spins
- to specified order in perturbation theory, using the formalism
- described in https://doi.org/10.1063/1.4928978 Syntax:
- Hr=rotframe(spin_system,H0,H,isotope,order)
- H0 -carrier Hamiltonian with respect to which the
- rotating frame transformation is to be done
- H -laboratory frame Hamiltonian H0+H1 that is to
- be transformed into the rotating frame
- isotope -string, such as '1H', specifying the spins
- with respect to which the transformation is
- being computed
- order -perturbation theory order in the rotating

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `intrep()`, `ischar()`, `isfield()`, `assume()`, `ismember()`, `isotope()`, `ishermitian()`, `isinf()`.
