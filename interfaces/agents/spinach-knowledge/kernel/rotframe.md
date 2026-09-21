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
