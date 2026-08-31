# experiments/fieldscan_enlev.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/fieldscan_enlev.m`
- Signature: `fieldscan_enlev(spin_system,parameters)`
- Total lines: 113

## Purpose

Plots a user-specified number of the lowest energy levels of the system as a function of the applied magnetic field. The energies are obtained using the Arnoldi method. Syntax: fieldscan_enlev(spin_system,parameters)

## Physical / mathematical content

- This file belongs to the `experiments` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 29-30: Check consistency; implemented by `grumble(spin_system,parameters)`.
- Lines 32-35: Generate the field grid; implemented by `B0=linspace(parameters.fields(1), parameters.fields(2), parameters.npoints)`.
- Lines 37-38: Zeeman Hamiltonian at the current orientation; implemented by `[Hz,Qz]=hamiltonian(assume(spin_system,'labframe','zeeman'))`.
- Lines 41-42: Coupling Hamiltonian at the current orientation; implemented by `[Hc,Qc]=hamiltonian(assume(spin_system,'labframe','couplings'))`.
- Lines 45-46: Loop over the fields; implemented by `parfor n=1:numel(B0)`.
- Lines 48-49: Build the Hamiltonian; implemented by `H=full(B0(n)*Hz+Hc)`.
- Lines 51-52: Get level energies; implemented by `E(:,n)=eigs(H,parameters.nstates,'sr')`.
- Lines 54-55: Sort level energies; implemented by `E(:,n)=sort(real(E(:,n)))`.
- Lines 59-60: Convert the energies to cm^-1; implemented by `E=hz2icm(E/(2*pi))`.
- Lines 62-63: Plot the energy levels; implemented by `kfigure(); plot(B0,E,'r-'); kgrid`.

### Control flow inferred from the code

- Line 46: `parfor` loop over `n=1:numel(B0)`.

### Key state/data transformations

- Lines 33-35: computes `B0` using `B0=linspace(parameters.fields(1), parameters.fields(2), parameters.npoints)`.
- Lines 38: computes `[Hz,Qz]` using `[Hz,Qz]=hamiltonian(assume(spin_system,'labframe','zeeman'))`.
- Lines 39: computes `Hz` using `Hz=Hz+orientation(Qz,parameters.orientation)`.
- Lines 42: computes `[Hc,Qc]` using `[Hc,Qc]=hamiltonian(assume(spin_system,'labframe','couplings'))`.
- Lines 43: computes `Hc` using `Hc=Hc+orientation(Qc,parameters.orientation)`.
- Lines 49: computes `H` using `H=full(B0(n)*Hz+Hc)`.
- Lines 52: computes `E(:,n)` using `E(:,n)=eigs(H,parameters.nstates,'sr')`.
- Lines 60: computes `E` using `E=hz2icm(E/(2*pi))`.

### Local helper functions

- Line 70: `grumble()` — `function grumble(spin_system,parameters)`.
  - Representative operation: `if ~strcmp(spin_system.bas.formalism,'zeeman-hilb')`.
  - Representative operation: `error('this function is only available in zeeman-hilb formalism.')`.

## Parameters / inputs

- parameters.fields -two-element vector in Tesla,
- ordered as [from to]
- parameters.npoints -number of points in the scan
- parameters.orientation -system orientation, three-
- element vector containing
- Euler angles in radians,
- ordered as [alp bet gam]
- parameters.nstates -number of lowest energy states
- to solve for

## Implementation structure

- Plots a user-specified number of the lowest energy levels of
- the system as a function of the applied magnetic field. The
- energies are obtained using the Arnoldi method. Syntax:
- fieldscan_enlev(spin_system,parameters)
- parameters.fields -two-element vector in Tesla,
- ordered as [from to]
- parameters.npoints -number of points in the scan
- parameters.orientation -system orientation, three-
- element vector containing
- Euler angles in radians,
- ordered as [alp bet gam]
- parameters.nstates -number of lowest energy states

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `hamiltonian()`, `assume()`, `orientation()`, `hz2icm()`, `kfigure()`, `kxlabel()`, `kylabel()`, `strcmp()`, `isfield()`, `isrow()`.
