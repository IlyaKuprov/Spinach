# experiments/fieldscan_magn.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/fieldscan_magn.m`
- Signature: `[fields,z_magn]=fieldscan_magn(spin_system,parameters)`
- Total lines: 162

## Purpose

Z magnetization of the sample as a function of magnetic field in a finite-speed magnetic field sweep experiment. Syntax: [fields,z_magn]=fieldscan_magn(spin_system,parameters)

## Physical / mathematical content

- This file belongs to the `experiments` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 39-40: Check consistency; implemented by `grumble(spin_system,parameters)`.
- Lines 42-43: Get the rotation matrix; implemented by `R=euler2dcm(parameters.orientation)`.
- Lines 45-46: Get total magnetic moment operator; implemented by `mz=sparse(0)`.
- Lines 49-50: Get the rotated g-tensor; implemented by `g=R*gtensorof(spin_system,n)*R'`.
- Lines 52-53: Biild elementary operators; implemented by `Sp=operator(spin_system,{'L+'},{n}); Sm=Sp'`.
- Lines 56-57: Build magnetic moment operator; implemented by `mz=mz-g(3,1)*Sx-g(3,2)*Sy-g(3,3)*Sz`.
- Lines 61-64: Generate the field grid; implemented by `fields=linspace(parameters.fields(1), parameters.fields(2), parameters.npoints)`.
- Lines 66-67: Zeeman Hamiltonian; implemented by `[Hz,Qz]=hamiltonian(assume(spin_system,'labframe','zeeman'))`.
- Lines 69-70: Coupling Hamiltonian; implemented by `[Hc,Qc]=hamiltonian(assume(spin_system,'labframe','couplings'))`.
- Lines 72-75: Thermal equilibrium at the starting field; implemented by `rho=equilibrium(spin_system,parameters.fields(1)*Hz+Hc, parameters.fields(1)*Qz+Qc, parameters.orientation)`.
- Lines 78-79: Preallocate the answer; implemented by `z_magn=zeros(size(fields),'like',1i)`.
- Lines 81-82: Hamiltonians at the current orientation; implemented by `Hz=Hz+orientation(Qz,parameters.orientation)`.
- Lines 85-86: Project into the active space; implemented by `if isfield(parameters,'nstates')`.
- Lines 91-92: Time stepping; implemented by `dt=parameters.sweep_time/(parameters.npoints-1)`.
- Lines 94-95: Evolution and acquisition; implemented by `for n=1:parameters.npoints`.
- Lines 97-98: Compute the observable; implemented by `z_magn(n)=hdot(rho,mz)`.
- Lines 100-101: Compute the propagator; implemented by `P=propagator(spin_system,fields(n)*Hz+Hc,dt)`.
- Lines 103-104: Apply the propagator; implemented by `rho=P*rho*P'`.

### Control flow inferred from the code

- Line 47: `for` loop over `n=1:spin_system.comp.nspins`.
- Line 86: conditional branch on `isfield(parameters,'nstates')`.
- Line 95: `for` loop over `n=1:parameters.npoints`.

### Key state/data transformations

- Lines 43: computes `R` using `R=euler2dcm(parameters.orientation)`.
- Lines 46: computes `mz` using `mz=sparse(0)`.
- Lines 50: computes `g` using `g=R*gtensorof(spin_system,n)*R'`.
- Lines 53: computes `Sp` using `Sp=operator(spin_system,{'L+'},{n}); Sm=Sp'`.
- Lines 54: computes `Sx` using `Sx=(Sp+Sm)/2; Sy=(Sp-Sm)/2i; Sz=-1i*(Sx*Sy-Sy*Sx)`.
- Lines 62-64: computes `fields` using `fields=linspace(parameters.fields(1), parameters.fields(2), parameters.npoints)`.
- Lines 67: computes `[Hz,Qz]` using `[Hz,Qz]=hamiltonian(assume(spin_system,'labframe','zeeman'))`.
- Lines 70: computes `[Hc,Qc]` using `[Hc,Qc]=hamiltonian(assume(spin_system,'labframe','couplings'))`.
- Lines 73-75: computes `rho` using `rho=equilibrium(spin_system,parameters.fields(1)*Hz+Hc, parameters.fields(1)*Qz+Qc, parameters.orientation)`.
- Lines 79: computes `z_magn` using `z_magn=zeros(size(fields),'like',1i)`.
- Lines 82: computes `Hz` using `Hz=Hz+orientation(Qz,parameters.orientation)`.
- Lines 83: computes `Hc` using `Hc=Hc+orientation(Qc,parameters.orientation)`.
- Lines 87: computes `[V,~]` using `[V,~]=eigs(mean(fields)*Hz+Hc,parameters.nstates,'sr')`.
- Lines 92: computes `dt` using `dt=parameters.sweep_time/(parameters.npoints-1)`.
- Lines 98: computes `z_magn(n)` using `z_magn(n)=hdot(rho,mz)`.
- Lines 101: computes `P` using `P=propagator(spin_system,fields(n)*Hz+Hc,dt)`.

### Local helper functions

- Line 114: `grumble()` — `function grumble(spin_system,parameters)`.
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
- parameters.sweep_time -sweep time, seconds
- parameters.nstates -(optional) number of lowest energy
- states to use for the effective
- Hamiltonian in the time domain

## Outputs

- fields -magnetic fields in Tesla at each point in time
- z_magn -total sample magnetisation at each point in time
- Note: this function requires Hilbert space formalism.

## Implementation structure

- Z magnetization of the sample as a function of magnetic field in a
- finite-speed magnetic field sweep experiment. Syntax:
- [fields,z_magn]=fieldscan_magn(spin_system,parameters)
- parameters.fields -two-element vector in Tesla,
- ordered as [from to]
- parameters.npoints -number of points in the scan
- parameters.orientation -system orientation, three-
- element vector containing
- Euler angles in radians,
- ordered as [alp bet gam]
- parameters.sweep_time -sweep time, seconds
- parameters.nstates -(optional) number of lowest energy

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `euler2dcm()`, `gtensorof()`, `operator()`, `hamiltonian()`, `assume()`, `equilibrium()`, `orientation()`, `isfield()`, `z_magn()`, `hdot()`, `propagator()`, `fields()`, `strcmp()`, `isrow()`.
