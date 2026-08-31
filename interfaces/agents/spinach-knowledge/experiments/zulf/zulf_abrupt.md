# experiments/zulf/zulf_abrupt.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/zulf/zulf_abrupt.m`
- Signature: `fid=zulf_abrupt(spin_system,parameters,H,R,K)`
- Total lines: 182

## Purpose

Zero-field magnetometry experiment that propagates the initial condition through an exponential drop in the external magnetic field and then runs the detection at zero field. Syntax: fid=zulf_abrupt(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Zero-field experiment implementations. They propagate J-coupled systems in the absence of strong carrier terms and often model abrupt field switching.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 45-46: Check consistency; implemented by `grumble(parameters,H,R,K)`.
- Lines 48-49: Build the coupling Hamiltonian; implemented by `H_z=hamiltonian(assume(spin_system,'labframe','zeeman'))`.
- Lines 52-53: Normalise Zeeman to 1 Tesla; implemented by `H_z=H_z/spin_system.inter.magnet`.
- Lines 55-57: Generate an exponential drop; implemented by `drop=expdrop(spin_system.inter.magnet,parameters.drop_field,parameters.drop_time, parameters.drop_npoints,parameters.drop_rate)`.
- Lines 60-61: Isotropic thermal equilibrium; implemented by `rho=equilibrium(spin_system)`.
- Lines 63-64: Propagate the system through the drop; implemented by `for n=1:numel(drop)`.
- Lines 66-67: Take a propagation step; implemented by `report(spin_system,['magnetic field drop step ' num2str(n) '/' num2str(numel(drop)) ' '])`.
- Lines 72-73: Magnetogyric ratio weights relative to 1H; implemented by `weights=spin_system.inter.gammas/spin('1H')`.
- Lines 75-76: Get gamma-weighted detection state; implemented by `coil=sparse(0)`.
- Lines 81-82: Get gamma-weighted pulse operator; implemented by `Sy=sparse(0)`.
- Lines 88-89: Apply the pulse; implemented by `rho=step(spin_system,Sy,rho,parameters.flip_angle)`.
- Lines 91-92: Compute the digitization parameter; implemented by `timestep=1/parameters.sweep`.
- Lines 94-95: Run the detection at zero field; implemented by `fid=evolution(spin_system,H_c+1i*R+1i*K,coil,rho,timestep,parameters.npoints-1,'observable')`.
- Lines 97-98: Emulate detection hardware; implemented by `switch parameters.detection`.
- Lines 102-104: Do nothing; implemented by `case 'uniaxial'`.
- Lines 106-107: Destroy imaginary part; implemented by `fid=real(fid)`.

### Control flow inferred from the code

- Line 64: `for` loop over `n=1:numel(drop)`.
- Line 77: `for` loop over `n=1:spin_system.comp.nspins`.
- Line 83: `for` loop over `n=1:spin_system.comp.nspins`.
- Line 98: dispatches on `parameters.detection`; cases `'quadrature'`, `'uniaxial'`.

### Key state/data transformations

- Lines 49: computes `H_z` using `H_z=hamiltonian(assume(spin_system,'labframe','zeeman'))`.
- Lines 50: computes `H_c` using `H_c=hamiltonian(assume(spin_system,'labframe','couplings'))`.
- Lines 56-57: computes `drop` using `drop=expdrop(spin_system.inter.magnet,parameters.drop_field,parameters.drop_time, parameters.drop_npoints,parameters.drop_rate)`.
- Lines 58: computes `time_step` using `time_step=parameters.drop_time/parameters.drop_npoints`.
- Lines 61: computes `rho` using `rho=equilibrium(spin_system)`.
- Lines 73: computes `weights` using `weights=spin_system.inter.gammas/spin('1H')`.
- Lines 76: computes `coil` using `coil=sparse(0)`.
- Lines 82: computes `Sy` using `Sy=sparse(0)`.
- Lines 92: computes `timestep` using `timestep=1/parameters.sweep`.
- Lines 95: computes `fid` using `fid=evolution(spin_system,H_c+1i*R+1i*K,coil,rho,timestep,parameters.npoints-1,'observable')`.

### Local helper functions

- Line 114: `grumble()` — `function grumble(parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.

## Parameters / inputs

- .drop_field -the magnetic field that the sample should
- be dropped to, starting from the field spe-
- cified in sys.magnet, Tesla
- .drop_time -drop time, seconds
- .drop_npoints -number of discretisation points in the drop
- .drop_rate -field drop rate, Hz
- .sweep -sweep width during acquisition
- .npoints -number of points during acquisition
- .detection -'uniaxial' to emulate common ZULF
- hardware, 'quadrature' for proper
- frequency sign discrimination
- .flip_angle -pulse flip angle in radians for
- protons; for other nuclei, this
- will be scaled by the gamma ratio

## Outputs

- fid -the free induction decay detected on the
- internally generated gamma-weighted state
- Note: this function ignores the offset parameter and makes its
- own Hamiltonian.

## Implementation structure

- Zero-field magnetometry experiment that propagates the initial condition
- through an exponential drop in the external magnetic field and then runs
- the detection at zero field. Syntax:
- fid=zulf_abrupt(spin_system,parameters,H,R,K)
- .drop_field -the magnetic field that the sample should
- be dropped to, starting from the field spe-
- cified in sys.magnet, Tesla
- .drop_time -drop time, seconds
- .drop_npoints -number of discretisation points in the drop
- .drop_rate -field drop rate, Hz
- .sweep -sweep width during acquisition
- .npoints -number of points during acquisition

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `hamiltonian()`, `assume()`, `expdrop()`, `equilibrium()`, `report()`, `num2str()`, `evolution()`, `drop()`, `spin()`, `weights()`, `state()`, `operator()`, `step()`, `ismatrix()`, `all()`.
