# kernel/states/equilibrium.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/states/equilibrium.m`
- Signature: `rho=equilibrium(spin_system,I,Q,euler_angles)`
- Total lines: 190

## Purpose

Returns the thermal equilibrium state at the current temperature. If the anisotropic part and the orientation parameters are not given, uses the isotropic Hamiltonian, otherwise uses the full Hamiltonian at the speci- fied orientation. Syntax: rho=equilibrium(spin_system,I,Q,euler_angles)

## Physical / mathematical content

- State-construction utilities. These routines build equilibrium states, singlets, triplets, partner-state expansions, and physically meaningful density operators in the active basis.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 49-50: Account for the orientation; implemented by `if nargin==4`.
- Lines 52-53: Check consistency; implemented by `grumble(spin_system,I,Q,euler_angles)`.
- Lines 55-56: Build anisotropic Hamiltonian; implemented by `I=I+orientation(Q,euler_angles)`.
- Lines 60-61: Check consistency; implemented by `grumble(spin_system,I)`.
- Lines 65-66: Build the isotropic Hamiltonian; implemented by `I=hamiltonian(assume(spin_system,'labframe'),'left')`.
- Lines 73-74: Complain and bomb out; implemented by `error('incorrect number of input arguments.')`.
- Lines 78-79: Get the temperature factor; implemented by `beta_factor=spin_system.tols.hbar/(spin_system.tols.kbol*spin_system.rlx.temperature)`.
- Lines 81-82: Decide how to proceed; implemented by `switch spin_system.bas.formalism`.
- Lines 84-85: Liouville space formalisms; implemented by `case {'sphten-liouv','zeeman-liouv'}`.
- Lines 87-88: Get the thermodynamic unit state; implemented by `switch spin_system.bas.formalism`.
- Lines 92-95: Unit population of T(0,0) state, normalisation is such because prod(spin_system.comp.mults) can be- come too large for double precision arithmetic; implemented by `unit=sparse(1,1,1,size(I,2),1)`.
- Lines 99-101: Stretched unit matrix, normalisation matched to the Hilbert space because systems are small; implemented by `unit=speye(prod(spin_system.comp.mults)); unit=unit(:)`.
- Lines 105-106: Catch silly calls; implemented by `if norm(I*unit,1)<1e-10`.
- Lines 110-111: Propagate unit state in imaginary time; implemented by `rho=step(spin_system,I,unit,-1i*beta_factor)`.
- Lines 113-114: Return to CPU if appropriate; implemented by `if isa(rho,'gpuArray'), rho=gather(rho); end`.
- Lines 116-117: Catch failed calls; implemented by `if any(isnan(rho(:)))`.
- Lines 121-122: Divide by partition function; implemented by `rho=rho/dot(unit,rho)`.
- Lines 124-125: Hilbert space; implemented by `case {'zeeman-hilb'}`.

### Control flow inferred from the code

- Line 50: conditional branch on `nargin==4`.
- Line 82: dispatches on `spin_system.bas.formalism`; cases `{'sphten-liouv','zeeman-liouv'}`, `'sphten-liouv'`, `'zeeman-liouv'`, `{'zeeman-hilb'}`.
- Line 88: dispatches on `spin_system.bas.formalism`; cases `'sphten-liouv'`, `'zeeman-liouv'`, `{'zeeman-hilb'}`.
- Line 106: conditional branch on `norm(I*unit,1)<1e-10`.
- Line 114: conditional branch on `isa(rho,'gpuArray'), rho=gather(rho); end`.
- Line 117: conditional branch on `any(isnan(rho(:)))`.
- Line 134: conditional branch on `scaling_factor>1, beta_factor=beta_factor/scaling_factor; end`.
- Line 140: conditional branch on `ismember('gpu',spin_system.sys.enable), rho=gpuArray(rho); end`.
- Line 143: conditional branch on `scaling_factor>1`.
- Line 144: `for` loop over `n=1:n_squarings`.
- Line 151: conditional branch on `isa(rho,'gpuArray'), rho=gather(rho); end`.

### Key state/data transformations

- Lines 56: computes `I` using `I=I+orientation(Q,euler_angles)`.
- Lines 79: computes `beta_factor` using `beta_factor=spin_system.tols.hbar/(spin_system.tols.kbol*spin_system.rlx.temperature)`.
- Lines 95: computes `unit` using `unit=sparse(1,1,1,size(I,2),1)`.
- Lines 111: computes `rho` using `rho=step(spin_system,I,unit,-1i*beta_factor)`.
- Lines 128: computes `mat_norm` using `mat_norm=abs(beta_factor)*cheap_norm(I)`.
- Lines 131: computes `n_squarings` using `n_squarings=max([0 ceil(log2(mat_norm))]); scaling_factor=2^n_squarings`.

### Local helper functions

- Line 158: `grumble()` — `function grumble(spin_system,I,Q,euler_angles)`.
  - Representative operation: `if ~isnumeric(I)`.
  - Representative operation: `error('isotropic Hamiltonian I must be numeric.')`.

## Parameters / inputs

- I -isotropic part of the Hamiltonian left side pro-
- duct superoperator (in Liouville space) or Hamil-
- tonian (in Hilbert space). If this argument is
- omitted, the Hamiltonian is built here and used
- to compute the thermal equilibrium state.
- Q -irreducible components of the anisotropic part
- of the Hamiltonian left side product superopera-
- tor (in Liouville space) or Hamiltonian (in Hil-
- bert space), as returned by hamiltonian.m; this
- is needed when the thermal equilibrium state de-
- pends on the system orientation.
- euler_angles -a row vector of Euler angles (in radians) speci-
- fying the system orientation relative to the in-
- put orientation. If the angles are not supplied,
- only isotropic part of the Hamiltonian is used.

## Outputs

- rho -thermal equilibrium density matrix (Hilbert spa-
- ce) or state vector (Liouville space).
- WARNING: Liouville space calculations must supply left side product su-
- peroperators, not commutation superoperators.
- WARNING: assumptions supplied to the hamiltonian.m call that generates
- I and Q must be 'labframe'.
- WARNING: spin system ground states are commonly degenerate; absolute
- zero temperatures are not supported.

## Implementation structure

- Returns the thermal equilibrium state at the current temperature. If the
- anisotropic part and the orientation parameters are not given, uses the
- isotropic Hamiltonian, otherwise uses the full Hamiltonian at the speci-
- fied orientation. Syntax:
- rho=equilibrium(spin_system,I,Q,euler_angles)
- I -isotropic part of the Hamiltonian left side pro-
- duct superoperator (in Liouville space) or Hamil-
- tonian (in Hilbert space). If this argument is
- omitted, the Hamiltonian is built here and used
- to compute the thermal equilibrium state.
- Q -irreducible components of the anisotropic part
- of the Hamiltonian left side product superopera-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `orientation()`, `hamiltonian()`, `assume()`, `speye()`, `unit()`, `step()`, `gather()`, `any()`, `isnan()`, `rho()`, `dot()`, `cheap_norm()`, `log2()`, `propagator()`, `ismember()`.
