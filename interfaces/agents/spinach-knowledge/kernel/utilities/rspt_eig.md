# kernel/utilities/rspt_eig.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/rspt_eig.m`
- Signature: `[E,V,dE,T,LP]=rspt_eig(spin_system,parameters,Hz,Hc,Hmw,B)`
- Total lines: 196

## Purpose

Eigensystem of sparse Hamiltonians to user-specified order in RSPT with careful handling of diagonal dominance and an opti- on to do exact diagonalisation (expensive). The function also returns eigenvalue derivatives and transition moments between eigenvectors under a user-specified operator. Parametrisation matches use cases in field-swept EPR spectroscopy. Syntax: [E,V,dE,T,LP]=rspt_eig(spin_system,parameters,Hz,Hc

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 60-61: Check consistency; implemented by `grumble(parameters,Hz,Hc,Hmw,B)`.
- Lines 63-64: Recursive call for symmetry; implemented by `if isfield(spin_system.bas,'irrep')`.
- Lines 66-67: Preallocate irrep output blocks; implemented by `n_irreps=numel(spin_system.bas.irrep)`.
- Lines 70-71: Loop over irreps; implemented by `for n=1:n_irreps`.
- Lines 73-74: Extract irrep projector; implemented by `P=spin_system.bas.irrep(n).projector`.
- Lines 76-77: Project Hamiltonian components; implemented by `HzIrr=P'*Hz*P; HzIrr=(HzIrr+HzIrr')/2`.
- Lines 80-81: Project microwave operator; implemented by `HmwIrr=P'*Hmw*P; HmwIrr=(HmwIrr+HmwIrr')/2`.
- Lines 83-84: Issue a recursive call; implemented by `spin_system_nosym=spin_system`.
- Lines 89-90: Project back; implemented by `V{n}=P*V{n}`.
- Lines 94-95: Concatenate irrep blocks; implemented by `E=cell2mat(E); V=cell2mat(V)`.
- Lines 99-100: Hamiltonian; implemented by `H=B*Hz+Hc; H=full((H+H')/2)`.
- Lines 102-103: Decide the method; implemented by `switch parameters.rspt_order`.
- Lines 107-108: Diag and off-diag Hamiltonian; implemented by `H_diag=diag(H); H_offd=H-diag(diag(H))`.
- Lines 110-111: Call perturbation theory; implemented by `[E,V]=rspert(H_diag,H_offd,parameters.rspt_order)`.
- Lines 115-116: Full diagonalisation; implemented by `[V,E]=eig(H,'vector')`.
- Lines 120-121: Complain and bomb out; implemented by `error('unsupported perturbation theory order.')`.
- Lines 127-128: Sort the energies in ascending order; implemented by `[E,idx]=sort(E,'ascend'); V=V(:,idx)`.
- Lines 130-131: Hellmann-Feynman dE/dB, if required; implemented by `if nargout>2, dE=real(diag(V'*Hz*V)); end`.

### Control flow inferred from the code

- Line 64: conditional branch on `isfield(spin_system.bas,'irrep')`.
- Line 71: `for` loop over `n=1:n_irreps`.
- Line 103: dispatches on `parameters.rspt_order`; cases `{1,2,3,4}`, `Inf`.
- Line 131: conditional branch on `nargout>2, dE=real(diag(V'*Hz*V)); end`.
- Line 134: conditional branch on `nargout>3, T=abs(V'*Hmw*V).^2; end`.
- Line 137: conditional branch on `(nargout>4)&&isfield(parameters,'rho0')`.
- Line 140: conditional branch on `isa(parameters.rho0,'function_handle')`.

### Key state/data transformations

- Lines 67: computes `n_irreps` using `n_irreps=numel(spin_system.bas.irrep)`.
- Lines 68: computes `E` using `E=cell(n_irreps,1); V=cell(1,n_irreps)`.
- Lines 74: computes `P` using `P=spin_system.bas.irrep(n).projector`.
- Lines 77: computes `HzIrr` using `HzIrr=P'*Hz*P; HzIrr=(HzIrr+HzIrr')/2`.
- Lines 78: computes `HcIrr` using `HcIrr=P'*Hc*P; HcIrr=(HcIrr+HcIrr')/2`.
- Lines 81: computes `HmwIrr` using `HmwIrr=P'*Hmw*P; HmwIrr=(HmwIrr+HmwIrr')/2`.
- Lines 84: computes `spin_system_nosym` using `spin_system_nosym=spin_system`.
- Lines 85: computes `spin_system_nosym.bas` using `spin_system_nosym.bas=rmfield(spin_system.bas,'irrep')`.
- Lines 86-87: computes `[E{n},V{n}]` using `[E{n},V{n}]=rspt_eig(spin_system_nosym,parameters, HzIrr,HcIrr,HmwIrr,B)`.
- Lines 90: computes `V{n}` using `V{n}=P*V{n}`.
- Lines 100: computes `H` using `H=B*Hz+Hc; H=full((H+H')/2)`.
- Lines 108: computes `H_diag` using `H_diag=diag(H); H_offd=H-diag(diag(H))`.
- Lines 111: computes `[E,V]` using `[E,V]=rspert(H_diag,H_offd,parameters.rspt_order)`.
- Lines 116: computes `[V,E]` using `[V,E]=eig(H,'vector')`.
- Lines 128: computes `[E,idx]` using `[E,idx]=sort(E,'ascend'); V=V(:,idx)`.
- Lines 143-145: computes `rho0` using `rho0=parameters.rho0(B,parameters.orientation(1), parameters.orientation(2), parameters.orientation(3))`.
- Lines 155: computes `LP` using `LP=real(diag(V'*rho0*V))`.

### Local helper functions

- Line 170: `grumble()` — `function grumble(parameters,Hz,Hc,Hmw,B)`.
  - Representative operation: `if ~isfield(parameters,'rspt_order')`.
  - Representative operation: `error('parameters.rspt_order subfield must be present.')`.

## Parameters / inputs

- Hz -laboratory frame Hamiltonian, containing only
- Zeeman terms at 1 Tesla
- Hc -laboratory frame Hamiltonian, containing all
- spin-spin couplings, but no Zeeman terms
- Hmw -observable operator without the amplitude pre-
- factor (2-norm should be around 1)
- B -magnetic field, Tesla
- parameters.rspt_order -perturbation theory order to use
- to account for the off-diagonal
- part of the Hamiltonian, Inf for
- exact diagonalisation
- parameters.rho0 -[optional] when a matrix, sets a user-
- specified thermal equilibrium state;
- when a function handle f(B,alp,bet,gam)
- sets the function to call to obtain
- the thermal equilibrium at each orien-
- tation and magnetic field; if not pro-
- vided, the thermal equilibrium is com-
- puted at the current temperature, ori-
- entation, and magnetic field

## Outputs

- E -a column vector of energies, sorted in ascen-
- ding order (rad/s)
- V -a matrix with eigenvectors in columns, sorted
- left to right in the same order as the energies
- dE -a column vector of dE/dB derivatives, sorted in
- the same order as the energies
- T -a matrix of transition moments under Hmw
- LP -a column vector of energy level populations, sor-
- ted in the same order as the energies

## Implementation structure

- Eigensystem of sparse Hamiltonians to user-specified order in
- RSPT with careful handling of diagonal dominance and an opti-
- on to do exact diagonalisation (expensive). The function also
- returns eigenvalue derivatives and transition moments between
- eigenvectors under a user-specified operator. Parametrisation
- matches use cases in field-swept EPR spectroscopy. Syntax:
- [E,V,dE,T,LP]=rspt_eig(spin_system,parameters,Hz,Hc,Hmw,B)
- Hz - laboratory frame Hamiltonian, containing only
- Zeeman terms at 1 Tesla
- Hc - laboratory frame Hamiltonian, containing all
- spin-spin couplings, but no Zeeman terms
- Hmw - observable operator without the amplitude pre-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isfield()`, `rmfield()`, `cell2mat()`, `rspert()`, `equilibrium()`, `isscalar()`.
