# kernel/correlation.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/correlation.m`
- Signature: `rho=correlation(spin_system,rho,orders,spins)`
- Total lines: 168

## Purpose

Correlation order selection function -keeps only the specified orders of spin correlation in the state vector. This is useful as an analyti- cal replacement for complicated phase cycles. Syntax: rho=correlation(spin_system,rho,correlation_orders,spins)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 39-40: Set the default to all spins; implemented by `if ~exist('spins','var'), spins='all'; end`.
- Lines 42-43: Check consistency; implemented by `grumble(spin_system,rho,orders,spins)`.
- Lines 45-46: Store dimension statistics; implemented by `spn_dim=size(spin_system.bas.basis,1)`.
- Lines 53-54: Fold indirect dimensions; implemented by `rho=reshape(rho,[spn_dim spc_dim])`.
- Lines 56-57: Parse the spin specification; implemented by `if ~isnumeric(spins)`.
- Lines 65-66: Decide how to proceed; implemented by `switch spin_system.bas.formalism`.
- Lines 70-71: Compute the order of correlation for each basis state; implemented by `orders_present=sum(logical(spin_system.bas.basis(:,spins)),2)`.
- Lines 73-74: Wipe all correlation orders except those specified by the user; implemented by `state_mask=false(size(spin_system.bas.basis,1),1)`.
- Lines 79-80: Apply the mask; implemented by `rho(~state_mask,:)=0`.
- Lines 84-85: Hilbert space dimension and multiplicities; implemented by `dim=prod(spin_system.comp.mults)`.
- Lines 88-89: Identity component channels of the selected spins; implemented by `idch=cell(1,numel(spins))`.
- Lines 102-103: Generating function samples on the roots of unity; implemented by `nsel=numel(spins); gsam=cell(1,nsel+1)`.
- Lines 112-113: Discrete Fourier weights of the correlation orders to keep; implemented by `wts=zeros(nsel+1,1); okords=orders(ismember(orders,0:nsel))`.
- Lines 118-119: Assemble the projected state; implemented by `rho=wts(1)*gsam{1}`.
- Lines 126-127: Unfold indirect dimensions; implemented by `rho=reshape(rho,problem_dims)`.
- Lines 129-130: Report overly destructive calls; implemented by `if norm(rho,1)<1e-10`.

### Control flow inferred from the code

- Line 40: conditional branch on `~exist('spins','var'), spins='all'; end`.
- Line 47: conditional branch on `strcmp(spin_system.bas.formalism,'zeeman-hilb')`.
- Line 57: conditional branch on `~isnumeric(spins)`.
- Line 58: conditional branch on `strcmp(spins,'all')`.
- Line 66: dispatches on `spin_system.bas.formalism`; cases `'sphten-liouv'`, `{'zeeman-liouv','zeeman-hilb'}`.
- Line 75: `for` loop over `n=orders`.
- Line 90: `for` loop over `n=1:numel(spins)`.
- Line 92: `for` loop over `p=1:mults(spins(n))`.
- Line 93: `for` loop over `q=1:mults(spins(n))`.
- Line 104: `for` loop over `j=0:nsel`.
- Line 106: `for` loop over `n=1:nsel`.
- Line 114: `for` loop over `j=0:nsel`.
- Line 120: `for` loop over `j=1:nsel`.
- Line 130: conditional branch on `norm(rho,1)<1e-10`.

### Key state/data transformations

- Lines 46: computes `spn_dim` using `spn_dim=size(spin_system.bas.basis,1)`.
- Lines 50: computes `spc_dim` using `spc_dim=numel(rho)/spn_dim`.
- Lines 51: computes `problem_dims` using `problem_dims=size(rho)`.
- Lines 54: computes `rho` using `rho=reshape(rho,[spn_dim spc_dim])`.
- Lines 59: computes `spins` using `spins=1:length(spin_system.comp.isotopes)`.
- Lines 71: computes `orders_present` using `orders_present=sum(logical(spin_system.bas.basis(:,spins)),2)`.
- Lines 74: computes `state_mask` using `state_mask=false(size(spin_system.bas.basis,1),1)`.
- Lines 80: computes `rho(~state_mask,:)` using `rho(~state_mask,:)=0`.
- Lines 85: computes `dim` using `dim=prod(spin_system.comp.mults)`.
- Lines 86: computes `mults` using `mults=spin_system.comp.mults`.
- Lines 89: computes `idch` using `idch=cell(1,numel(spins))`.
- Lines 91: computes `idch{n}` using `idch{n}=sparse(dim^2,dim^2)`.
- Lines 94-96: computes `A` using `A=kron(kron(speye(prod(mults(1:(spins(n)-1)))), sparse(p,q,1,mults(spins(n)),mults(spins(n)))), speye(prod(mults((spins(n)+1):end))))`.
- Lines 103: computes `nsel` using `nsel=numel(spins); gsam=cell(1,nsel+1)`.
- Lines 105: computes `gs` using `gs=rho; x=exp(2i*pi*j/(nsel+1))`.
- Lines 109: computes `gsam{j+1}` using `gsam{j+1}=gs`.
- Lines 113: computes `wts` using `wts=zeros(nsel+1,1); okords=orders(ismember(orders,0:nsel))`.
- Lines 115: computes `wts(j+1)` using `wts(j+1)=sum(exp(-2i*pi*okords(:)*j/(nsel+1)))/(nsel+1)`.

### Local helper functions

- Line 137: `grumble()` — `function grumble(spin_system,rho,correlation_orders,spins)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv','zeeman-hilb'})`.
  - Representative operation: `error('analytical correlation order selection is only available for sphten-liouv, zeeman-liouv, and zeeman-hilb formalisms.')`.

## Parameters / inputs

- rho -a state vector or a horizontal stack thereof;
- in zeeman-hilb, a density matrix or a horizon-
- tal stack thereof
- correlation_orders -a row vector of correlation
- orders to keep
- spins -which spins to consider (e.g.
- '1H', '13C', 'all')

## Outputs

- rho -the state vector with the undesired orders of
- spin correlations zeroed out
- Note: this function requires sphten-liouv, zeeman-liouv, or zeeman-
- hilb formalism; Fokker-Planck direct products are supported in
- the Liouville space formalisms. In the Zeeman formalisms the
- selection is an exact projection built from per-spin identity
- component channels because correlation order is not diagonal
- in the Zeeman basis; zeeman-hilb density matrices are stretch-
- ed into Liouville space, filtered there, and folded back.

## Implementation structure

- Correlation order selection function -keeps only the specified orders
- of spin correlation in the state vector. This is useful as an analyti-
- cal replacement for complicated phase cycles. Syntax:
- rho=correlation(spin_system,rho,correlation_orders,spins)
- rho - a state vector or a horizontal stack thereof;
- in zeeman-hilb, a density matrix or a horizon-
- tal stack thereof
- correlation_orders - a row vector of correlation
- orders to keep
- spins - which spins to consider (e.g.
- '1H', '13C', 'all')
- rho -the state vector with the undesired orders of

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `strcmp()`, `logical()`, `false()`, `rho()`, `mults()`, `spins()`, `speye()`, `conj()`, `orders()`, `ismember()`, `wts()`, `okords()`, `report()`, `vector()`.
