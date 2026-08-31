# kernel/states/deut_pair.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/states/deut_pair.m`
- Signature: `[S,T,Q,Tc,Qc]=deut_pair(spin_system,spin_a,spin_b,options)`
- Total lines: 198

## Purpose

All possible states of a spin-1 pair, classified by the total spin into singlet, triplet, and quartet. Syntax: [S,T,Q,Tc,Qc]=deut_pair(spin_system,spin_a,spin_b,options)

## Physical / mathematical content

- State-construction utilities. These routines build equilibrium states, singlets, triplets, partner-state expansions, and physically meaningful density operators in the active basis.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 49-50: Set defaults; implemented by `if ~exist('options','var'), options.dephasing=0; end`.
- Lines 52-53: Check consistency; implemented by `grumble(spin_system,spin_a,spin_b,options)`.
- Lines 55-57: Build component vectors in Hilbert space as per Eq 1; implemented by `alp=[1; 0; 0]; bet=[0; 1; 0]; gam=[0; 0; 1]`.
- Lines 68-69: Obtain spherical tensor expansions; implemented by `S=0; T={0,0,0}; Q={0,0,0,0,0}; Tc={0,0,0,0}`.
- Lines 74-75: Get spherical tensor indices; implemented by `[L1,M1]=lin2lm(n-1); [L2,M2]=lin2lm(k-1)`.
- Lines 77-78: Skip non-stationary states on user request; implemented by `if ((M1~=0)||(M2~=0))&&(options.dephasing==1), continue; end`.
- Lines 80-81: Build Spinach state descriptor; implemented by `descr_a=['T' int2str(L1) ',' int2str(M1)]`.
- Lines 84-86: Get the Spinach state; implemented by `rho=state(spin_system,{descr_a,descr_b}, {spin_a, spin_b })`.
- Lines 88-89: Get and normalise the two-spin state; implemented by `tss=kron(IST{n},IST{k}); tss=tss/norm(tss,'fro')`.
- Lines 91-92: Normalise the Spinach state; implemented by `switch spin_system.bas.formalism`.
- Lines 96-97: Hilbert space; implemented by `rho=rho/norm(rho,'fro')`.
- Lines 101-102: Liouville space; implemented by `rho=rho/norm(full(rho),2)`.
- Lines 106-107: Complain and bomb out; implemented by `error('unsupported formalism.')`.
- Lines 111-112: Singlet state population; implemented by `S=S+trace(S0'*tss'*S0)*rho`.
- Lines 114-115: Triplet state populations; implemented by `if nargout>1`.
- Lines 121-122: Quintet state populations; implemented by `if nargout>2`.
- Lines 130-131: Triplet state coherences; implemented by `if nargout>3`.
- Lines 138-139: Quintet state coherences; implemented by `if nargout>4`.

### Control flow inferred from the code

- Line 50: conditional branch on `~exist('options','var'), options.dephasing=0; end`.
- Line 71: `for` loop over `n=1:numel(IST)`.
- Line 72: `for` loop over `k=1:numel(IST)`.
- Line 78: conditional branch on `((M1~=0)||(M2~=0))&&(options.dephasing==1), continue; end`.
- Line 92: dispatches on `spin_system.bas.formalism`; cases `'zeeman-hilb'`, `{'zeeman-liouv','sphten-liouv'}`.
- Line 115: conditional branch on `nargout>1`.
- Line 122: conditional branch on `nargout>2`.
- Line 131: conditional branch on `nargout>3`.
- Line 139: conditional branch on `nargout>4`.

### Key state/data transformations

- Lines 57: computes `alp` using `alp=[1; 0; 0]; bet=[0; 1; 0]; gam=[0; 0; 1]`.
- Lines 58: computes `S0` using `S0=(1/sqrt(3))*(kron(alp,gam)-kron(bet,bet)+kron(gam,alp))`.
- Lines 59: computes `Tp` using `Tp=(1/sqrt(2))*(kron(alp,bet)-kron(bet,alp))`.
- Lines 60: computes `T0` using `T0=(1/sqrt(2))*(kron(alp,gam)-kron(gam,alp))`.
- Lines 61: computes `Tm` using `Tm=(1/sqrt(2))*(kron(bet,gam)-kron(gam,bet))`.
- Lines 62: computes `Qpp` using `Qpp=kron(alp,alp)`.
- Lines 63: computes `Qp` using `Qp=(1/sqrt(2))*(kron(alp,bet)+kron(bet,alp))`.
- Lines 64: computes `Q0` using `Q0=(1/sqrt(6))*(kron(alp,gam)+2*kron(bet,bet)+kron(gam,alp))`.
- Lines 65: computes `Qm` using `Qm=(1/sqrt(2))*(kron(bet,gam)+kron(gam,bet))`.
- Lines 66: computes `Qmm` using `Qmm=kron(gam,gam)`.
- Lines 69: computes `S` using `S=0; T={0,0,0}; Q={0,0,0,0,0}; Tc={0,0,0,0}`.
- Lines 70: computes `Qc` using `Qc={0,0,0,0,0,0,0,0}; IST=irr_sph_ten(3)`.
- Lines 75: computes `[L1,M1]` using `[L1,M1]=lin2lm(n-1); [L2,M2]=lin2lm(k-1)`.
- Lines 81: computes `descr_a` using `descr_a=['T' int2str(L1) ',' int2str(M1)]`.
- Lines 82: computes `descr_b` using `descr_b=['T' int2str(L2) ',' int2str(M2)]`.
- Lines 85-86: computes `rho` using `rho=state(spin_system,{descr_a,descr_b}, {spin_a, spin_b })`.
- Lines 89: computes `tss` using `tss=kron(IST{n},IST{k}); tss=tss/norm(tss,'fro')`.
- Lines 116: computes `T{1}` using `T{1}=T{1}+trace(Tp' *tss'*Tp )*rho`.

### Local helper functions

- Line 157: `grumble()` — `function grumble(spin_system,spin_a,spin_b,options)`.
  - Representative operation: `if (~isnumeric(spin_a))||(~isnumeric(spin_b))|| (~isscalar(spin_a))||(~isscalar(spin_b))|| (mod(spin_a,1)~=0)||(mod(spin_b,1)~=0)|| (spin_a<1)||(spin_b<1)||(spin_a==spin…`.
  - Representative operation: `(~isscalar(spin_a))||(~isscalar(spin_b))|| (mod(spin_a,1)~=0)||(mod(spin_b,1)~=0)|| (spin_a<1)||(spin_b<1)||(spin_a==spin_b)`.

## Parameters / inputs

- spin_a -the number of the first spin
- spin_b -the number of the second spin
- options.dephasing -set to 1 to eliminate the states that are
- not stationary under the Zeeman Hamilto-
- nian of two inequivalent spins, i.e. to
- keep only zero-projection products; the
- default is to keep everything

## Outputs

- S -singet state density matrix (Hilbert space)
- or state vector (Liouville space)
- T -triplet state density matrices (Hilbert space)
- or state vectors (Liouville space), ordered in
- a cell array as {T+,T0,T-}
- Q -quintet state density matrices (Hilbert space)
- or state vectors (Liouville space), ordered in
- a cell array as {Q++,Q+,Q0,Q-,Q--}
- Tc -coherences between triplet states:
- {T0 -> T-, T+ -> T0, T--> T0, T0 -> T+}
- Qc -coherences between quintet states:
- {Q--> Q--, Q0 -> Q-, Q+ -> Q0, Q++ -> Q+, ...
- Q---> Q-, Q--> Q0, Q0 -> Q+, Q+ -> Q++ }
- WARNING: the states above are NOT irreducible spherical tensors -
- Bargon just kroneckered up some Zeeman states and gave
- them what looked to him like reasonable labels.

## Implementation structure

- All possible states of a spin-1 pair, classified by the total spin
- into singlet, triplet, and quartet. Syntax:
- [S,T,Q,Tc,Qc]=deut_pair(spin_system,spin_a,spin_b,options)
- spin_a -the number of the first spin
- spin_b -the number of the second spin
- options.dephasing -set to 1 to eliminate the states that are
- not stationary under the Zeeman Hamilto-
- nian of two inequivalent spins, i.e. to
- keep only zero-projection products; the
- default is to keep everything
- S -singet state density matrix (Hilbert space)
- or state vector (Liouville space)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `irr_sph_ten()`, `lin2lm()`, `int2str()`, `state()`, `isscalar()`, `isfield()`, `ismember()`.
