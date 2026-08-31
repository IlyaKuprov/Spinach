# kernel/optimcon/aux_mat.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/aux_mat.m`
- Signature: `[auxm_l,auxm_r]=aux_mat(drifts,controls,cc_comm_idx,...`
- Total lines: 210

## Purpose

Builds auxiliary matrices for the calculation of the directional derivatives of the trapezium product quadrature propagator: expm(-1i*((HL+HR)/2+(1i*dt/12)*[HL,HR])*dt) with respect to control coefficients in the evolution generators HL and HR on the left and the right edge of the interval. The de- rivatives are calculated using Eq 16 of Goodwin and Kuprov:

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 63-64: If no second index, define as 0; implemented by `if nargin<9, j=0; end`.
- Lines 66-67: Check consistency; implemented by `grumble(drifts,controls,cc_comm_idx,cc_comm,dt,cL,cR,k,j)`.
- Lines 69-70: Build left and right generators; implemented by `G{1}=drifts{1}; G{2}=drifts{2}`.
- Lines 76-77: Build the overall generator and a blank matrix; implemented by `G=(G{1}+G{2})/2+1i*dt*(1/12)*(G{1}*G{2}-G{2}*G{1})`.
- Lines 80-82: Start right directional derivative of the generator; implemented by `Hk_dir_R=(1/2)*controls{k}+1i*dt*(1/12)*(drifts{1}*controls{k}- controls{k}*drifts{1})`.
- Lines 84-86: Start left directional derivative of the generator; implemented by `Hk_dir_L=(1/2)*controls{k}+1i*dt*(1/12)*(controls{k}*drifts{2}- drifts{2}*controls{k})`.
- Lines 88-89: Contributions from other controls; implemented by `for n=1:numel(controls)`.
- Lines 91-92: If the commutator is non-zero; implemented by `if ~cc_comm_idx(n,k)`.
- Lines 94-95: Add right side contribution; implemented by `Hk_dir_R=Hk_dir_R+1i*dt*(1/12)*cL(n)*cc_comm{n,k}`.
- Lines 97-98: Add left side contribution; implemented by `Hk_dir_L=Hk_dir_L+1i*dt*(1/12)*cR(n)*cc_comm{k,n}`.
- Lines 104-105: For 3x3 block matrices; implemented by `if j~=0`.
- Lines 107-109: Start right directional derivative of the generator; implemented by `Hj_dir_R=(1/2)*controls{j}+1i*dt*(1/12)*(drifts{1}*controls{j}- controls{j}*drifts{1})`.
- Lines 111-113: Start left directional derivative of the generator; implemented by `Hj_dir_L=(1/2)*controls{j}+1i*dt*(1/12)*(controls{j}*drifts{2}- drifts{2}*controls{j})`.
- Lines 118-119: If the commutator is non-zero; implemented by `if ~cc_comm_idx(n,j)`.
- Lines 121-122: Add right side contribution; implemented by `Hj_dir_R=Hj_dir_R+1i*dt*(1/12)*cL(n)*cc_comm{n,j}`.
- Lines 124-125: Add left side contribution; implemented by `Hj_dir_L=Hj_dir_L+1i*dt*(1/12)*cR(n)*cc_comm{j,n}`.
- Lines 133-134: Build auxiliary matrices; implemented by `if j==0`.
- Lines 136-137: 2x2 block matrix; implemented by `auxm_l=[G, Hk_dir_L; B, G]`.

### Control flow inferred from the code

- Line 64: conditional branch on `nargin<9, j=0; end`.
- Line 71: `for` loop over `n=1:numel(controls)`.
- Line 89: `for` loop over `n=1:numel(controls)`.
- Line 92: conditional branch on `~cc_comm_idx(n,k)`.
- Line 105: conditional branch on `j~=0`.
- Line 116: `for` loop over `n=1:numel(controls)`.
- Line 119: conditional branch on `~cc_comm_idx(n,j)`.
- Line 134: conditional branch on `j==0`.

### Key state/data transformations

- Lines 70: computes `G{1}` using `G{1}=drifts{1}; G{2}=drifts{2}`.
- Lines 73: computes `G{2}` using `G{2}=G{2}+cR(n)*controls{n}`.
- Lines 77: computes `G` using `G=(G{1}+G{2})/2+1i*dt*(1/12)*(G{1}*G{2}-G{2}*G{1})`.
- Lines 78: computes `B` using `B=spalloc(size(G,1),size(G,2),0)`.
- Lines 81-82: computes `Hk_dir_R` using `Hk_dir_R=(1/2)*controls{k}+1i*dt*(1/12)*(drifts{1}*controls{k}- controls{k}*drifts{1})`.
- Lines 85-86: computes `Hk_dir_L` using `Hk_dir_L=(1/2)*controls{k}+1i*dt*(1/12)*(controls{k}*drifts{2}- drifts{2}*controls{k})`.
- Lines 108-109: computes `Hj_dir_R` using `Hj_dir_R=(1/2)*controls{j}+1i*dt*(1/12)*(drifts{1}*controls{j}- controls{j}*drifts{1})`.
- Lines 112-113: computes `Hj_dir_L` using `Hj_dir_L=(1/2)*controls{j}+1i*dt*(1/12)*(controls{j}*drifts{2}- drifts{2}*controls{j})`.
- Lines 137: computes `auxm_l` using `auxm_l=[G, Hk_dir_L; B, G]`.
- Lines 138: computes `auxm_r` using `auxm_r=[G, Hk_dir_R; B, G]`.

### Local helper functions

- Line 151: `grumble()` — `function grumble(drifts,controls,cc_comm_idx,cc_comm,dt,cL,cR,k,j)`.
  - Representative operation: `if ~iscell(drifts)`.
  - Representative operation: `error('drifts must be a cell array of matrices.')`.

## Syntax

```matlab
[auxm_l,auxm_r]=aux_mat(drifts,controls,cc_comm_idx,...
cc_comm,dt,cL,cR,k,j)
```

## Parameters / inputs

- drifts -a cell array of two matrices containing drift
- generators at the left (first element) and the
- right (second element) edge of the interval
- controls -a cell array of K control generators
- cc_comm_idx -a KxK matrix of logicals indicating non-zero
- commutation of controls
- cc_comm -a KxK cell array control commutation relarions
- dt -interval duration, seconds
- cL -control generator coefficients at the left
- edge of the interval
- cR -control generator coefficients at the right
- edge of the interval
- k -the index of the generator inside controls
- array that the differentiation refers to
- j -(optional) the index of the 2nd generator
- inside controls array that the differentiation
- refers to. Required for 3x3 block auxiliary
- matrices

## Outputs

- auxm_l -auxilary matrix for the derivative of the
- interval propagator with respect to cL
- auxm_r -auxilary matrix for the derivative of the
- interval propagator with respect to cR

## Implementation structure

- Builds auxiliary matrices for the calculation of the directional
- derivatives of the trapezium product quadrature propagator:
- expm(-1i*((HL+HR)/2+(1i*dt/12)*[HL,HR])*dt)
- with respect to control coefficients in the evolution generators
- HL and HR on the left and the right edge of the interval. The de-
- rivatives are calculated using Eq 16 of Goodwin and Kuprov:
- [auxm_l,auxm_r]=aux_mat(drifts,controls,cc_comm_idx,...
- cc_comm,dt,cL,cR,k,j)
- drifts -a cell array of two matrices containing drift
- generators at the left (first element) and the
- right (second element) edge of the interval
- controls -a cell array of K control generators

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spalloc()`, `cc_comm_idx()`, `iscell()`, `islogical()`, `isequal()`, `isscalar()`.
