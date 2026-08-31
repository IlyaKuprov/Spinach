# kernel/utilities/corrfun.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/corrfun.m`
- Signature: `[weights,rates,states]=corrfun(spin_system,n,k,m,p,q)`
- Total lines: 162

## Purpose

Wigner matrix element correlation function under isotropic, axial, and rhombic rotational diffusion. Syntax: [weights,rates,states]=corrfun(spin_system,n,k,m,p,q)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 56-57: Check consistency; implemented by `grumble(spin_system,n,k,m,p,q)`.
- Lines 59-60: Preallocate outputs; implemented by `weights=cell(1,numel(spin_system.chem.parts))`.
- Lines 64-65: Index basis states for different chemical species; implemented by `for s=1:numel(spin_system.chem.parts)`.
- Lines 69-70: Loop over chemical species; implemented by `for s=1:numel(spin_system.chem.parts)`.
- Lines 72-73: Select rotational diffusion model; implemented by `switch numel(spin_system.rlx.tau_c{s})`.
- Lines 77-78: Assume second rank correlation time is supplied; implemented by `D=1/(6*spin_system.rlx.tau_c{s})`.
- Lines 80-81: Use isotropic rotational diffusion model; implemented by `weights{s}=(1/(2*n+1))*krondelta(k,p)*krondelta(m,q)`.
- Lines 86-87: Assume second rank correlation time is supplied; implemented by `D_ax=1/(6*spin_system.rlx.tau_c{s}(1))`.
- Lines 90-91: Use axial rotational diffusion model; implemented by `weights{s}=(1/(2*n+1))*krondelta(k,p)*krondelta(m,q)`.
- Lines 96-97: Use anisotropic rotational diffusion model, only L=2 permitted; implemented by `Dxx=1/(6*spin_system.rlx.tau_c{s}(1))`.
- Lines 101-104: Refuse to process degenerate cases and ranks other than 2; implemented by `if (abs(Dxx-Dyy)<1e-6*mean([Dxx Dyy Dzz]))|| (abs(Dyy-Dzz)<1e-6*mean([Dxx Dyy Dzz]))|| (abs(Dzz-Dxx)<1e-6*mean([Dxx Dyy Dzz]))`.
- Lines 111-112: Compute decay rates; implemented by `delta=sqrt(Dxx^2+Dyy^2+Dzz^2-Dxx*Dyy-Dxx*Dzz-Dyy*Dzz)`.
- Lines 119-120: Compute coefficients; implemented by `lambda_p=sqrt(2/3)*(Dxx+Dyy-2*Dzz+2*delta)/(Dxx-Dyy)`.
- Lines 127-128: Compute weights; implemented by `for j=1:5, weights{s}(j)=(1/5)*krondelta(k,p)*h(j,m)*h(j,q); end`.

### Control flow inferred from the code

- Line 65: `for` loop over `s=1:numel(spin_system.chem.parts)`.
- Line 70: `for` loop over `s=1:numel(spin_system.chem.parts)`.
- Line 73: dispatches on `numel(spin_system.rlx.tau_c{s})`; cases `1`, `2`, `3`.
- Line 102: conditional branch on `(abs(Dxx-Dyy)<1e-6*mean([Dxx Dyy Dzz]))||`.
- Line 107: conditional branch on `n~=2`.
- Line 128: `for` loop over `j=1:5, weights{s}(j)=(1/5)*krondelta(k,p)*h(j,m)*h(j,q); end`.

### Key state/data transformations

- Lines 60: computes `weights` using `weights=cell(1,numel(spin_system.chem.parts))`.
- Lines 61: computes `rates` using `rates=cell(1,numel(spin_system.chem.parts))`.
- Lines 62: computes `states` using `states=cell(1,numel(spin_system.chem.parts))`.
- Lines 66: computes `states{s}` using `states{s}=(sum(spin_system.bas.basis(:,spin_system.chem.parts{s}),2)>0)`.
- Lines 78: computes `D` using `D=1/(6*spin_system.rlx.tau_c{s})`.
- Lines 81: computes `weights{s}` using `weights{s}=(1/(2*n+1))*krondelta(k,p)*krondelta(m,q)`.
- Lines 82: computes `rates{s}` using `rates{s}=-n*(n+1)*D`.
- Lines 87: computes `D_ax` using `D_ax=1/(6*spin_system.rlx.tau_c{s}(1))`.
- Lines 88: computes `D_eq` using `D_eq=1/(6*spin_system.rlx.tau_c{s}(2))`.
- Lines 97: computes `Dxx` using `Dxx=1/(6*spin_system.rlx.tau_c{s}(1))`.
- Lines 98: computes `Dyy` using `Dyy=1/(6*spin_system.rlx.tau_c{s}(2))`.
- Lines 99: computes `Dzz` using `Dzz=1/(6*spin_system.rlx.tau_c{s}(3))`.
- Lines 108: computes `error('rhombic rotational diffusion only implemented for n` using `error('rhombic rotational diffusion only implemented for n=2.')`.
- Lines 112: computes `delta` using `delta=sqrt(Dxx^2+Dyy^2+Dzz^2-Dxx*Dyy-Dxx*Dzz-Dyy*Dzz)`.
- Lines 113: computes `rates{s}(1)` using `rates{s}(1)=-(4*Dxx+Dyy+Dzz)`.
- Lines 114: computes `rates{s}(2)` using `rates{s}(2)=-(Dxx+4*Dyy+Dzz)`.
- Lines 115: computes `rates{s}(3)` using `rates{s}(3)=-(Dxx+Dyy+4*Dzz)`.
- Lines 116: computes `rates{s}(4)` using `rates{s}(4)=-(2*Dxx+2*Dyy+2*Dzz-2*delta)`.

### Local helper functions

- Line 137: `grumble()` — `function grumble(spin_system,n,k,m,p,q)`.
  - Representative operation: `if ~strcmp(spin_system.bas.formalism,'sphten-liouv')`.
  - Representative operation: `error('Redfield relaxation theory is only available for sphten-liouv formalism.')`.

## Parameters / inputs

- spin_system -the output of [[create.m]] to which ro-
- tational correlation time should have
- been supplied. For a single correlation
- time, the isotropic rotational diffusi-
- on model is used; a vector with two
- correlation times is assumed to be cor-
- relation times for rotation around and
- perpendicularly to the main axis res-
- pectively); a vector with three corre-
- lation times is assumed to be the cor-
- relation times for the rotation around
- the XX, YY and ZZ direction respecti-
- vely of the rotational diffusion tensor.
- n,k,m,p,q -the five indices found in the ensemble-
- averaged Wigner function product:
- <D{n}{k,m}(0)*D{n}{p,q}(t)'>

## Outputs

- weights -a cell array (one element for each che-
- mical species) of vectors listing the
- weights of the exponential components
- of the decays
- rates -a cell array (one element for each che-
- mical species) of vectors listing the
- decay rates (negative numbers) of the
- exponential components of the decays
- states -a cell array (one element for each che-
- mical species) of logical vectors indi-
- cating which states in the basis set
- belong to which chemical species
- Note: Wigner function indices are sorted in descending order, that is,
- k=[1 2 3 4 5] in the input represents [2 1 0 -1 -2] for n=2.
- Note: second rank rotational correlation times (as per Spinach input)
- will be updated automatically if other ranks are specified.

## Implementation structure

- Wigner matrix element correlation function under isotropic, axial,
- and rhombic rotational diffusion. Syntax:
- [weights,rates,states]=corrfun(spin_system,n,k,m,p,q)
- spin_system -the output of [[create.m]] to which ro-
- tational correlation time should have
- been supplied. For a single correlation
- time, the isotropic rotational diffusi-
- on model is used; a vector with two
- correlation times is assumed to be cor-
- relation times for rotation around and
- perpendicularly to the main axis res-
- pectively); a vector with three corre-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `krondelta()`, `strcmp()`, `any()`.
