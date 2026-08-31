# kernel/utilities/rlx_scalar.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/rlx_scalar.m`
- Signature: `R=rlx_scalar(spin_system,H0,H1,tau_c_array)`
- Total lines: 90

## Purpose

Scalar relaxation superoperator using Redfield theory. Syntax: R=rlx_scalar(spin_system,H0,H1,tau_c_array)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 33-34: Check consistency; implemented by `grumble(H0,H1,tau_c_array)`.
- Lines 36-37: Get the ouput started; implemented by `R=sparse(0)`.
- Lines 39-40: Loop over correlation function components; implemented by `for n=1:numel(tau_c_array)`.
- Lines 42-43: Extract component weight and correlation time; implemented by `weight=tau_c_array{n}(1); tau_c=tau_c_array{n}(2)`.
- Lines 45-46: Ignore insignificant components; implemented by `if (weight~=0)&&(tau_c~=0)`.
- Lines 48-49: Set the upper integration limit according to the accuracy goal; implemented by `upper_limit=2*tau_c*log(1/spin_system.tols.rlx_integration)`.
- Lines 51-52: Remove inconsequential non-zeroes from a copy of H0; implemented by `H0c=clean_up(spin_system,H0,1e-2/upper_limit)`.
- Lines 54-55: Take the integral using the auxiliary matrix exponential technique; implemented by `R=R-weight*H1*expmint(spin_system,H0c,H1',H0c+(1i/tau_c)*speye(size(H0c)),upper_limit)`.

### Control flow inferred from the code

- Line 40: `for` loop over `n=1:numel(tau_c_array)`.
- Line 46: conditional branch on `(weight~=0)&&(tau_c~=0)`.

### Key state/data transformations

- Lines 37: computes `R` using `R=sparse(0)`.
- Lines 43: computes `weight` using `weight=tau_c_array{n}(1); tau_c=tau_c_array{n}(2)`.
- Lines 49: computes `upper_limit` using `upper_limit=2*tau_c*log(1/spin_system.tols.rlx_integration)`.
- Lines 52: computes `H0c` using `H0c=clean_up(spin_system,H0,1e-2/upper_limit)`.

### Local helper functions

- Line 64: `grumble()` — `function grumble(H0,H1,tau_c_array)`.
  - Representative operation: `if (~isnumeric(H0))||(~isnumeric(H1))||(~ismatrix(H0))|| (~ismatrix(H1))||(~ishermitian(H0))||(~ishermitian(H1))`.
  - Representative operation: `(~ismatrix(H1))||(~ishermitian(H0))||(~ishermitian(H1))`.

## Parameters / inputs

- H0 -background Hamiltonian
- H1 -the stochastically modulated interaction operator
- multiplied by its root mean square modulation depth
- tau_c_array -a cell array of the following format:
- {[weight_a,tau_a],[weight_b,tau_b],...}
- giving weights of the exponential components
- of the correlation function and the associa-
- ted correlation times, e.g. {[1.0,1e-12]}

## Outputs

- R -relaxation superoperator, a negative definite matrix
- Note: if H1(t) has a non-zero time or ensemble average value,
- that average must be subtracted out and placed into H0

## Implementation structure

- Scalar relaxation superoperator using Redfield theory. Syntax:
- R=rlx_scalar(spin_system,H0,H1,tau_c_array)
- H0 -background Hamiltonian
- H1 -the stochastically modulated interaction operator
- multiplied by its root mean square modulation depth
- tau_c_array -a cell array of the following format:
- {[weight_a,tau_a],[weight_b,tau_b],...}
- giving weights of the exponential components
- of the correlation function and the associa-
- ted correlation times, e.g. {[1.0,1e-12]}
- R -relaxation superoperator, a negative definite matrix
- Note: if H1(t) has a non-zero time or ensemble average value,

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `clean_up()`, `expmint()`, `speye()`, `ismatrix()`, `ishermitian()`, `all()`, `iscell()`.
