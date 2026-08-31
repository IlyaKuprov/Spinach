# kernel/pulses/wave_basis.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/wave_basis.m`
- Signature: `basis_waves=wave_basis(basis_type,n_func,n_points)`
- Total lines: 105

## Purpose

Common basis sets for the expansion of pulse waveforms. Returns the wave- form basis functions as columns of a matrix. Syntax: basis_waves=wave_basis(basis_type,n_functions,n_steps)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 39-40: Check consistency; implemented by `grumble(basis_type,n_func,n_points)`.
- Lines 42-43: Preallocate the array; implemented by `basis_waves=zeros(n_func,n_points)`.
- Lines 45-46: Generate the functions; implemented by `switch basis_type`.
- Lines 50-51: Fill the array; implemented by `parfor n=1:n_func`.
- Lines 72-73: Complain and bomb out; implemented by `error('unrecognized basis function type.')`.
- Lines 77-78: Orthogonalize the array; implemented by `basis_waves=orth(basis_waves')`.
- Lines 80-81: Check outgoing dimension; implemented by `if size(basis_waves,2)~=n_func`.

### Control flow inferred from the code

- Line 46: dispatches on `basis_type`; cases `'sine_waves'`, `'cosine_waves'`, `'legendre'`.
- Line 51: `parfor` loop over `n=1:n_func`.
- Line 58: `parfor` loop over `n=1:n_func`.
- Line 65: `parfor` loop over `n=1:n_func`.
- Line 81: conditional branch on `size(basis_waves,2)~=n_func`.

### Key state/data transformations

- Lines 43: computes `basis_waves` using `basis_waves=zeros(n_func,n_points)`.
- Lines 52: computes `basis_waves(n,:)` using `basis_waves(n,:)=sin(n*linspace(-pi,pi,n_points))`.

### Local helper functions

- Line 88: `grumble()` — `function grumble(basis_type,n_func,n_points)`.
  - Representative operation: `if ~ischar(basis_type)`.
  - Representative operation: `error('basis_type parameter must be a character string.')`.

## Parameters / inputs

- basis_type -may be set to 'sine_waves', 'cosine_waves',
- and 'legendre'. The sine and the cosine op-
- tions return the corresponding functions in
- the [-pi,pi] interval, legendre option re-
- turns legendre polynomials in the [-1,1] in-
- terval.
- n_func -the number of functions to return (integer
- frequencies starting from zero on the case
- of cosines, integer frequencies starting
- from 1 inthe case of sines, legendre poly-
- nomial ranks in the case of legendre func-
- tion basis set.
- n_points -number of discretization points.

## Outputs

- basis_waves -a matrix with the basis waves in columns
- Note: because the resulting waveforms are discretised, they are not pre-
- cisely orthogonal under the standard scalar multiplication. An ex-
- tra orthogonalisation step is therefore applied to make them ortho-
- gonal as vectors. As a result, some functions may be upside-down.

## Implementation structure

- Common basis sets for the expansion of pulse waveforms. Returns the wave-
- form basis functions as columns of a matrix. Syntax:
- basis_waves=wave_basis(basis_type,n_functions,n_steps)
- basis_type -may be set to 'sine_waves', 'cosine_waves',
- and 'legendre'. The sine and the cosine op-
- tions return the corresponding functions in
- the [-pi,pi] interval, legendre option re-
- turns legendre polynomials in the [-1,1] in-
- terval.
- n_func -the number of functions to return (integer
- frequencies starting from zero on the case
- of cosines, integer frequencies starting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `basis_waves()`, `legendreP()`, `orth()`, `ischar()`.
