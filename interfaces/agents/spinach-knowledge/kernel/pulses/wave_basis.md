# kernel/pulses/wave_basis.m

- Signature: `basis_waves=wave_basis(basis_type,n_func,n_points)`

## Purpose

Common basis sets for the expansion of pulse waveforms. Returns the wave- form basis functions as columns of a matrix. Syntax: basis_waves=wave_basis(basis_type,n_functions,n_steps)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

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
