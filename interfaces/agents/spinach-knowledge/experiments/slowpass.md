# experiments/slowpass.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/slowpass.m`
- Signature: `spectrum=slowpass(spin_system,parameters,H,R,K)`
- Total lines: 199

## Purpose

Slow passage detection -calculates spectrum values at the user- specified frequency positions using the Fourier transform of the Liouville -von Neumann equation. The biggest advantage over the fid+fft style detection is easy parallelization and the possibi- lity of getting spectrum values at specific frequencies without recalculating the entire free induction decay. Syntax: spectrum=slowpass(spin_system,parameters,H,

## Physical / mathematical content

- This file belongs to the `experiments` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters.sweep vector with two elements giving
- the spectrum frequency extents
- in Hz
- parameters.npoints number of points in the spectrum
- parameters.rho0 initial state
- parameters.coil detection state
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- spectrum -the spectrum of the system with the specified
- starting state detected on the specified coil
- state within the frequency interval requested
- Note: relaxation must be present in the system dynamics, or the
- matrix inversion operation would fail to converge. The re-
- laxation matrix R must *not* be thermalized.

## Implementation structure

- Slow passage detection -calculates spectrum values at the user-
- specified frequency positions using the Fourier transform of the
- Liouville -von Neumann equation. The biggest advantage over the
- fid+fft style detection is easy parallelization and the possibi-
- lity of getting spectrum values at specific frequencies without
- recalculating the entire free induction decay. Syntax:
- spectrum=slowpass(spin_system,parameters,H,R,K)
- parameters.sweep vector with two elements giving
- the spectrum frequency extents
- in Hz
- parameters.npoints number of points in the spectrum
- parameters.rho0 initial state

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `sim2liouv()`, `reduce()`, `speye()`, `ismember()`, `report()`, `gpuArray()`, `dot()`, `freq_grid()`, `spectrum()`, `gather()`, `ilu()`, `nnz()`, `num2str()`, `gmres()`, `ismatrix()`.
