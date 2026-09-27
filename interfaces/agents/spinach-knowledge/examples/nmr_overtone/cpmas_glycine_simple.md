# examples/nmr_overtone/cpmas_glycine_simple.m

- Signature: `cpmas_glycine_simple()`

## Purpose

Simulates a 14N-overtone/proton cross-polarisation spectrum for glycine under MAS. The source estimates a calculation time of hours and credits Ilya Kuprov, M. Carravetta, and M. Concistre.

## Physical / mathematical content

The source attributes the glycine quadrupolar tensor data to O'Dell and Ratcliffe ([DOI](http://dx.doi.org/10.1016/j.cplett.2011.08.030)). It specifies 14N and 1H at 14.1 T, 14N quadrupolar parameters of 1.18 MHz and η=0.53, and a 1H shift of 32.4. Damping relaxation is used with diagonal retention, zero equilibrium, and rate 300. The basis is sphten-liouv without approximation.

## Numerical / algorithmic content

The source disables Krylov and trajectory-level options and uses the 6400-point spherical powder grid `rep_2ang_6400pts_sph`, rank 7, and a MAS rate of −19.840 kHz. The spectrum spans [44, 52] kHz with 256 acquired and zero-filled points; the RF powers are 55.0 and 35.1 kHz, the RF frequency is 48 kHz, and the RF duration is 100 μs. The spectrum is computed with `singlerot` and `@overtone_cp`.

## Implementation structure

The function builds the system and basis, sets the MAS and cross-polarisation operators, runs the single spectrum simulation, and plots its real part.
