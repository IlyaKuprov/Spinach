# examples/nmr_overtone/cpmas_glycine_accum.m

- Signature: `cpmas_glycine_accum()`

## Purpose

Shows a 14N-overtone/proton cross-polarisation accumulation profile in glycine under MAS as the RF contact duration is varied. The source estimates minutes of computation and credits Ilya Kuprov, M. Carravetta, and M. Concistre.

## Physical / mathematical content

The source attributes the glycine quadrupolar tensor data to O'Dell and Ratcliffe ([DOI](http://dx.doi.org/10.1016/j.cplett.2011.08.030)). It specifies a 14.10220742 T field, 14N and 1H spins, a 14N quadrupolar tensor from `eeqq2nqi(1.18e6,0.53,1,[0 0 0])`, and a 1H shift of 32.4. The two spin coordinates are separated by 1 Å. Relaxation is damped, with diagonal retention, zero equilibrium, and rate 300.

## Numerical / algorithmic content

The simulation uses the sphten-liouv basis without approximation, disables Krylov and trajectory-level options, and uses a 200-point octahedral powder grid (`rep_2ang_200pts_oct`) with rank 7. The spectrum spans [44, 52] kHz with 256 acquired and zero-filled points. The RF contact duration is stepped from 10 to 100 μs in ten increments; each step calls `singlerot` with `@overtone_cp`.

## Implementation structure

The code constructs the spin system and basis, defines the MAS axis and CP preparation/detection operators, then loops over contact durations, simulates and plots each spectrum in a panel.
