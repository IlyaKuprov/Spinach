# experiments/relaxan.m

- Signature: `[r1,r2,t1,t2,R]=relaxan(spin_system,euler_angles)`

## Purpose

Automated relaxation theory analysis. Prints longitudinal and transverse relaxation rates and times for all spins in the system. Syntax: [r1,r2,t1,t2,R]=relaxan(spin_system,euler_angles)

## Physical / mathematical content

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Parameters / inputs

- euler_angles -optional euler angles for situations
- when relaxation properties are orien-
- tation-dependent

## Outputs

- r1 -a vector of longitudinal relaxation rates
- for each spin
- r2 -a vector of transverse relaxation rates
- for each spin
- t1 -a vector of longitudinal relaxation times
- for each spin
- t2 -a vector of transverse relaxation times
- for each spin
- R -complete relaxation superoperator
- Note: dynamic frequency shifts are dropped.

## Implementation structure

- Automated relaxation theory analysis. Prints longitudinal
- and transverse relaxation rates and times for all spins in
- the system. Syntax:
- [r1,r2,t1,t2,R]=relaxan(spin_system,euler_angles)
- euler_angles -optional euler angles for situations
- when relaxation properties are orien-
- tation-dependent
- r1 -a vector of longitudinal relaxation rates
- for each spin
- r2 -a vector of transverse relaxation rates
- t1 -a vector of longitudinal relaxation times
- t2 -a vector of transverse relaxation times
