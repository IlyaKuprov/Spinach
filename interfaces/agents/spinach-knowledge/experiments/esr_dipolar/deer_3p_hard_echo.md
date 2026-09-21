# experiments/esr_dipolar/deer_3p_hard_echo.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/esr_dipolar/deer_3p_hard_echo.m`
- Signature: `echo=deer_3p_hard_echo(spin_system,parameters,H,R,K)`
- Total lines: 137

## Purpose

Samples the spin echo in the three-pulse DEER experiment to determine its precise location -in simulations, DEER spin echoes can be narrow and easy to miss. For high-spin electrons the echo may not be in the expected place. Ideal hard pulses are used, each pulse only hits its specific electron or transition, depending on the pulse operator supp- lied. Syntax: echo=deer_3p_hard_echo(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Dipolar ESR experiment implementations. The pulse logic resolves dipolar couplings by echo modulation, with selective excitation and time-domain accumulation.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters.ex_prob -probe pulse operator
- parameters.ex_pump -pump pulse operator
- parameters.ta -time between first and third pulse, s
- parameters.tb -time between first and second pulse, s
- parameters.tc -time interval to sample around the
- echo, s
- parameters.rho0 -initial condition
- parameters.nsteps -number of sampling steps in the tc
- interval
- parameters.coil -detection state
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- echo -the signal detected during the parameters.tc interval
- with parameters.nsteps points in it

## Implementation structure

- Samples the spin echo in the three-pulse DEER experiment to determine
- its precise location -in simulations, DEER spin echoes can be narrow
- and easy to miss. For high-spin electrons the echo may not be in the
- expected place. Ideal hard pulses are used, each pulse only hits its
- specific electron or transition, depending on the pulse operator supp-
- lied. Syntax:
- echo=deer_3p_hard_echo(spin_system,parameters,H,R,K)
- parameters.ex_prob -probe pulse operator
- parameters.ex_pump -pump pulse operator
- parameters.ta -time between first and third pulse, s
- parameters.tb -time between first and second pulse, s
- parameters.tc -time interval to sample around the

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `step()`, `evolution()`, `ismatrix()`, `all()`, `isfield()`.
